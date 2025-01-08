#!/bin/bash
#
# run_atac.sh
#
# 自动遍历输入文件夹，识别 R1/R2 FASTQ 文件并执行 ATAC-seq 分析流程

set -e  # 若任意命令出错，则脚本退出

############################################
#                参数设置部分
############################################

# 如果没有输入任何参数，或第一个参数不是文件夹，则退出
if [ $# -lt 1 ] || [ ! -d "\$1" ]; then
  echo "用法: \$0 <fastq_folder>"
  echo "其中 <fastq_folder> 是包含FASTQ文件的文件夹"
  exit 1
fi

# 1) 输入FASTQ所在目录
fastq_dir="\$1"

# 2) 参考基因组的Bowtie2索引
genome="/home/maolp/mao/Ref/AllnewstarRef/Homo/HG19/HG19BT/hg19"

# 3) 输出目录
mkdir -p trimmed bams stats beds peaks bws

############################################
#  辅助函数：从R1文件名推断R2文件名
############################################
replace_R1_with_R2() {
  # 处理下列几种R1后缀：
  # 1) .R1.fastq.gz
  # 2) .R1.fq.gz
  # 3) _1.fastq.gz
  # 4) _1.fq.gz
  # 转换为相应的 R2 后缀
  local r1_file="\$1"
  local r2_file

  # 用 sed 做简单的模式替换：
  # 将所有包含 R1 或者 _1 的部分替换为 R2 或者 _2
  # 例如： xxx.R1.fastq.gz -> xxx.R2.fastq.gz
  #        xxx_1.fq.gz    -> xxx_2.fq.gz
  # 注意这里用正则匹配多种可能:
  # ([._]) 代表匹配 '.' 或 '_'
  # (R|_)1 代表匹配 R1 或 _1
  # (\.fastq\.gz|\.fq\.gz) 代表匹配结尾的 .fastq.gz 或 .fq.gz
  r2_file=$(echo "$r1_file" | sed -E 's/([._])(R|_)1(\.fastq\.gz|\.fq\.gz)$/\1\22\3/')

  echo "$r2_file"
}

############################################
#  收集所有 R1 文件，并自动推断 R2
############################################
# 通配符: *{.R1.fastq.gz,.R1.fq.gz,_1.fastq.gz,_1.fq.gz}
# 意思是：匹配目录下以 .R1.fastq.gz、.R1.fq.gz、_1.fastq.gz 或 _1.fq.gz 结尾的文件
shopt -s nullglob  # 若通配符找不到文件，则数组为空，而不是返回原字符串
R1_FILES=("$fastq_dir"/*{.R1.fastq.gz,.R1.fq.gz,_1.fastq.gz,_1.fq.gz})

if [ ${#R1_FILES[@]} -eq 0 ]; then
  echo "在目录 $fastq_dir 下没有检测到 R1 FASTQ 文件"
  exit 1
fi

############################################
#   主循环：针对每个 R1 文件进行处理
############################################
for R1 in "${R1_FILES[@]}"; do
  # 推断 R2 文件路径
  R2=$(replace_R1_with_R2 "$R1")

  # 如果没有找到对应的R2，就给一个警告，跳过该样本
  if [ ! -f "$R2" ]; then
    echo "Warning: 对应的R2文件不存在: $R2"
    echo "跳过此样本: $R1"
    continue
  fi

  # 根据 R1 文件名提取样本名。例如:
  #   /path/xxx.R1.fastq.gz -> xxx
  #   /path/xxx_1.fq.gz     -> xxx
  # 具体实现可以灵活调整，这里简单示例：
  fname=$(basename "$R1")                # 取文件名
  sample=${fname%%.R1.fastq.gz}         # 先尝试去掉 .R1.fastq.gz
  sample=${sample%%.R1.fq.gz}           # 再尝试去掉 .R1.fq.gz
  sample=${sample%%_1.fastq.gz}         # 再尝试去掉 _1.fastq.gz
  sample=${sample%%_1.fq.gz}            # 最后尝试去掉 _1.fq.gz

  echo "============================================="
  echo "  Processing Sample: ${sample}"
  echo "  R1: $R1"
  echo "  R2: $R2"
  echo "============================================="

  ###################################################
  # Run1: 质量控制和过滤
  ###################################################
  echo "[Run1] Quality control and trimming..."
  trim_galore --paired --fastqc --cores 8 --output_dir trimmed/ "$R1" "$R2"

  # 更新FASTQ文件路径为修剪后的文件
  fastq1_trimmed="./trimmed/${sample}.R1_val_1.fq.gz"
  fastq2_trimmed="./trimmed/${sample}.R2_val_2.fq.gz"

  ###################################################
  # Run2: 比对及排序
  ###################################################
  echo "[Run2] Alignment and sorting..."
  bowtie2 \
    -p 10 \
    -X 1000 \
    -x "$genome" \
    -1 "$fastq1_trimmed" \
    -2 "$fastq2_trimmed" \
    | samtools sort -O bam -@ 5 -o "./bams/${sample}_Run2.bam"

  ###################################################
  # Run3: 去除PCR重复
  ###################################################
  echo "[Run3] Removing PCR duplicates..."
  sambamba markdup -r "./bams/${sample}_Run2.bam" "./bams/${sample}_Run3.sambamba.rmdup.bam"

  ###################################################
  # Run4: 去除线粒体、低质量序列
  ###################################################
  echo "[Run4] Removing mitochondrial reads and low-quality reads..."
  samtools view -h -f 2 -q 30 "./bams/${sample}_Run3.sambamba.rmdup.bam" \
    | grep -v chrM \
    | samtools sort -O bam -@ 5 -o "./bams/${sample}_Run4.last.bam"

  ###################################################
  # Run5: 测序深度、覆盖度、比对率、重复率统计
  ###################################################
  echo "[Run5] Calculating stats..."
  samtools index "./bams/${sample}_Run3.sambamba.rmdup.bam"
  samtools flagstat "./bams/${sample}_Run3.sambamba.rmdup.bam" \
    > "./stats/${sample}_Run5.rmdup.stat"

  ###################################################
  # Run6: BAM 转 BED
  ###################################################
  echo "[Run6] Converting BAM to BED..."
  bedtools bamtobed -i "./bams/${sample}_Run4.last.bam" \
    > "./beds/${sample}_Run6.last.bed"

  ###################################################
  # Run7: MACS2 峰调用
  ###################################################
  echo "[Run7] Calling peaks with MACS2..."
  macs2 callpeak \
    -t "./beds/${sample}_Run6.last.bed" \
    -g hs \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    -n "${sample}_Run7" \
    --outdir "./peaks/"

  ###################################################
  # Run8: 使用deeptools进行可视化
  ###################################################
  echo "[Run8] Visualizing with deeptools..."
  bamCoverage \
    --normalizeUsing CPM \
    -b "./bams/${sample}_Run4.last.bam" \
    -o "./bws/${sample}_Run8.last.bw" &

  echo "Done with sample: ${sample}"
  echo
done

###################################################
# Run9: 汇总质量控制结果 (MultiQC)
###################################################
echo "[Run9] Summarizing quality control results with MultiQC..."
multiqc ./

echo "All samples processed successfully!"
