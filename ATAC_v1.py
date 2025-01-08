#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import subprocess
from pathlib import Path

def parse_args():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(
        description="ATAC-seq分析流程：自动识别输入文件夹中配对的fastq文件，并执行一系列分析步骤。"
    )
    parser.add_argument(
        "-i", "--input_dir", 
        required=True, 
        help="输入fastq文件所在文件夹路径。"
    )
    parser.add_argument(
        "-g", "--genome_index",
        default="/home/maolp/mao/Ref/AllnewstarRef/Homo/HG19/HG19BT/hg19",
        help="Bowtie2的参考基因组索引前缀(默认: /home/maolp/mao/Ref/AllnewstarRef/Homo/HG19/HG19BT/hg19)"
    )
    parser.add_argument(
        "-f", "--fastqc_cores",
        type=int,
        default=8,
        help="执行trim_galore时FastQC使用的线程数(默认: 8)"
    )
    parser.add_argument(
        "-b", "--bowtie2_cores",
        type=int,
        default=10,
        help="执行bowtie2时使用的线程数(默认: 10)"
    )
    parser.add_argument(
        "-s", "--samtools_cores",
        type=int,
        default=5,
        help="执行samtools sort时使用的线程数(默认: 5)"
    )
    parser.add_argument(
        "-p", "--peak_genome_size",
        default="hs",
        help="MACS2物种基因组大小参数(默认: hs)"
    )
    parser.add_argument(
        "--no_run",
        action="store_true",
        help="如果指定，则只打印命令而不实际执行，用于调试。"
    )
    return parser.parse_args()

def find_paired_fastq_files(input_dir):
    """
    在 input_dir 中查找所有支持的fastq文件，自动配对R1和R2。
    支持的命名后缀：
      - .R1.fastq.gz / .R2.fastq.gz
      - .R1.fq.gz   / .R2.fq.gz
      - _1.fastq.gz / _2.fastq.gz
      - _1.fq.gz    / _2.fq.gz
    
    返回：{
        sample_name1: (r1_file_path, r2_file_path),
        sample_name2: (r1_file_path, r2_file_path),
        ...
    }
    """
    # 匹配模式，用来捕获 (prefix) + (R1/R2 或 1/2) + (fastq/fq) + (gz)
    pattern = re.compile(
        r"^(?P<sample>.+?)"                 # 样本前缀，非贪婪匹配
        r"(?:[._](?:R|_)?(?P<read>[12]))"   # .R1  / _1 / .R2 / _2等
        r"\.(?:fastq|fq)\.gz$"             # .fastq.gz 或 .fq.gz
    )

    paired_files = {}

    input_dir = Path(input_dir)
    for file_path in input_dir.iterdir():
        if not file_path.is_file():
            continue
        match = pattern.match(file_path.name)
        if match:
            sample = match.group("sample")
            read = match.group("read")  # '1' or '2'
            if sample not in paired_files:
                paired_files[sample] = ["", ""]  # 占位

            # 根据 read 来放置
            if read == "1":
                paired_files[sample][0] = str(file_path.resolve())
            elif read == "2":
                paired_files[sample][1] = str(file_path.resolve())

    # 过滤掉没有配对成功的样本
    paired_files = {
        k: tuple(v)
        for k, v in paired_files.items()
        if v[0] and v[1]
    }

    return paired_files

def run_command(cmd, no_run=False):
    """
    执行给定的命令行。如果 no_run=True，则仅打印而不执行
    """
    print(f"[CMD] {cmd}")
    if not no_run:
        subprocess.run(cmd, shell=True, check=True)

def main():
    args = parse_args()

    # 1. 找到所有的R1-R2成对文件
    sample_dict = find_paired_fastq_files(args.input_dir)
    if not sample_dict:
        print("未找到任何可配对的FASTQ文件，请检查输入文件夹。")
        return

    # 创建输出文件夹
    for outdir in ["trimmed", "bams", "stats", "beds", "peaks", "bws"]:
        os.makedirs(outdir, exist_ok=True)

    # 2. 循环处理每个样本
    for sample, (fastq1, fastq2) in sample_dict.items():
        print(f"\n>>> Processing sample: {sample}")
        print(f"    R1: {fastq1}")
        print(f"    R2: {fastq2}")

        # ================ Run1: 质量控制和过滤 ================
        print("Run1: Quality control and trimming...")
        # trim_galore --paired --fastqc --cores X --output_dir trimmed/ R1 R2
        cmd_trim = (
            f"trim_galore --paired --fastqc --cores {args.fastqc_cores} "
            f"--output_dir trimmed/ {fastq1} {fastq2}"
        )
        run_command(cmd_trim, no_run=args.no_run)

        # trimmed文件名：默认Trim Galore会输出 {basename}_val_1.fq.gz / val_2.fq.gz
        # 例如：sample.R1_val_1.fq.gz。我们只需要把原始前缀替换成trimmed目录下。
        base1 = os.path.basename(fastq1)
        base2 = os.path.basename(fastq2)
        # Trim Galore 通常的输出：base1 -> base1的前缀+"_val_1.fq.gz"
        # 举例： O1-1.R1.fastq.gz -> O1-1.R1_val_1.fq.gz
        # 这里直接拼出新的文件名（也可以用更严格的方法检测）
        fastq1_trimmed = os.path.join("trimmed", base1.replace(".R1.", ".R1_val_1.").replace("_1.", "_1_val_1."))
        fastq2_trimmed = os.path.join("trimmed", base2.replace(".R2.", ".R2_val_2.").replace("_2.", "_2_val_2."))

        # ================ Run2: 比对及排序 ================
        print("Run2: Alignment and sorting...")
        # bowtie2 -p 10 -X 1000 -x {genome} -1 R1 -2 R2 | samtools sort -O bam -@ 5 -o bams/sample_Run2.bam
        out_bam_run2 = f"./bams/{sample}_Run2.bam"
        cmd_align_sort = (
            f"bowtie2 -p {args.bowtie2_cores} -X 1000 -x {args.genome_index} "
            f"-1 {fastq1_trimmed} -2 {fastq2_trimmed} | "
            f"samtools sort -O bam -@ {args.samtools_cores} -o {out_bam_run2}"
        )
        run_command(cmd_align_sort, no_run=args.no_run)

        # ================ Run3: 去除PCR重复 ================
        print("Run3: Removing PCR duplicates...")
        out_bam_run3 = f"./bams/{sample}_Run3.sambamba.rmdup.bam"
        cmd_rmdup = (
            f"sambamba markdup -r {out_bam_run2} {out_bam_run3}"
        )
        run_command(cmd_rmdup, no_run=args.no_run)

        # ================ Run4: 去除线粒体基因 & 低质量序列 ================
        print("Run4: Removing mitochondrial genes and low-quality sequences...")
        out_bam_run4 = f"./bams/{sample}_Run4.last.bam"
        cmd_filter_mt = (
            f"samtools view -h -f 2 -q 30 {out_bam_run3} | "
            f"grep -v chrM | "
            f"samtools sort -O bam -@ {args.samtools_cores} -o {out_bam_run4}"
        )
        run_command(cmd_filter_mt, no_run=args.no_run)

        # ================ Run5: 统计测序深度、覆盖度、比对率、重复率等 ================
        print("Run5: Calculating sequencing metrics (depth, coverage, alignment rate, duplication rate)...")
        # samtools index & samtools flagstat
        cmd_index = f"samtools index {out_bam_run3}"
        run_command(cmd_index, no_run=args.no_run)

        stat_file = f"./stats/{sample}_Run5.rmdup.stat"
        cmd_flagstat = f"samtools flagstat {out_bam_run3} > {stat_file}"
        run_command(cmd_flagstat, no_run=args.no_run)

        # ================ Run6: 将bam文件转换为bed文件 ================
        print("Run6: Converting BAM to BED...")
        out_bed_run6 = f"./beds/{sample}_Run6.last.bed"
        cmd_bamtobed = f"bedtools bamtobed -i {out_bam_run4} > {out_bed_run6}"
        run_command(cmd_bamtobed, no_run=args.no_run)

        # ================ Run7: 使用MACS2进行peak calling ================
        print("Run7: Calling peaks with MACS2...")
        # macs2 callpeak -t bedfile -g hs --nomodel --shift -100 --extsize 200 -n sample_Run7 --outdir peaks/
        cmd_macs2 = (
            f"macs2 callpeak "
            f"-t {out_bed_run6} "
            f"-g {args.peak_genome_size} "
            f"--nomodel --shift -100 --extsize 200 "
            f"-n {sample}_Run7 --outdir ./peaks/"
        )
        run_command(cmd_macs2, no_run=args.no_run)

        # ================ Run8: 使用deeptools进行可视化 ================
        print("Run8: Visualizing with deeptools (bamCoverage)...")
        out_bw = f"./bws/{sample}_Run8.last.bw"
        cmd_bamcoverage = (
            f"bamCoverage --normalizeUsing CPM "
            f"-b {out_bam_run4} "
            f"-o {out_bw} &"
        )
        # 注意：后面加 & 代表在后台执行，如果希望严格顺序执行，可去掉 &。
        run_command(cmd_bamcoverage, no_run=args.no_run)

    # ================ Run9: 使用multiqc进行质量控制结果的汇总 ================
    print("\nRun9: Summarizing quality control results with MultiQC...")
    cmd_multiqc = "multiqc ./"
    run_command(cmd_multiqc, no_run=args.no_run)

    print("\nAll samples processed.")

if __name__ == "__main__":
    main()
