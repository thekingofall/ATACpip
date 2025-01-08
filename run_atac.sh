#!/bin/bash

# ---------------------------------------------------------
# Script to process ATAC-seq data
# ---------------------------------------------------------

# Usage: ./process_atac.sh -i <fastq_directory> -g <genome_bowtie2_index>

# Required parameters
fastq_dir=""
genome=""

# Output directories
trimmed_dir="trimmed"
bams_dir="bams"
stats_dir="stats"
beds_dir="beds"
peaks_dir="peaks"
bws_dir="bws"

# ---------------------------------------------------------
# Helper function to print usage instructions
# ---------------------------------------------------------
print_usage() {
  echo "Usage: ./process_atac.sh -i <fastq_directory> -g <genome_bowtie2_index>"
  echo ""
  echo "Required parameters:"
  echo "  -i <fastq_directory>      Path to the directory containing the FASTQ files."
  echo "  -g <genome_bowtie2_index> Path to the Bowtie2 index of the reference genome (without the .bt2 extension)."
  echo ""
  exit 1
}

# ---------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------
while getopts "i:g:" opt; do
  case "$opt" in
    i) fastq_dir="$OPTARG" ;;
    g) genome="$OPTARG" ;;
    \?) print_usage ;;
  esac
done

# Check if required parameters are provided
if [ -z "$fastq_dir" ] || [ -z "$genome" ]; then
  echo "Error: Missing required parameters."
  print_usage
fi

# Check if fastq directory exists
if [ ! -d "$fastq_dir" ]; then
  echo "Error: FASTQ directory not found: $fastq_dir"
  exit 1
fi

# Check if genome index exists (basic check, assuming .1.bt2 exists)
if [ ! -f "$genome.1.bt2" ]; then
  echo "Error: Genome Bowtie2 index not found: $genome"
  echo "Make sure you provide the path without the .bt2 extension."
  exit 1
fi

# ---------------------------------------------------------
# Create output directories
# ---------------------------------------------------------
mkdir -p "$trimmed_dir" "$bams_dir" "$stats_dir" "$beds_dir" "$peaks_dir" "$bws_dir"

# ---------------------------------------------------------
# Function to extract sample name from FASTQ file name
# ---------------------------------------------------------
extract_sample_name() {
  local filename="\$1"
  local base=$(basename "$filename")
  # Try different patterns to extract the sample name
  if [[ "$base" =~ (.*)\.R[12]\.fastq\.gz ]]; then
    echo "${BASH_REMATCH[1]}"
  elif [[ "$base" =~ (.*)\.R[12]\.fq\.gz ]]; then
    echo "${BASH_REMATCH[1]}"
  elif [[ "$base" =~ (.*)_[12]\.fastq\.gz ]]; then
    echo "${BASH_REMATCH[1]}"
  elif [[ "$base" =~ (.*)_[12]\.fq\.gz ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "Unknown"
  fi
}

# ---------------------------------------------------------
# Find all potential FASTQ files
# ---------------------------------------------------------
find "$fastq_dir" -maxdepth 1 \( -name "*.R1.fastq.gz" -o -name "*.R1.fq.gz" -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) -print0 | while IFS= read -r -d $'\0' f1; do
  # Extract sample name from R1 file
  sample=$(extract_sample_name "$f1")

  # Skip if sample name couldn't be extracted
  if [ "$sample" == "Unknown" ]; then
    echo "Warning: Could not extract sample name from $f1. Skipping."
    continue
  fi

  # Construct potential R2 file names based on different patterns
  r2_pattern1="${fastq_dir}/${sample}.R2.fastq.gz"
  r2_pattern2="${fastq_dir}/${sample}.R2.fq.gz"
  r2_pattern3="${fastq_dir}/${sample}_2.fastq.gz"
  r2_pattern4="${fastq_dir}/${sample}_2.fq.gz"

  # Find the corresponding R2 file
  f2=""
  if [ -f "$r2_pattern1" ]; then
    f2="$r2_pattern1"
  elif [ -f "$r2_pattern2" ]; then
    f2="$r2_pattern2"
  elif [ -f "$r2_pattern3" ]; then
    f2="$r2_pattern3"
  elif [ -f "$r2_pattern4" ]; then
    f2="$r2_pattern4"
  fi

  # If corresponding R2 file is found, process the sample
  if [ -n "$f2" ]; then
    echo "Processing sample: ${sample}"

    # Define FASTQ file paths
    fastq1="$f1"
    fastq2="$f2"

    # Run1: 质量控制和过滤
    echo "Run1: Quality control and trimming..."
    trim_galore --paired --fastqc --cores 8 --output_dir "$trimmed_dir" "$fastq1" "$fastq2"

    # Update FASTQ file paths for trimmed files
    fastq1_trimmed="$trimmed_dir/${sample}.R1_val_1.fq.gz"
    fastq2_trimmed="$trimmed_dir/${sample}.R2_val_2.fq.gz"

    # Run2: 比对及排序
    echo "Run2: Alignment and sorting..."
    bowtie2 -p 10 -X 1000 -x "$genome" -1 "$fastq1_trimmed" -2 "$fastq2_trimmed" | samtools sort -O bam -@ 5 -o "$bams_dir/${sample}_Run2.bam"

    # Run3: 去除PCR重复
    echo "Run3: Removing PCR duplicates..."
    sambamba markdup -r "$bams_dir/${sample}_Run2.bam" "$bams_dir/${sample}_Run3.sambamba.rmdup.bam"

    # Run4: 去除线粒体基因，去除低质量序列
    echo "Run4: Removing mitochondrial genes and low-quality sequences..."
    samtools view -h -f 2 -q 30 "$bams_dir/${sample}_Run3.sambamba.rmdup.bam" | grep -v chrM | samtools sort -O bam -@ 5 -o "$bams_dir/${sample}_Run4.last.bam"

    # Run5: 对测序深度、覆盖度、比对率、重复率进行统计
    echo "Run5: Calculating sequencing depth, coverage, alignment rate, and duplication rate..."
    samtools index "$bams_dir/${sample}_Run3.sambamba.rmdup.bam"
    samtools flagstat "$bams_dir/${sample}_Run3.sambamba.rmdup.bam" > "$stats_dir/${sample}_Run5.rmdup.stat"

    # Run6: 将bam文件转换为bed文件
    echo "Run6: Converting BAM to BED..."
    bedtools bamtobed -i "$bams_dir/${sample}_Run4.last.bam" > "$beds_dir/${sample}_Run6.last.bed"

    # Run7: 使用MACS2进行peak calling
    echo "Run7: Calling peaks with MACS2..."
    macs2 callpeak -t "$beds_dir/${sample}_Run6.last.bed" -g hs --nomodel --shift -100 --extsize 200 -n "${sample}_Run7" --outdir "$peaks_dir"

    # Run8: 使用deeptools进行可视化
    echo "Run8: Visualizing with deeptools..."
    bamCoverage --normalizeUsing CPM -b "$bams_dir/${sample}_Run4.last.bam" -o "$bws_dir/${sample}_Run8.last.bw" &
  else
    echo "Warning: Could not find matching R2 file for sample ${sample} (R1 file: $f1). Skipping."
  fi
done

# Wait for all background bamCoverage processes to finish
wait

# Run9: 使用multiqc进行质量控制结果的汇总
echo "Run9: Summarizing quality control results with MultiQC..."
multiqc ./

echo "All samples processed."
