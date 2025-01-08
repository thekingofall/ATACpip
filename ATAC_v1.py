import argparse
import glob
import os
import subprocess

def identify_samples(fastq_dir):
    """
    Identifies sample names from a directory of FASTQ files.

    Args:
        fastq_dir: The directory containing the FASTQ files.

    Returns:
        A list of sample names.
    """
    samples = set()
    for fastq_file in glob.glob(os.path.join(fastq_dir, "*.fastq.gz")) + glob.glob(os.path.join(fastq_dir, "*.fq.gz")):
        basename = os.path.basename(fastq_file)
        if ".R1." in basename:
            sample_name = basename.split(".R1.")[0]
        elif "_1." in basename:
            sample_name = basename.split("_1.")[0]
        elif ".R2." in basename:
            sample_name = basename.split(".R2.")[0]
        elif "_2." in basename:
            sample_name = basename.split("_2.")[0]
        else:
            continue
        samples.add(sample_name)
    return list(samples)

def run_atac_seq_pipeline(fastq_dir, genome, output_dir, threads):
    """
    Runs the ATAC-seq pipeline.

    Args:
        fastq_dir: The directory containing the FASTQ files.
        genome: The path to the Bowtie2 index files.
        output_dir: The directory to store the output files.
        threads: The number of threads to use.
    """
    samples = identify_samples(fastq_dir)

    # Create output directories
    os.makedirs(os.path.join(output_dir, "trimmed"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "bams"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "stats"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "beds"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "peaks"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "bws"), exist_ok=True)

    for sample in samples:
        print(f"Processing {sample}...")

        # Define FASTQ file paths
        fastq1 = ""
        fastq2 = ""

        for ext in [".R1.fastq.gz", "_1.fastq.gz", ".R1.fq.gz", "_1.fq.gz"]:
            potential_fastq1 = os.path.join(fastq_dir, sample + ext)
            if os.path.exists(potential_fastq1):
                fastq1 = potential_fastq1
                break
        for ext in [".R2.fastq.gz", "_2.fastq.gz", ".R2.fq.gz", "_2.fq.gz"]:
            potential_fastq2 = os.path.join(fastq_dir, sample + ext)
            if os.path.exists(potential_fastq2):
                fastq2 = potential_fastq2
                break

        if not fastq1 or not fastq2:
            print(f"Error: Could not find paired FASTQ files for sample {sample}. Skipping.")
            continue

        # Run1: Quality control and filtering
        print("Run1: Quality control and trimming...")
        subprocess.run([
            "trim_galore",
            "--paired",
            "--fastqc",
            "--cores", str(threads),
            "--output_dir", os.path.join(output_dir, "trimmed"),
            fastq1,
            fastq2
        ])

        # Update FASTQ file paths to trimmed files
        fastq1_trimmed = os.path.join(output_dir, "trimmed", f"{sample}_R1_val_1.fq.gz")
        fastq2_trimmed = os.path.join(output_dir, "trimmed", f"{sample}_R2_val_2.fq.gz")

        # Run2: Alignment and sorting
        print("Run2: Alignment and sorting...")
        with open(os.path.join(output_dir, "bams", f"{sample}_Run2.bam"), "w") as outfile:
            subprocess.run([
                "bowtie2",
                "-p", str(threads),
                "-X", "1000",
                "-x", genome,
                "-1", fastq1_trimmed,
                "-2", fastq2_trimmed
            ], stdout=subprocess.PIPE, check=True)
            subprocess.run([
                "samtools",
                "sort",
                "-O", "bam",
                "-@", str(threads//2 if threads//2 > 0 else 1)
            ], stdin=subprocess.PIPE, stdout=outfile, check=True)

        # Run3: Remove PCR duplicates
        print("Run3: Removing PCR duplicates...")
        subprocess.run([
            "sambamba",
            "markdup",
            "-r",
            os.path.join(output_dir, "bams", f"{sample}_Run2.bam"),
            os.path.join(output_dir, "bams", f"{sample}_Run3.sambamba.rmdup.bam")
        ])

        # Run4: Remove mitochondrial genes and low-quality sequences
        print("Run4: Removing mitochondrial genes and low-quality sequences...")
        with open(os.path.join(output_dir, "bams", f"{sample}_Run4.last.bam"), "w") as outfile:
            p1 = subprocess.Popen([
                "samtools",
                "view",
                "-h",
                "-f", "2",
                "-q", "30",
                os.path.join(output_dir, "bams", f"{sample}_Run3.sambamba.rmdup.bam")
            ], stdout=subprocess.PIPE)
            p2 = subprocess.Popen([
                "grep",
                "-v",
                "chrM"
            ], stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()
            subprocess.run([
                "samtools",
                "sort",
                "-O", "bam",
                "-@", str(threads//2 if threads//2 > 0 else 1)
            ], stdin=p2.stdout, stdout=outfile)
            p2.stdout.close()

        # Run5: Calculate sequencing depth, coverage, alignment rate, and duplication rate
        print("Run5: Calculating sequencing depth, coverage, alignment rate, and duplication rate...")
        subprocess.run([
            "samtools",
            "index",
            os.path.join(output_dir, "bams", f"{sample}_Run3.sambamba.rmdup.bam")
        ])
        with open(os.path.join(output_dir, "stats", f"{sample}_Run5.rmdup.stat"), "w") as outfile:
            subprocess.run([
                "samtools",
                "flagstat",
                os.path.join(output_dir, "bams", f"{sample}_Run3.sambamba.rmdup.bam")
            ], stdout=outfile)

        # Run6: Convert BAM to BED
        print("Run6: Converting BAM to BED...")
        subprocess.run([
            "bedtools",
            "bamtobed",
            "-i", os.path.join(output_dir, "bams", f"{sample}_Run4.last.bam"),
            ">", os.path.join(output_dir, "beds", f"{sample}_Run6.last.bed")
        ], shell=True)

        # Run7: Call peaks with MACS2
        print("Run7: Calling peaks with MACS2...")
        subprocess.run([
            "macs2",
            "callpeak",
            "-t", os.path.join(output_dir, "beds", f"{sample}_Run6.last.bed"),
            "-g", "hs",
            "--nomodel",
            "--shift", "-100",
            "--extsize", "200",
            "-n", f"{sample}_Run7",
            "--outdir", os.path.join(output_dir, "peaks")
        ])

        # Run8: Visualize with deeptools
        print("Run8: Visualizing with deeptools...")
        subprocess.Popen([
            "bamCoverage",
            "--normalizeUsing", "CPM",
            "-b", os.path.join(output_dir, "bams", f"{sample}_Run4.last.bam"),
            "-o", os.path.join(output_dir, "bws", f"{sample}_Run8.last.bw")
        ])

    # Run9: Summarize quality control results with MultiQC
    print("Run9: Summarizing quality control results with MultiQC...")
    subprocess.run(["multiqc", output_dir])

    print("All samples processed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ATAC-seq data processing pipeline.")
    parser.add_argument("fastq_dir", help="Directory containing the FASTQ files.")
    parser.add_argument("genome", help="Path to the Bowtie2 index files.")
    parser.add_argument("output_dir", help="Directory to store the output files.")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads to use (default: 8).")

    args = parser.parse_args()

    run_atac_seq_pipeline(args.fastq_dir, args.genome, args.output_dir, args.threads)
