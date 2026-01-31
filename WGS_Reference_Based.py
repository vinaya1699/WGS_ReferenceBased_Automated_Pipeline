#!/usr/bin/env python3

# Example Usage :
# python Genotype_By_Sequencing.py -d /path/to/workdir -org Phaseolus_vulgaris --qual 40 --min_dp 10 --f_missing 0.7 --min_af 0.1

# Please make sure to have all the required tools installed and accessible in your PATH.
# All the snp and indel filtration parameters are as per GATK pipeline. You can make changes if required.


import pandas as pd
import numpy as np
import os
import sys
import glob
import subprocess
import argparse
import datetime
import time
import shutil

print(r"""
__        __   ____    _____
\ \      / /  / ___|  / ____|
 \ \ /\ / /  | |  _  | (___
  \ V  V /   | |_| |  \___ \
   \_/\_/     \____|  |____/

        Reference Based Whole Genome Sequencing (WGS) Automation Pipeline
        : Vinaya Kadam
""")

script_start = time.time()
print(f"Script started at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# ---------------- ARGUMENTS ---------------- #
parser = argparse.ArgumentParser(description="GBS Automation Script")
parser.add_argument('-d', '--Working_Directory', type=str, required=True, help='Input working directory containing raw FASTQ files')
parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads')
parser.add_argument('-org', '--organism', type=str, required=True, help='Organism name (reference fasta prefix)')

args = parser.parse_args()
threads = args.threads
organism = args.organism

# ---------------- SETUP ---------------- #
os.chdir(args.Working_Directory)
print(f"Changed working directory to: {args.Working_Directory}")

ref_dir = "0_Reference_Genome"
ref_fasta = f"{ref_dir}/{organism}.fasta"
picard = "/Analysis3/Vinaya/picard.jar"
gatk = "/apps/gatk-4.2.6.1/gatk"
vcf_dir = "4_Variant_Calling"
snpEff = "/apps/snpEff5.0/snpEff.jar"
snpEff_data = "/apps/snpEff5.0/data"

Input_Raw_files = glob.glob("1_RawData/*_R1.fastq.gz")
files_series = pd.Series(Input_Raw_files)
samples = files_series.str.replace(r'.*/|_R1\.fastq\.gz$', '', regex=True).unique()
print("Sample names extracted:", samples)

folders = [
    "2_Clean_data",
    "3_Alignment",
    "1_RawData/Fastqc_Output",
    "2_Clean_data/Fastqc_Output",
    "4_Variant_Calling",
    "5_Variant_Filtration",
    "6_Variant_Annotation",
    "logs"
]
for folder in folders:
    os.makedirs(folder, exist_ok=True)
    os.chmod(folder, 0o777)

log_dir = "logs"

# ---------------- HELPER FUNCTION ---------------- #
def run_command(cmd, log_file):
    """Run a command and log stdout+stderr to file and console."""
    with open(log_file, "w") as logfile:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        for line in process.stdout:
            print(line, end="")       # Console
            logfile.write(line)       # Log file
        process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Command failed! Check log: {log_file}")

# ---------------- FASTQC ---------------- #
print("############################### fastqc ##################################")
for sample in samples:
    r1 = f"1_RawData/{sample}_R1.fastq.gz"
    r2 = f"1_RawData/{sample}_R2.fastq.gz"
    log = f"{log_dir}/{sample}_fastqc_raw.log"
    cmd = ["fastqc", "-t", str(threads), r1, r2, "-o", "1_RawData/Fastqc_Output/"]
    run_command(cmd, log)

# ---------------- FASTP ---------------- #
print("############################### fastp ##################################")
for sample in samples:
    r1 = f"1_RawData/{sample}_R1.fastq.gz"
    r2 = f"1_RawData/{sample}_R2.fastq.gz"
    r1_clean = f"2_Clean_data/{sample}_R1.fq.gz"
    r2_clean = f"2_Clean_data/{sample}_R2.fq.gz"
    html = f"2_Clean_data/{sample}_fastp.html"
    json = f"2_Clean_data/{sample}_fastp.json"
    log = f"{log_dir}/{sample}_fastp.log"
    cmd = [
        "fastp", "-i", r1, "-I", r2, "-o", r1_clean, "-O", r2_clean,
        "-h", html, "-j", json, "-w", str(threads), "-c", "--detect_adapter_for_pe"
    ]
    run_command(cmd, log)

# ---------------- CLEAN DATA QC ---------------- #
print("############################### Clean Data QC ##################################")
for sample in samples:
    r1 = f"2_Clean_data/{sample}_R1.fq.gz"
    r2 = f"2_Clean_data/{sample}_R2.fq.gz"
    log = f"{log_dir}/{sample}_fastqc_clean.log"
    cmd = ["fastqc", "-t", str(threads), r1, r2, "-o", "2_Clean_data/Fastqc_Output/"]
    run_command(cmd, log)

# ---------------- INDEX REFERENCE ---------------- #
print("############################### Indexing Ref. Genome ##################################")
index_extensions = [".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"]
missing_files = []
for ext in index_extensions:
    file_path = os.path.join(ref_dir, os.path.splitext(os.path.basename(ref_fasta))[0] + ext)
    if not os.path.exists(file_path):
        missing_files.append(file_path)

if missing_files:
    print(f"Missing index files: {missing_files}")
    log = f"{log_dir}/bwa_index.log"
    run_command(["bwa", "index", ref_fasta], log)
else:
    print("All index files found. Reference genome already indexed.")

# ---------------- ALIGNMENT ---------------- #
print("############################### Alignment & MarkDuplicates ##################################")
alignment_dir = "3_Alignment"
os.makedirs(alignment_dir, exist_ok=True)

for sample in samples:
    r1 = f"2_Clean_data/{sample}_R1.fq.gz"
    r2 = f"2_Clean_data/{sample}_R2.fq.gz"
    bam = f"{alignment_dir}/{sample}.bam"
    marked_bam = f"{alignment_dir}/{sample}_sorted_rdgrp_dup_marked.bam"
    metrics = f"{alignment_dir}/{sample}_dup-metrics.txt"
    read_group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA"

    # ---------------- Mapping + Sorting ---------------- #
    log_align = f"{log_dir}/{sample}_alignment.log"
    print(f"Running BWA MEM + samtools sort for sample: {sample}")
    with open(log_align, "w") as logfile:
        bwa_cmd = ["bwa", "mem", "-R", read_group, "-M", "-t", str(threads), ref_fasta, r1, r2]
        samtools_cmd = ["samtools", "sort", "-@", str(threads), "-o", bam, "-"]

        # Run BWA
        bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = bwa_proc.communicate()
        logfile.write("### BWA stdout ###\n")
        logfile.write(stdout or "")
        logfile.write("\n### BWA stderr ###\n")
        logfile.write(stderr or "")
        if bwa_proc.returncode != 0 or not stdout:
            raise RuntimeError(f"BWA MEM failed for {sample}. Check {log_align}")

        # Run samtools sort
        samtools_proc = subprocess.Popen(samtools_cmd, stdin=subprocess.PIPE, stdout=logfile, stderr=logfile, text=True)
        samtools_stdout, samtools_stderr = samtools_proc.communicate(input=stdout)
        if samtools_proc.returncode != 0:
            raise RuntimeError(f"samtools sort failed for {sample}. Check {log_align}")

    # ---------------- MarkDuplicates ---------------- #
    log_markdup = f"{log_dir}/{sample}_markduplicates.log"
    markdup_cmd = [
        "java", "-Xmx16g", "-jar", picard, "MarkDuplicates",
        f"I={bam}", f"O={marked_bam}", f"M={metrics}", "QUIET=TRUE"
    ]
    run_command(markdup_cmd, log_markdup)
    print(f"Final BAM with duplicates marked: {marked_bam}")

# ---------------- GATK HaplotypeCaller ---------------- #
print("######################## GATK HaplotypeCaller : Variant Calling########################")
os.makedirs(vcf_dir, exist_ok=True)

# Index fasta if missing
fai_file = f"{ref_fasta}.fai"
if not os.path.exists(fai_file):
    run_command(["samtools", "faidx", ref_fasta], f"{log_dir}/samtools_faidx.log")

dict_file = os.path.splitext(ref_fasta)[0] + ".dict"
if not os.path.exists(dict_file):
    run_command(["java", "-jar", picard, "CreateSequenceDictionary", f"R={ref_fasta}", f"O={dict_file}"], f"{log_dir}/picard_dict.log")

for sample in samples:
    marked_bam = f"3_Alignment/{sample}_sorted_rdgrp_dup_marked.bam"
    run_command(["samtools", "index", marked_bam], f"{log_dir}/{sample}_bam_index.log")
    vcf = f"{vcf_dir}/{sample}.vcf"
    log = f"{log_dir}/{sample}_haplotypecaller.log"
    cmd = [gatk, "HaplotypeCaller", "-I", marked_bam, "-R", ref_fasta, "-O", vcf]
    run_command(cmd, log)

print("############################### Variant Filtration ##################################")

# ...after variant calling...
filtration_dir = "5_Variant_Filtration"
os.makedirs(filtration_dir, exist_ok=True)
for sample in samples:
    vcf = f"{vcf_dir}/{sample}.vcf"
    # Separate SNPs
    run_command([gatk, "SelectVariants", "-V", vcf, "-select-type", "SNP", "-O", f"{filtration_dir}/{sample}_snps.vcf.gz"], f"{log_dir}/{sample}_select_snps.log")
    # Separate INDELs
    run_command([gatk, "SelectVariants", "-V", vcf, "-select-type", "INDEL", "-O", f"{filtration_dir}/{sample}_indels.vcf.gz"], f"{log_dir}/{sample}_select_indels.log")
    # Filter SNPs
    run_command([
        gatk, "VariantFiltration", "-V", f"{filtration_dir}/{sample}_snps.vcf.gz",
        "-filter", "QD < 2.0", "--filter-name", "QD2",
        "-filter", "QUAL < 30.0", "--filter-name", "QUAL30",
        "-filter", "SOR > 3.0", "--filter-name", "SOR3",
        "-filter", "FS > 60.0", "--filter-name", "FS60",
        "-filter", "MQ < 40.0", "--filter-name", "MQ40",
        "-filter", "MQRankSum < -12.5", "--filter-name", "MQRankSum-12.5",
        "-filter", "ReadPosRankSum < -8.0", "--filter-name", "ReadPosRankSum-8",
        "-O", f"{filtration_dir}/{sample}_snps_filtered.vcf.gz"
    ], f"{log_dir}/{sample}_filter_snps.log")
    # Filter INDELs
    run_command([
        gatk, "VariantFiltration", "-V", f"{filtration_dir}/{sample}_indels.vcf.gz",
        "-filter", "QD < 2.0", "--filter-name", "QD2",
        "-filter", "QUAL < 30.0", "--filter-name", "QUAL30",
        "-filter", "FS > 200.0", "--filter-name", "FS200",
        "-filter", "ReadPosRankSum < -20.0", "--filter-name", "ReadPosRankSum-20",
        "-O", f"{filtration_dir}/{sample}_indels_filtered.vcf.gz"
    ], f"{log_dir}/{sample}_filter_indels.log")
    # Select PASS SNPs
    run_command([
        gatk, "SelectVariants", "-R", ref_fasta, "-V", f"{filtration_dir}/{sample}_snps_filtered.vcf.gz",
        "--exclude-filtered", "-O", f"{filtration_dir}/{sample}_snps_PASS.vcf.gz"
    ], f"{log_dir}/{sample}_pass_snps.log")
    # Select PASS INDELs
    run_command([
        gatk, "SelectVariants", "-R", ref_fasta, "-V", f"{filtration_dir}/{sample}_indels_filtered.vcf.gz",
        "--exclude-filtered", "-O", f"{filtration_dir}/{sample}_indels_PASS.vcf.gz"
    ], f"{log_dir}/{sample}_pass_indels.log")
# ...existing code...

print("########################## Merging of Variant Calls #############################")

# ...after filtration...
annotation_dir = "6_Variant_Annotation"
os.makedirs(annotation_dir, exist_ok=True)
for sample in samples:
    run_command([
        "java", "-jar", picard, "MergeVcfs",
        f"I={filtration_dir}/{sample}_snps_PASS.vcf.gz",
        f"I={filtration_dir}/{sample}_indels_PASS.vcf.gz",
        f"O={annotation_dir}/{sample}_Passed_Variants.vcf"
    ], f"{log_dir}/{sample}_mergevcf.log")

print("############################### SNP Annotation ##################################")

def ensure_snpeff_config_entry(organism, snpEff_data, snpEff_path):
    config_path = os.path.join(os.path.dirname(snpEff_path), "snpEff.config")
    entry = f"{organism}.genome : {organism}\n"
    # Check if entry exists
    with open(config_path, "r") as f:
        config = f.read()
    if f"{organism}.genome" not in config:
        print(f"Adding {organism} entry to snpEff.config")
        with open(config_path, "a") as f:
            f.write(f"\n# Added by automation\n{organism}.genome : {organism}\n")
    else:
        print(f"snpEff.config already contains entry for {organism}")

# Ensure snpEff config is correct before building DB
ensure_snpeff_config_entry(organism, snpEff_data, snpEff)

organism_snpEff_dir = os.path.join(snpEff_data, organism)

# Create organism directory if it does not exist
if not os.path.exists(organism_snpEff_dir):
    print(f"Creating directory for snpEff organism data at {organism_snpEff_dir}")
    os.makedirs(organism_snpEff_dir, exist_ok=True)

# Copy fasta to snpEff data directory, rename to sequences.fa
ref_fasta_dest = os.path.join(organism_snpEff_dir, "sequences.fa")
if not os.path.exists(ref_fasta_dest):
    print(f"Copying reference fasta {ref_fasta} to {ref_fasta_dest}")
    shutil.copy2(ref_fasta, ref_fasta_dest)
else:
    print(f"Reference fasta already exists at {ref_fasta_dest}")

# Copy GTF annotation, assuming it is alongside fasta with .gtf extension
ref_gtf = os.path.splitext(ref_fasta)[0] + ".gtf"
ref_gtf_dest = os.path.join(organism_snpEff_dir, "genes.gtf")
if not os.path.exists(ref_gtf_dest):
    print(f"Copying annotation GTF {ref_gtf} to {ref_gtf_dest}")
    shutil.copy2(ref_gtf, ref_gtf_dest)
else:
    print(f"GTF annotation already exists at {ref_gtf_dest}")

print("Reference genome and annotation copied and renamed for snpEff.")

# Build snpEff database if not already built
snpEff_db_file = os.path.join(organism_snpEff_dir, "snpEffectPredictor.bin")
if not os.path.exists(snpEff_db_file):
    print(f"Building snpEff database for organism: {organism}")
    build_cmd = [
        "java", "-jar", snpEff,
        "build", "-gtf22", "-v", organism
    ]
    process = subprocess.Popen(build_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    for line in process.stdout:
        print(line, end="")
    process.wait()
    if process.returncode != 0:
        raise RuntimeError("snpEff database build failed!")
    print("snpEff database build completed.")
else:
    print(f"snpEff database already exists for organism: {organism}")

# Annotate merged variants with stats
for sample in samples:
    vcf_input = f"{annotation_dir}/{sample}_Passed_Variants.vcf"
    vcf_output = f"{annotation_dir}/{sample}_Annotated_Variants.vcf"
    csv_stats = f"{annotation_dir}/{sample}_variants_ann.csv.stats"
    html_stats = f"{annotation_dir}/{sample}_variants_summary.html"
    log = f"{log_dir}/{sample}_Variant_annotation.log"
    cmd = [
        "java", "-jar", snpEff,
        "-v", organism,
        "-noLog", "-no-downstream", "-no-upstream", "-no-utr",
        "-o", "vcf",
        vcf_input,
        "-csvStats", csv_stats,
        "-htmlStats", html_stats
    ]
    with open(log, "w") as logfile, open(vcf_output, "w") as vcf_out:
        process = subprocess.Popen(cmd, stdout=vcf_out, stderr=subprocess.PIPE, text=True)
        for line in process.stderr:
            print(line, end="")
            logfile.write(line)
        process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"snpEff annotation failed for {sample}. Check {log}")
    print(f"Annotated VCF for sample {sample} written to {vcf_output}")
    print(f"Stats: {csv_stats}, HTML summary: {html_stats}")

# Annotate SNPs with stats
for sample in samples:
    Snp_vcf_input = f"{filtration_dir}/{sample}_snps_PASS.vcf.gz"
    Snp_vcf_output = f"{annotation_dir}/{sample}_Annotated_SNPs.vcf"
    Snp_csv_stats = f"{annotation_dir}/{sample}_snps_ann.csv.stats"
    Snp_html_stats = f"{annotation_dir}/{sample}_snps_summary.html"
    log = f"{log_dir}/{sample}_snp_annotation.log"
    cmd = [
        "java", "-jar", snpEff,
        "-v", organism,
        "-noLog", "-no-downstream", "-no-upstream", "-no-utr",
        "-o", "vcf",
        Snp_vcf_input,
        "-csvStats", Snp_csv_stats,
        "-htmlStats", Snp_html_stats
    ]
    with open(log, "w") as logfile, open(Snp_vcf_output, "w") as SNP_vcf_out:
        process = subprocess.Popen(cmd, stdout=SNP_vcf_out, stderr=subprocess.PIPE, text=True)
        for line in process.stderr:
            print(line, end="")
            logfile.write(line)
        process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Snp annotation failed for {sample}. Check {log}")
    print(f"Annotated SNPs for sample {sample} written to {Snp_vcf_output}")
    print(f"SNP stats: {Snp_csv_stats}, HTML summary: {Snp_html_stats}")

# Annotate INDELs with stats
for sample in samples:
    Indels_vcf_input = f"{filtration_dir}/{sample}_indels_PASS.vcf.gz"
    Indels_vcf_output = f"{annotation_dir}/{sample}_Annotated_Indels.vcf"
    Indels_csv_stats = f"{annotation_dir}/{sample}_indels_ann.csv.stats"
    Indels_html_stats = f"{annotation_dir}/{sample}_indels_summary.html"
    log = f"{log_dir}/{sample}_Indels_annotation.log"
    cmd = [
        "java", "-jar", snpEff,
        "-v", organism,
        "-noLog", "-no-downstream", "-no-upstream", "-no-utr",
        "-o", "vcf",
        Indels_vcf_input,
        "-csvStats", Indels_csv_stats,
        "-htmlStats", Indels_html_stats
    ]
    with open(log, "w") as logfile, open(Indels_vcf_output, "w") as Indels_vcf_out:
        process = subprocess.Popen(cmd, stdout=Indels_vcf_out, stderr=subprocess.PIPE, text=True)
        for line in process.stderr:
            print(line, end="")
            logfile.write(line)
        process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"Indels annotation failed for {sample}. Check {log}")
    print(f"Annotated Indels for sample {sample} written to {Indels_vcf_output}")
    print(f"Indels stats: {Indels_csv_stats}, HTML summary: {Indels_html_stats}")


#####################################################################################################
print(f"Script ended at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"Total script duration: {time.time() - script_start:.2f} seconds")