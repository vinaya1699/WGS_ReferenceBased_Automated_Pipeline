**Reference-Based Whole Genome Sequencing (WGS) Automated Pipeline**

🚀 An end-to-end reference-based Whole Genome Sequencing (WGS) analysis pipeline written in Python.

This pipeline is not suitable for Human data as BQSR step is not included.

Author: Vinaya Kadam

**🔬 Pipeline Workflow**

<img width="1578" height="3846" alt="Workflow" src="https://github.com/user-attachments/assets/9397d584-f44b-46f0-9342-34770b60fed3" />

**📌 Features**

Raw data quality control using FastQC

Read trimming and filtering using fastp

Alignment to reference genome using BWA-MEM

BAM processing and duplicate marking using Samtools & Picard ; Calculation of genome coverage

Variant calling using GATK HaplotypeCaller

SNP & INDEL filtration following GATK Best Practices

Variant annotation using snpEff

Automatic logging for each step

Supports multiple samples in one run

**📂 Directory Structure**

The pipeline expects and generates the following structure:

**Working_Directory**

 0_Reference_Genome
   ├── organism.fasta
   └── organism.gtf

1_RawData
   ├── sample1_R1.fastq.gz
   ├── sample1_R2.fastq.gz
   └── Fastqc_Output

2_Clean_data
   └── Fastqc_Output

3_Alignment

4_Variant_Calling

5_Variant_Filtration

6_Variant_Annotation

logs

**🧬 Input Requirements**
1. Raw FASTQ files

Paired-end reads

Naming format:

<sample_name>_R1.fastq.gz

<sample_name>_R2.fastq.gz

2. Reference files

Located inside 0_Reference_Genome/:

Reference genome FASTA:
organism.fasta


Annotation file (GTF format):
organism.gtf


**⚠️ organism name must match the filename prefix exactly.**

**🛠️ Software Dependencies**

Make sure the following tools are installed and available in your $PATH:

Python ≥ 3.7
FastQC
fastp
BWA
Samtools
Picard
GATK (≥ 4.2)
snpEff (≥ 5.0)
Java (≥ 8)
Python libraries:
pip install pandas numpy

**▶️ Usage**

python WGS_Reference_Based.py -d /path/to/Working_Directory -org organism_name -t 20

Arguments
Argument	Description
-d / --Working_Directory	Working directory containing raw data
-org / --organism	Organism name (reference FASTA prefix)
-t / --threads	Number of threads (default: 4)




**📊 Output Files**

BAM files (duplicate-marked)

Raw & filtered VCF files

Annotated SNP & INDEL VCFs

snpEff HTML summary reports

Detailed log files for each step

**⚠️ Notes**

This pipeline follows GATK Best Practices for hard filtering.

Filtration thresholds can be modified inside the script if required.

Designed for Linux-based HPC or server environments.
