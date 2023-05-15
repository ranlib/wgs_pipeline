#!/usr/bin/env python
"""
WGS analysis pipeline
"""

import os
import argparse
import vcfpy
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Set up argument parser
parser = argparse.ArgumentParser(description='Whole genome analysis pipeline')
parser.add_argument('--fastq1', type=str, help='Path to first paired FASTQ file')
parser.add_argument('--fastq2', type=str, help='Path to secnd paired FASTQ file')
parser.add_argument('--sample', type=str, help='Sample name')
parser.add_argument('--genbank', type=str, help='Path to reference genome GenBank file')
parser.add_argument('--threads', type=int, default=4, help='Number of threads to use (default: 4)')
parser.add_argument('--output_dir', type=str, help='Path to output directory')
parser.add_argument('--contaminant-db', type=str, help='Path to contaminant database')
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok =True)

# Create FASTA file from GenBank file
fasta_file = os.path.splitext(args.genbank)[0] + ".fasta"
if not os.path.exists(fasta_file):
    with open(fasta_file, "w") as f:
        for record in SeqIO.parse(args.genbank, "genbank"):
            seq = SeqRecord(record.seq, record.id, "", "")
            SeqIO.write(seq, f, "fasta")

# Create BWA index for reference genome if it doesn't exist
if not os.path.exists(f"{fasta_file}.bwt"):
    os.system(f"bwa index {fasta_file}")

# Decontaminate FASTQ files with bbduk
output_base = os.path.basename(args.fastq1).replace(".fastq", "")
output_path = os.path.join(args.output_dir, "bbduk")
os.makedirs(output_path, exist_ok =True)
fq_clean_1 = os.path.join(output_path, args.sample + "_1.fq.gz")
fq_clean_2 = os.path.join(output_path, args.sample + "_2.fq.gz")
fq_clean_3 = os.path.join(output_path, args.sample + "_stats.txt")
os.system(f"bbduk.sh -Xmx6g in1={args.fastq1} in2={args.fastq2} out1={fq_clean_1} out2={fq_clean_2} ref={args.contaminant_db} k=31 hdist=1 stats={fq_clean_3}")

# Perform FastQC quality control
fastqc_dir = os.path.join(args.output_dir, "fastqc")
os.makedirs(fastqc_dir, exist_ok =True)
os.system(f"fastqc {fq_clean_1} {fq_clean_2} --outdir {fastqc_dir} --extract -t {args.threads}")

# Clean up reads with Trimmomatic
output_path = os.path.join(args.output_dir, "trimmomatic")
os.makedirs(output_path, exist_ok =True)
fq_trimd_ispaired_1 =  os.path.join(output_path, args.sample + "_ispaired_1.fq.gz")
fq_trimd_unpaired_1 =  os.path.join(output_path, args.sample + "_unpaired_1.fq.gz")
fq_trimd_ispaired_2 =  os.path.join(output_path, args.sample + "_ispaired_2.fq.gz")
fq_trimd_unpaired_2 =  os.path.join(output_path, args.sample + "_unpaired_2.fq.gz")
os.system(f"trimmomatic PE -threads {args.threads} {fq_clean_1} {fq_clean_2} {fq_trimd_ispaired_1} {fq_trimd_unpaired_1} {fq_trimd_ispaired_2} {fq_trimd_unpaired_2} -phred33 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50")

# Map reads with BWA
bwa_path = os.path.join(args.output_dir,"bwa")
os.makedirs(bwa_path, exist_ok =True)
bwa_output = os.path.join(bwa_path, args.sample + ".sam")
read_group=f"\"@RG\\tID:{args.sample}\\tSM:{args.sample}\\tPL:ILLUMINA\""
command = f"bwa mem -t {args.threads} -R {read_group} {fasta_file} {fq_trimd_ispaired_1} {fq_trimd_ispaired_2} > {bwa_output}"
os.system(command)

# Convert SAM to BAM and sort with SAMtools
samtools_path = os.path.join(args.output_dir, "samtools")
os.makedirs(samtools_path, exist_ok =True)
bam_path_unsorted = os.path.join(samtools_path, f"{args.sample}.unsorted.bam")
bam_path_issorted = os.path.join(samtools_path, f"{args.sample}.issorted.bam")
os.system(f"samtools view -@ {args.threads} -h -o {bam_path_unsorted} -Sb {bwa_output}")
os.system(f"samtools sort -@ {args.threads} -o {bam_path_issorted} {bam_path_unsorted}")

# Mark duplicates and perform SNP calling with GATK
# Run GATK Mutect2 in microbe mode
gatk_path = os.path.join(args.output_dir, "gatk")
os.makedirs(gatk_path, exist_ok =True)

gatk_dedup_bam = os.path.join(gatk_path, f"{args.sample}.dedup.bam")
gatk_metrics = os.path.join(gatk_path, f"{args.sample}.metrics.txt")
os.system(f"gatk MarkDuplicates -I {bam_path_issorted} -O {gatk_dedup_bam} -M {gatk_metrics} --REMOVE_DUPLICATES true")

# index all input files for gatk
os.system(f"samtools faidx {fasta_file}")
os.system(f"samtools index {gatk_dedup_bam}")
vcf_file = os.path.join(gatk_path, f"{args.sample}.vcf.gz")
os.system(f"gatk Mutect2 -R {fasta_file} -I {gatk_dedup_bam} -O {vcf_file} --max-mnp-distance 2")

vcf_file_split = os.path.join(gatk_path, f"{args.sample}_split.vcf.gz")
os.system(f"gatk LeftAlignAndTrimVariants -R {fasta_file} -V {vcf_file} -O {vcf_file_split} --split-multi-allelics")

vcf_file_filtered = os.path.join(gatk_path, f"{args.sample}_filtered.vcf.gz")
os.system(f"gatk FilterMutectCalls -R {fasta_file} -V {vcf_file} --min-reads-per-strand 1 --min-median-read-position 10 --min-allele-fraction 0.01 --microbial-mode true -O {vcf_file_filtered}")

# Convert vcf file to csv
csv_file = os.path.join(gatk_path, f"{args.sample}_filtered.csv.gz")
os.system(f"gatk VariantsToTable -V {vcf_file_filtered} -F CHROM -F POS -F REF -F ALT -GF GT -GF AD -GF DP -O {csv_file}")


# Get the gene IDs you want to filter for
gene_ids = [
"rpoB",
"katG",
"inhA",
"ahpC",
"rrs",
"rpsL",
"gidB",
"embB",
"embA",
"embC",
"pncA",
]

gb_record = SeqIO.parse(args.genbank, "genbank")

# Get the positions of the genes from the GenBank file
gene_positions = []
for feature in gb_record.features:
    if feature.type == "gene" and feature.qualifiers.get("gene_id") in gene_ids:
        start = feature.location.start.position
        end = feature.location.end.position
        gene_positions.append((start, end))

print("> gene_positions:")
print(gene_positions)
vcf_reader = vcfpy.Reader.from_path(vcf_file_filtered)

# Filter the VCF file for the positions of the genes
filtered_vcf_records = []
for record in vcf_reader:
    for pos in gene_positions:
        if record.POS >= pos[0] and record.POS <= pos[1]:
            filtered_vcf_records.append(record)

# Do something with the filtered VCF records
for record in filtered_vcf_records:
    print(record)
