#!/usr/bin/env python
import csv
import vcfpy
from Bio import SeqIO

# Read in the GenBank file
gb_file = "references/NC_00932.3/NC_00932.3.gb"
gb_record = SeqIO.parse(gb_file, "genbank")

# Get the gene IDs you want to filter for
gene_names = [
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
# Read in the list of gene names to filter for from a text file
gene_names_file = "genes.txt"
with open(gene_names_file, "r") as f:
    gene_names = [line.strip() for line in f.readlines()]

# Get the positions of the genes from the GenBank file
gene_positions = {}
for record in gb_record:
    for feature in record.features:
        if feature.type == "gene":
            genes = feature.qualifiers.get("gene")
            if genes is not None:
                for gene in genes:
                    if gene in gene_names:
                        start = feature.location.start.position 
                        end = feature.location.end.position
                        gene_positions[gene] = (start, end)

# Write the gene positions to a CSV file
csv_output_file = "gene_positions.csv"
with open(csv_output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["gene_name", "start_pos", "end_pos"])
    for gene_name, positions in gene_positions.items():
        writer.writerow([gene_name, positions[0], positions[1]])

# Read in the VCF file
vcf_file = "ERR552797/gatk/ERR552797_filtered.vcf.gz"
vcf_reader = vcfpy.Reader.from_path(vcf_file)

# Filter the VCF file for the positions of the genes
filtered_vcf_records = []
for record in vcf_reader:
    for gene_name, positions in gene_positions.items():
        if record.POS >= positions[0] and record.POS <= positions[1]:
            filtered_vcf_records.append(record)

# Write the filtered VCF records to a new VCF file
vcf_output_file = "filtered.vcf"
vcf_writer = vcfpy.Writer.from_path(vcf_output_file, vcf_reader.header)
for record in filtered_vcf_records:
    vcf_writer.write_record(record)

# Close the VCF writer
vcf_writer.close()

# Do something with the filtered VCF records
#for record in filtered_vcf_records:
#    print(record)
