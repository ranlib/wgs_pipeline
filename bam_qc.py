#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser(description='Perform WGS QC with Picard.')
parser.add_argument('--input', type=str, help='Path to input BAM file')
parser.add_argument('--reference', type=str, help='Path to reference genome fasta file')
parser.add_argument('--output', type=str, help='Path to output directory')
args = parser.parse_args()

# Make output directory if it doesn't exist
if not os.path.exists(args.output):
    os.makedirs(args.output)

# Define paths for intermediate files
markdups_bam = os.path.join(args.output, 'markdups.bam')
metrics_file = os.path.join(args.output, 'metrics.txt')

# Step 1: Mark duplicates
os.system(f"gatk MarkDuplicates I={args.input} O={markdups_bam} M={metrics_file}")

# Step 2: Calculate coverage statistics
cov_output = os.path.join(args.output, 'coverage')
os.system(f"gatk CollectWgsMetrics I={markdups_bam} O={cov_output} R={args.reference}")

# # Step 3: Generate GC bias metrics
# gc_output = os.path.join(args.output, 'gc_bias')
# os.system(f"gatk CollectGcBiasMetrics I={markdups_bam} O={gc_output} R={args.reference} CHART={gc_output}.pdf S={gc_output}.summary")

# # Step 4: Generate quality control metrics
# qc_output = os.path.join(args.output, 'qc')
# os.system(f"gatk CollectMultipleMetrics I={markdups_bam} O={qc_output} R={args.reference}")

# # Step 5: Generate insert size metrics
# is_output = os.path.join(args.output, 'insert_size')
# os.system(f"gatk CollectInsertSizeMetrics I={markdups_bam} O={is_output} HISTOGRAM_FILE={is_output}.pdf")
