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

# Step 1: Mark duplicates
output_mark_duplicates  = os.path.join( args.output, "output_mark_duplicates")
os.makedirs( output_mark_duplicates, exist_ok = True) 
markdups_bam = os.path.join(output_mark_duplicates, 'markdups.bam')
metrics_file = os.path.join(output_mark_duplicates, 'metrics.txt')
os.system(f"gatk MarkDuplicates -I {args.input} -O {markdups_bam} -M {metrics_file}")

# Step 2: Calculate coverage statistics
output_collect_wgs_metrics = os.path.join(args.output, 'output_collect_wgs_metrics')
os.makedirs( output_collect_wgs_metrics, exist_ok = True)
cov_output = os.path.join(output_collect_wgs_metrics, 'coverage.txt')
os.system(f"gatk CollectWgsMetrics -I {markdups_bam} -O {cov_output} -R {args.reference}")

# Step 3: Generate GC bias metrics
output_collect_gc_bias_metrics = os.path.join(args.output, 'output_collect_gc_bias_metrics')
os.makedirs( output_collect_gc_bias_metrics, exist_ok = True)
gc_output = os.path.join(output_collect_gc_bias_metrics, 'gc_output.txt')
gc_output_chart = os.path.join(output_collect_gc_bias_metrics, 'gc_output_chart.pdf')
gc_output_summary = os.path.join(output_collect_gc_bias_metrics, 'gc_output.summary')
os.system(f"gatk CollectGcBiasMetrics -I {markdups_bam} -O {gc_output} -R {args.reference} -CHART {gc_output_chart} -S {gc_output_summary}")

# STEP 4: Generate quality control metrics
output_collect_multiple_metrics = os.path.join(args.output, 'output_collect_multiple_metrics')
os.makedirs( output_collect_multiple_metrics, exist_ok = True)
qc_output = os.path.join(output_collect_multiple_metrics, 'qc') # qc is the suffix
os.system(f"gatk CollectMultipleMetrics -I {markdups_bam} -O {qc_output} -R {args.reference}")

# Step 5: Generate insert size metrics
output_collect_insert_size_metrics = os.path.join(args.output, 'output_collect_insert_size_metrics')
os.makedirs( output_collect_insert_size_metrics, exist_ok = True)
is_output = os.path.join(output_collect_insert_size_metrics, 'insert_size.pdf')
os.system(f"gatk CollectInsertSizeMetrics -I {markdups_bam} -O {is_output} -H {is_output}")
