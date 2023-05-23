#!/bin/bash

time wgs_pipeline.py \
    --genbank references/NC_000962.3/NC_000962.3.gb \
    --fastq1 data/ERR552797_1.fastq.gz \
    --fastq2 data/ERR552797_2.fastq.gz \
    --sample ERR552797 \
    --threads 8 \
    --contaminant-db contamination/NC_001422.1.fa.gz \
    --output_dir ERR552797 \
    --genes TB_genes.txt \
    --intervals TB_genes_intervals.tsv \
    --vcf ERR552797.vcf

#run2> run.err

#time bam_qc.py --input ERR552797/samtools/ERR552797.issorted.bam --reference references/NC_000962.3/NC_000962.3.fasta --output ERR552797

#time multiqc ERR552797
