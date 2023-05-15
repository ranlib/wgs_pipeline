wgs_pipeline.py \
    --genbank references/NC_00932.3/NC_00932.3.gb \
    --fastq1 data/ERR552797_1.fastq.gz \
    --fastq2 data/ERR552797_2.fastq.gz \
    --sample ERR552797 \
    --threads 8 \
    --contaminant-db contamination/NC_001422.1.fa.gz \
    --output_dir ERR552797
