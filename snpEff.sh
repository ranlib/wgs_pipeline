#!/bin/bash

ORG=NC_000962.3

mkdir -p snpEff/data/$ORG
cp references/$ORG/$ORG.gb snpEff/data/$ORG/genes.gbk
cp /home/dieterbest/anaconda3/envs/wgs_pipeline/share/snpeff-5.1-2/snpEff.config snpEff/
echo "$ORG.genome: $ORG" >> snpEff/snpEff.config
snpEff build -c snpEff/snpEff.config -genbank -v $ORG

exit 0
