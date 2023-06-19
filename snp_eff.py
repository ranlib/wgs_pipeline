#!/usr/bin/env python
"""
test script for snpEff
"""
import os
import argparse

parser = argparse.ArgumentParser(description="Annotate vcf file with snpEff.")
parser.add_argument("--input", "-i", type=str, help="Input vcf file")
parser.add_argument("--reference", "-r", type=str, help="Name of reference genome")
parser.add_argument("--database", "-d", type=str, help="Path of snpEff config file")
parser.add_argument("--output", "-o", type=str, help="Output vcf file")
args = parser.parse_args()

# annotate vcf with snpEff
database = args.database
genome = args.reference
command = f"snpEff -c {database} -csvStats stats.csv {genome} {args.input} > {args.output}"
# print(command)
os.system(command)
