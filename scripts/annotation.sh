#!/bin/sh
echo "AnnotSV Annotations:"

# Inputs argumenets 
input_vcf=$2
output=$4

# AnnotSV annotations
export ANNOTSV=data/AnnotSV
$ANNOTSV/bin/AnnotSV -SVinputFile $input_vcf  -genomeBuild GRCh38 -outputFile $output 
