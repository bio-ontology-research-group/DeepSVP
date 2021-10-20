#!/bin/sh
echo "AnnotSV Annotations:"

# Inputs argumenets 
input_vcf=$1
output=$2

# AnnotSV annotations
export ANNOTSV=data/AnnotSV
$ANNOTSV/bin/AnnotSV -SVinputFile $input_vcf  -genomeBuild GRCh38 -outputFile $output 
