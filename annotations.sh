#!/bin/sh
echo "Annotations:"

# Inputs argumenets 
HPO=$1
input_vcf=$2
output=$3

#--------------------------------------------------------------------------------------------
# StrVCTVRE score

strVCTVRE_path="./annotation/"

python $strVCTVRE_path/StrVCTVRE.py -i $input_vcf -o $output_vcf
#--------------------------------------------------------------------------------------------

# GADO score

GADO_path="./annotation/GadoCommandline-1.0.1"

java -jar $GADO_path/GADO.jar \
 --mode PROCESS \
 --output hpoProcessed.txt \
 --caseHpo $HPO \
 --hpoOntology $GADO_path/hp.obo \
 --hpoPredictionsInfo $GADO_path/hpo_predictions_info_01_02_2018.txt

 java -jar $GADO_path/GADO.jar \
  --mode PRIORITIZE \
  --output $output \
  --caseHpoProcessed hpoProcessed.txt \
  --genes  $GADO_path/hpo_predictions_genes_01_02_2018.txt \
  --hpoPredictions $GADO_path/hpo_predictions_sigOnly_spiked_01_02_2018

rm hpoProcessed.txt
rm *.txtgado.log
rm $output/gado_resultgado.log
#--------------------------------------------------------------------------------------------

# AnnotSV annotations

export ANNOTSV=../annotation/AnnotSV
$ANNOTSV/bin/AnnotSV -SVinputFile $input_vcf  -genomeBuild GRCh38 -outputFile $output 

#--------------------------------------------------------------------------------------------

# Annovar annotations
path_annovar='./annotation/'
perl $path_annovar/table_annovar.pl  $input_vcf  $path_annovar/humandb/ -buildver hg38 -out $output -remove -protocol phastConsElements100way,genomicSuperDups -operation  r,r -nastring . -vcfinput
