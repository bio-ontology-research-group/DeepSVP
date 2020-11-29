#!/bin/bash
#######################################################################################################################################################
#  Trimmomatic + GATK 4  + bgzip ( the *.g.VCF files) & tabix-index based workflow 
#    Azza Althagafi, CBRC
#	   Nagarajan Kathiresan, KSL 
#       Modification required for any users:
#        1. INPUT (directory location)
#        2. PROJECT (directory location)
#
########################################################################################################################################################



## 1. Software module 

module load trimmomatic/0.38 bwa/0.7.17/gnu-6.4.0 samtools/1.8 gatk/4.0.1.1 tabix/0.2.6


## 2. Reference list 

export REF=/ibex/reference/KSL/human_g1k_v37_decoy/human_g1k_v37_decoy.fasta
export MILLS=/ibex/reference/KSL/human_ref/Mills_and_1000G_gold_standard.indels.b37.vcf
export DB_SNP=/ibex/reference/KSL/human_ref/dbsnp_138.b37.vcf
export GATK=/sw/csi/gatk/3.8/el7.6_binary/GATK-3.8/GenomeAnalysisTK.jar


## Directory Variables 
export PROJECT=/ibex/scratch/kathirn/work/NIST_data/test_naga/Trio
export INPUT=/ibex/scratch/kathirn/work/NIST_data/test_naga/Trio/data

export BAM=${PROJECT}/BAM;
export VCF=${PROJECT}/VCF;
export LOGS=${PROJECT}/LOGS;
mkdir -p $BAM;
mkdir -p $VCF;
mkdir -p $LOGS;

## For Sample count 
set COUNT=0;

## Workflow steps 
export EXTN="clean.fq.gz";
for SAMPLE in `ls $INPUT/*_1.$EXTN`;
 do 
    PREFIX=`basename $SAMPLE _1.$EXTN` ;
    LOCATION=${SAMPLE%/*};

 #### Step 1. trimming of reads
    MEM="32gb"
    CORES=4
    JOB1_NAME="Trimming"
    JOB1_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB1_NAME}.${PREFIX} --time=2:00:00 --output=$LOGS/${JOB1_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB1_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB1_CMD="time -p java -XX:+UseParallelGC -XX:ParallelGCThreads=${CORES} -jar $TRIMMOMATIC_JAR PE -phred33 $LOCATION/${PREFIX}_1.$EXTN $LOCATION/${PREFIX}_2.$EXTN $BAM/$PREFIX.trimmed.P1.fastq $BAM/$PREFIX.up.1.fast $BAM/$PREFIX.trimmed.P2.fastq $BAM/$PREFIX.up.2.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50" ;
    JOB1_ID=$(${JOB1_TYPE} --parsable --wrap="${JOB1_CMD}");
    echo "$PREFIX sample with the job id=$JOB1_ID and Job Name=$JOB1_NAME submitted"
   
 #### Step 2. BWA MEM
    MEM="115gb"
    CORES=16
    JOB2_NAME="bwa-mem"
    JOB2_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB2_NAME}.${PREFIX} --time=12:00:00 --output=$LOGS/${JOB2_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB2_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB2_CMD="time -p bwa mem -M -k 30 -t $CORES $REF $BAM/$PREFIX.trimmed.P1.fastq $BAM/$PREFIX.trimmed.P2.fastq | samtools view -@ $CORES -b -S -h -q 30 - | samtools sort - > $BAM/$PREFIX.sorted.bam"
    JOB2_ID=$(${JOB2_TYPE} --parsable --dependency=afterok:${JOB1_ID} --wrap="${JOB2_CMD}");
    echo "$PREFIX sample with the job id=$JOB2_ID and Job Name=$JOB2_NAME submitted"   

 #### Step 3. MarkDuplicates 
    MEM="64gb"
    CORES=1
    JOB3_NAME="MarkDuplicate"
    JOB3_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB3_NAME}.${PREFIX} --time=4:00:00 --output=$LOGS/${JOB3_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB3_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB3_CMD="time -p gatk MarkDuplicates --INPUT=$BAM/$PREFIX.sorted.bam --METRICS_FILE=$BAM/$PREFIX.metrics.txt --OUTPUT=$BAM/$PREFIX.rmdup.bam"
    JOB3_ID=$(${JOB3_TYPE} --parsable --dependency=afterok:${JOB2_ID} --wrap="${JOB3_CMD}");
    echo "$PREFIX sample with the job id=$JOB3_ID and Job Name=$JOB3_NAME submitted"   

 ##### 4. AddOrReplace
    ## Note: VALIDATION_STRINGENCY=LENIENT missing in GATK 4.0 
    MEM="32gb"
    CORES=1
    JOB4_NAME="AddOrReplace"
    JOB4_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB4_NAME}.${PREFIX} --time=3:00:00 --output=$LOGS/${JOB4_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB4_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB4_CMD="time -p gatk AddOrReplaceReadGroups --INPUT=$BAM/$PREFIX.rmdup.bam --OUTPUT=$BAM/$PREFIX.rgroup.bam --SORT_ORDER=coordinate --RGSM=$PREFIX --RGPU=none --RGID=1 --RGLB=lib1 --RGPL=Illumina"
    JOB4_ID=$(${JOB4_TYPE} --parsable --dependency=afterok:${JOB3_ID} --wrap="${JOB4_CMD}");
    echo "$PREFIX sample with the job id=$JOB4_ID and Job Name=$JOB4_NAME submitted"   
    
 ##### 5. Samtools Index
    MEM="32gb"
    CORES=1
    JOB5_NAME="Samtool-Index"
    JOB5_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB5_NAME}.${PREFIX} --time=90:00 --output=$LOGS/${JOB5_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB5_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB5_CMD="time -p samtools index $BAM/$PREFIX.rgroup.bam"
    JOB5_ID=$(${JOB5_TYPE} --parsable --dependency=afterok:${JOB4_ID} --wrap="${JOB5_CMD}");
    echo "$PREFIX sample with the job id=$JOB5_ID and Job Name=$JOB5_NAME submitted"   


##### 6. Before BaseRecalibrator -  First pass of the Base Quality Score Recalibration (BQSR)
    MEM="115gb"
    CORES=16
    JOB6_NAME="BaseRecalibrator.before"
    JOB6_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB6_NAME}.${PREFIX} --time=18:00:00 --output=$LOGS/${JOB6_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB6_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB6_CMD="time -p gatk BaseRecalibrator --input $BAM/$PREFIX.rgroup.bam --known-sites $MILLS --reference $REF --known-sites $DB_SNP --output $BAM/$PREFIX.recal.before.table"
    JOB6_ID=$(${JOB6_TYPE} --parsable --dependency=afterok:${JOB5_ID} --wrap="${JOB6_CMD}");
    echo "$PREFIX sample with the job id=$JOB6_ID and Job Name=$JOB6_NAME submitted"  


##### 7. PrintReads  
    MEM="115gb"
    CORES=16
    JOB7_NAME="PrintReads"
    JOB7_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB7_NAME}.${PREFIX} --time=10:00:00 --output=$LOGS/${JOB7_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB7_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB7_CMD="time -p gatk PrintReads --input $BAM/$PREFIX.rgroup.bam --reference $REF --output $BAM/$PREFIX.recal.bam"
    JOB7_ID=$(${JOB7_TYPE} --parsable --dependency=afterok:${JOB5_ID} --wrap="${JOB7_CMD}");
    echo "$PREFIX sample with the job id=$JOB7_ID and Job Name=$JOB7_NAME submitted" 



##### 8. After BaseRecalibrator -  First pass of the Base Quality Score Recalibration (BQSR)
    MEM="115gb"
    CORES=16
    JOB8_NAME="BaseRecalibrator.after"
    JOB8_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB8_NAME}.${PREFIX} --time=18:00:00 --output=$LOGS/${JOB8_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB8_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB8_CMD="time -p gatk BaseRecalibrator --input $BAM/$PREFIX.recal.bam --known-sites $MILLS --reference $REF --known-sites $DB_SNP --output $BAM/$PREFIX.recal.after.table"
    JOB8_ID=$(${JOB8_TYPE} --parsable --dependency=afterok:${JOB7_ID} --wrap="${JOB8_CMD}");
    echo "$PREFIX sample with the job id=$JOB8_ID and Job Name=$JOB8_NAME submitted"  

 
 ##### 9. HaplotypeCaller
    MEM="115gb"
    CORES=16
    JOB9_NAME="HaplotypeCaller"
    JOB9_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB9_NAME}.${PREFIX} --time=30:00:00 --output=$LOGS/${JOB9_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB9_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB9_CMD="time -p gatk HaplotypeCaller --input $BAM/$PREFIX.recal.bam --output $VCF/$PREFIX.g.vcf --reference $REF --dbsnp $DB_SNP --emit-ref-confidence GVCF --max-alternate-alleles 16 --lenient true --native-pair-hmm-threads 16"
    JOB9_ID=$(${JOB9_TYPE} --parsable --dependency=afterok:${JOB7_ID} --wrap="${JOB9_CMD}");
    echo "$PREFIX sample with the job id=$JOB9_ID and Job Name=$JOB9_NAME submitted"   

 ##### 10. Compress the g.VCF file using bgzip
    MEM="32gb"
    CORES=1
    JOB10_NAME="bgzip"
    JOB10_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB10_NAME}.${PREFIX} --time=2:00:00 --output=$LOGS/${JOB10_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB10_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB10_CMD="time -p bgzip -f $VCF/$PREFIX.g.vcf" ;
    #JOB10_ID=$(${JOB10_TYPE} --parsable --dependency=afterok:${JOB9_ID} --wrap="${JOB10_CMD}");
    JOB10_ID=$(${JOB10_TYPE} --parsable --wrap="${JOB10_CMD}");
    echo "$PREFIX sample with the job id=$JOB10_ID and Job Name=$JOB10_NAME submitted"   

 ##### 11. Create Tabix-Index for the g.VCF.gz file
    MEM="32gb"
    CORES=1
    JOB11_NAME="tabix"
    JOB11_TYPE="sbatch -A ibex-cs --partition=batch --job-name=${JOB11_NAME}.${PREFIX} --time=30:00 --output=$LOGS/${JOB11_NAME}.${PREFIX}.%J.out --error=$LOGS/${JOB11_NAME}.${PREFIX}.%J.err --nodes=1 --cpus-per-task=${CORES} --mem=${MEM} --constraint=intel" ;
    JOB11_CMD="time -p tabix $VCF/$PREFIX.g.vcf.gz" ;
    JOB11_ID=$(${JOB11_TYPE} --parsable --dependency=afterok:${JOB10_ID} --wrap="${JOB11_CMD}");
    echo "$PREFIX sample with the job id=$JOB11_ID and Job Name=$JOB11_NAME submitted"   

 
 COUNT=$((COUNT + 1))
done
 
echo "========================================"
echo "Total Number of samples submitted: ${COUNT} in 11 Steps"
echo "========================================" 
