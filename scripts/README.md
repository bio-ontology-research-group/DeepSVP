# DeepSVP - Data preprocessing and training

- Details for predicting gene-disease associations with DL2Vec can be found in the [experiment](https://github.com/bio-ontology-research--group/DL2Vec/tree/master/Experiment).
- Dataset for training: [dbVar:benign](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.clinical.benign_or_likely_benign.vcf.gz), and 
    [dbVar:pathogenic](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/GRCh38.variant_call.clinical.pathogenic_or_likely_pathogenic.vcf.gz)

- ``annotations.sh``: This script is used to annotate the varaints.
- ``data_preprocessing.py``: preprocessing the annotations and features.
- ``pheno_model.py``: script to get the DL2vec score using the trained model.
- ``deepsvp_training.py``: script to train and testing the model, with Hyperparameter optimization\
	In order to run use the comand: ```python deepsvp_training.py ONTOLOGY TYPE_OF_OPERATION```
	where ONTOLOGY is any of: ```(mp, go, cl, hp, uberon, union)``` and TYPE_OF_OPERATION: ```(max, mean)```
	
- ``BWA_GATK.sh`` : script to run GATK workflow for the input fastq files for the real samples, run using KAUST Supercomputing [IBEX](https://www.hpc.kaust.edu.sa/ibex).
- ``run_Manta.sh`` : script to generate VCF with the structural variants (SVs), we used [Manta](https://github.com/Illumina/manta) to identify the candidate SVs.  run using KAUST Supercomputing [IBEX](https://www.hpc.kaust.edu.sa/ibex).
