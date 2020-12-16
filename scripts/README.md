# DeepSVP - Data preprocessing and training

Details for predicting gene-disease associations with DL2Vec can be found in the [experiment](https://github.com/bio-ontology-research-group/DL2Vec/tree/master/Experiment).
- ``annotations.sh``: This script is used to annotate the varaints.
- ``data_preprocessing.py``: preprocessing the annotations and features.
- ``dl2vec_prediction.py``: script to get the DL2vec score using the trained model.
- ``model_training.py``: script to train and testing the model, with Hyperparameter optimization
- ``BWA_GATK.sh`` : script to run GATK workflow for the input fastq files for the real samples, run using KAUST Supercomputing [IBEX](https://www.hpc.kaust.edu.sa/ibex).
- ``run_Manta.sh`` : script to generate VCF with the structural variants (SVs), we used [Manta](https://github.com/Illumina/manta) to identify the candidate SVs.  run using KAUST Supercomputing [IBEX](https://www.hpc.kaust.edu.sa/ibex).
