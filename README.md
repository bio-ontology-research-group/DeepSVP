# DeepSVP
DeepSVP is a computational method to prioritize structural variants (SV) involved in genetic diseases by combining genomic information with information about gene functions. We incorporate phenotypes linked to genes, functions of gene products, gene expression in individual celltypes, and anatomical sites of expression. DeepSVP systematically relates them to their phenotypic consequences through ontologies and machine learning.
                                                                  
## Training dataset
We train and evaluate our method using human SV collected from [dbvar](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/) dataset.

## Annotation data sources (integrated in the candidate SV prediction workflow)
We integrated the annotations from different sources:
- Gene ontology ([GO](http://geneontology.org/docs/download-go-annotations/))
- Uber-anatomy ontology ([UBERON](https://www.ebi.ac.uk/ols/ontologies/uberon))
- Mammalian Phenotype ontology ([MP](http://www.informatics.jax.org/vocab/mp_ontology))
- Human Phenotype Ontology ([HPO](https://hpo.jax.org/app/download/annotation))

This work is done using [DL2vec](https://github.com/bio-ontology-research-group/DL2Vec). We convert different types of Description Logic axioms into graph representation, and then generate an embedding for each node and edge type.

We collected [genomics features](https://lbgi.fr/AnnotSV/annotations) using the [AnnotSV (v2.2)](https://lbgi.fr/AnnotSV/downloads) public tool. 


## Installation 
Using pip version 20.3.1:
```
pip install deepsvp
```

Or you can create a specific Conda Environments (e.g. named "deepsvp-py38-pip2031"):
```
conda create -n deepsvp-py38-pip2031 python=3.8 pip=20.3.1
conda activate deepsvp-py38-pip2031
pip3 install deepsvp
pip3 install networkx
pip3 install torch
pip3 list
conda deactivate
```

## Running the DeepSVP prediction model 
- Download all the files from [data](https://bio2vec.cbrc.kaust.edu.sa/data/DeepSVP/) and place them in the folder named "data":
```
mkdir DeepSVP/          ;# /path_of_your_DeepSVP_repository/
cd DeepSVP
wget "https://bio2vec.cbrc.kaust.edu.sa/data/DeepSVP/data.zip"
unzip data.zip
cd data                 ;# /path_of_your_DeepSVP_data_repository/
wget "https://bio2vec.cbrc.kaust.edu.sa/data/DeepSVP/experiments.zip"   # can be very long
unzip experiments.zip
```
- Download and install the required [AnnoSV (2.2)](https://lbgi.fr/AnnotSV/downloads) tool in the "data" folder:
```
cd /path_of_your_DeepSVP_data_repository/
wget "https://lbgi.fr/AnnotSV/Sources/AnnotSV_2.2.tar.gz"
gunzip -c AnnotSV_2.2.tar.gz | tar -xvf -
cd AnnotSV_2.2
make PREFIX=. install
cd ..
mv AnnotSV_2.2/ AnnotSV/
```

- Add AnnotSV (v2.2) annotation in your VCF input file ($your_input_vcf):
```
bash 
export ANNOTSV=/path_of_your_DeepSVP_data_repository/AnnotSV
$ANNOTSV/bin/AnnotSV/AnnotSV.tcl -SVinputFile $your_input_vcf  -genomeBuild GRCh38 -outputFile $your_annotsv_output.annotated.tsv

```
And place the annotated VCF file ($your_annotsv_output.annotated.tsv) in the data folder (/path_of_your_DeepSVP_data_repository/). 

- Run the command `deepsvp --help` to display help and parameters:
```
Usage: deepsvp [OPTIONS]
      
     DeepSVP: A phenotype-based tool to prioritize caustive CNV using WGS data
     and Phenotype/Gene Functional Similarity
  
Options:
    -d, --data-root TEXT      Data root folder  [required]
    -i, --in-file TEXT        Annotated Input file  [required]
    -p, --hpo TEXT            List of phenotype ids separated by commas
                              [required]
    -maf, --maf_filter FLOAT  Allele frequency filter using gnomAD and 1000G
                              default<=0.01
    -m, --model_type TEXT     Ontology model, one of the following (go , mp ,
                              hp, cl, uberon, union), default=mp
    -ag, --aggregation TEXT   Aggregation method for the genes within CNV (max
                              or mean) default=max
    -o, --outfile TEXT        Output result file
    --help                    Show this message and exit.        
```

- Run the example (with you own HPO terms):
```
    deepsvp -d data/ -i $your_annotsv_output.annotated.tsv -p HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218 -m cl -maf 0.01 -ag max -o example_output.txt
```    
Or run the example with the deepsvp-py38-pip2031 Conda Environment:
```
conda activate deepsvp-py38-pip2031
deepsvp -d data/ -i $your_annotsv_output.annotated.tsv -p HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218 -m cl -maf 0.01 -ag max -o example_output.txt
conda deactivate
```
Or by using [cwl-runner](https://github.com/common-workflow-language/cwltool), modify the input file in the input example yaml [deepsvp.yaml](https://github.com/bio-ontology-research-group/DeepSVP/blob/master/deepsvp.yaml) file and then run:

	cwl-runner deepsvp.cwl deepsvp.yaml 
    
 ```   
 |========                        | 25% Reading the input phenotypes...
 |================                | 50% Phenotype prediction... 
 |========================        | 75% CNV Prediction... 
 |================================| 100% DONE! You can find the prediction results in the output file: example_output.txt
```



#### Output:
The script will output a ranking a score for the candidate caustive CNV. 

## Scripts 
- Details for predicting pathogenic variants and comparison with other methods can be found in the [experiment](https://github.com/bio-ontology-research-group/DL2Vec/tree/master/Experiment) folder.
- ``annotations.sh``: This script is used to annotate the varaints.
- ``data_preprocessing.py``: preprocessing the annotations and features.
- ``pheno_model.py``: script to get the DL2vec score using the trained model.
- ``deepsvp_training.py``: script to train and testing the model, with Hyperparameter optimization
- ``BWA_GATK.sh`` : script to run GATK workflow for the input fastq files for the real samples, run using KAUST Supercomputing [IBEX](https://www.hpc.kaust.edu.sa/ibex).
- ``run_Manta.sh`` : script to generate VCF with the structural variants (SVs), we used [Manta](https://github.com/Illumina/manta) to identify the candidate SVs.  run using KAUST Supercomputing [IBEX](https://www.hpc.kaust.edu.sa/ibex).

## Final notes
For any questions or comments please contact: azza.althagafi@kaust.edu.sa
