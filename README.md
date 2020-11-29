# PhenoCNV - Prioritizing Copy Number Variants (CNV) using Phenotype and Gene Functional Similarity
                                                                  
## Dataset
We train and evaluate our method using human genomic Structural Variation collected from [dbvar](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/) dataset.

## Prediction the candidate CNVs workflow
We integrate the annotates from Gene ontology [GO](http://geneontology.org/docs/download-go-annotations/), Uber-anatomy ontology
 [UBERON](https://www.ebi.ac.uk/ols/ontologies/uberon), Mammalian Phenotype ontology [MP](http://www.informatics.jax.org/vocab/mp_ontology), and Human Phenotype Ontology [HPO](https://hpo.jax.org/app/download/annotation) using [DL2vec](https://github.com/bio-ontology-research-group/DL2Vec). We convert different types of Description Logic axioms into graph representation, and then generate an embedding for each node and edge type.
We collected genomics features using public tool [AnnotSV](https://lbgi.fr/AnnotSV/annotations). 

![workflow](images/workflow.png)

## Installation 
```
pip install phenocnv
```

## Running the prediction model
1. Download the required database AnnoSV, then run `bash annotation/download.sh` and place the annotated file into data folder. 
2. Run the command `phenocnv --help` to display help and parameters:
```
PhenoCNV: A phenotype-based tool to prioritize caustive CNV using WGS data and
Phenotype/Gene Functional Similarity

optional arguments:
  -h, --help            show this help message and exit
  -inputfile [INPUTFILE]
                        Path to VCF file
  -hpo [HPO]            List of phenotype ids separated by commas
  -outfile [OUTFILE]    Path to results file
  -model [MODEL]        Preferred model (go, mp, uberon, hp,
                        go_ppi,mp_ppi,uberon_ppi,hp_ppi) , default='hp'
  -operation [OPERATION]
                        Preferred operation for the gene annotation in big CNV
                        regions (max, mean) , default='max'
  --data                path to data folder
```

### Example:
    phenocnv -inputfile example.vcf -hpo HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218,HP:0000007  -outfile example_output.txt -model "hp" -operation mean --data ./data/

 ```   
 Annotate VCF file (example.vcf) with the phenotypes (HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218,HP:0000007)...
 |========                        | 25% Annotated files generated successfully.
 |================                | 50% Phenotype prediction...
 |========================        | 75% CNV Prediction...
 |================================| 100%
The analysis is Done. You can find the priortize list in the output file: example_output.txt 
```
#### Output:
The script will output a ranking a score for the candidate caustive CNV. 


## Scripts
- Details for predicting gene-disease associations with DL2Vec can be found in the [experiment](https://github.com/bio-ontology-research-group/DL2Vec/tree/master/Experiment).
- ``annotations.sh``: This script is used to annotate the varaints.
- ``data_preprocessing.py``: preprocessing the annotations and features selection. 
- ``dl2vec_score.py``: script to get the DL2vec score using the trained model.
- ``training.py``: script to train and testing the model, with Hyperparameter optimization

## Final notes
For any questions or comments please contact: azza.althagafi@kaust.edu.sa
