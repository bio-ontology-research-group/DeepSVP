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
Usage: phenocnv [OPTIONS]

  PredCNV: A phenotype-based tool to prioritize caustive  CNV using WGS data
  and Phenotype/Gene Functional Similarity

Options:
  -d, --data-root TEXT     Data root folder  [required]
  -i, --in-file TEXT       Annotated Input file  [required]
  -p, --hpo TEXT           List of phenotype ids separated by commas
                           [required]

  -m, --model TEXT         Ontology model, one of the following (go , mp , hp,
                           cl, uberon, combined), default=mp

  -ag, --aggregation TEXT  Aggregation method for the genes within CNV (max or
                           mean) default=max

  -o, --outfile TEXT       Output result file
  --help                   Show this message and exit.
```

### Example:
    phenocnv -d data/ -i example.tsv -p HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218,HP:0000007  -m 'go' -ag 'max' -oo example_output.txt

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
