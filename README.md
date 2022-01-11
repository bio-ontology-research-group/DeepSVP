# DeepSVP
DeepSVP is a computational method to prioritize structural variants involved in genetic diseases by combining genomic information with information about gene functions. We incorporate phenotypes linked to genes, functions
  of gene products, gene expression in individual celltypes, and
  anatomical sites of expression, and systematically relate them to
  their phenotypic consequences through ontologies and machine
  learning
                                                                  
## Dataset
We train and evaluate our method using human genomic Structural Variation collected from [dbvar](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/) dataset.

## Prediction the candidate CNVs workflow
We integrate the annotates from Gene ontology [GO](http://geneontology.org/docs/download-go-annotations/), Uber-anatomy ontology
 [UBERON](https://www.ebi.ac.uk/ols/ontologies/uberon), Mammalian Phenotype ontology [MP](http://www.informatics.jax.org/vocab/mp_ontology), and Human Phenotype Ontology [HPO](https://hpo.jax.org/app/download/annotation) using [DL2vec](https://github.com/bio-ontology-research-group/DL2Vec). We convert different types of Description Logic axioms into graph representation, and then generate an embedding for each node and edge type.
We collected genomics features using public tool [AnnotSV (v2.3 or 2.2)](https://lbgi.fr/AnnotSV/annotations). 


## Installation 
```
pip install deepsvp
```

## Running the prediction model 
- Download all the files in [data](https://bio2vec.cbrc.kaust.edu.sa/data/DeepSVP/) and place them into data folder.
- Download and install the required database [AnnoSV (v2.3 or 2.2)](https://lbgi.fr/AnnotSV/downloads), and then run:
    ```
    bash scripts/annotation.sh -i input.vcf -o annotated_file
    ```
    and place the annotated VCF file into data folder. 

- Run the command `deepsvp --help` to display help and parameters:
    ```
    Usage: main.py [OPTIONS]
      
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

### Example:

    deepsvp -d data/ -i example_annotsv.tsv -p HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218 -m cl -maf 0.01 -ag max -o example_output.txt

or by using [cwl-runner](https://github.com/common-workflow-language/cwltool), modify the input file in the input example yaml [deepsvp.yaml](https://github.com/bio-ontology-research-group/DeepSVP/blob/master/deepsvp.yaml) file and then run:

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
