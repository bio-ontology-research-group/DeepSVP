# PredCNV
PredCNV is a method that prioritizing Copy Number Variants (CNV) using Phenotype and Gene Functional Similarity. 

## Dataset
We train and evaluate our method using human genomic Structural Variation collected from [dbvar](https://ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_assembly/GRCh38/vcf/) dataset.

## Dependencies
The code was developed and tested using python 3.7. To install python dependencies run:  
 `pip install -r requirements.txt`

## Prediction the candidate CNVs workflow
We integrate the annotates from Gene ontology [GO](http://geneontology.org/docs/download-go-annotations/), Uber-anatomy ontology
 [UBERON](https://www.ebi.ac.uk/ols/ontologies/uberon), Mammalian Phenotype ontology [MP](http://www.informatics.jax.org/vocab/mp_ontology), and Human Phenotype Ontology [HPO](https://hpo.jax.org/app/download/annotation) using [DL2vec](https://github.com/bio-ontology-research-group/DL2Vec). We convert different types of Description Logic axioms into graph representation, and then generate an embedding for each node and edge type.
In addation we collected genomics features using public tools, [GADO](https://www.nature.com/articles/s41467-019-10649-4/), [StrVCTVRE](https://github.com/andrewSharo/StrVCTVRE), and databases [Annovar](https://annovar.openbioinformatics.org/), [AnnotSV](https://lbgi.fr/AnnotSV/annotations). 

## Scripts
- Details for predicting gene-disease associations with DL2Vec can be found in the [experiment](https://github.com/bio-ontology-research-group/DL2Vec/tree/master/Experiment).
- ``annotations.sh``: This script is used to annotate the varaints.
- ``data_preprocessing.py``: preprocessing the annotations and features selection. 
- ``dl2vec_score.py``: script to get the DL2vec score using the trained model.
- ``training.py``: script to train and testing the model, with Hyperparameter optimization

## Running PredCNV using pretrained models
1. Download the distribution file in [PredCNV.gz]()
2. Extract the distribution files PredCNV
3. Download the required database by: `cd PredCNV` then run:  `bash annotation/download.sh`
4. Run the command `python runPredCNV.py -h` to display help and parameters:
```
PredCNV: A phenotype-based tool to prioritize caustive CNV using WGS data and
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
```

### Example:
    python runPredCNV.py -inputfile example.vcf -hpo HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218,HP:0000007  -outfile example_output.txt -model "hp" -operation mean

 ```   
 Annotate VCF file (example.vcf) with the phenotypes (HP:0003701,HP:0001324,HP:0010628,HP:0003388,HP:0000774,HP:0002093,HP:0000508,HP:0000218,HP:0000007)...
 |========                        | 25% Annotated files generated successfully.
 |================                | 50% Phenotype prediction...
 |========================        | 75% CNV Prediction...
 |================================| 100%
The analysis is Done. You can find the priortize list in the output file: example_output.txt 
```
### Output:
The script will output a ranking a score for the candidate caustive CNV. 
``head example_output.txt``
```
+---------------------------------+-----------------+--------+------------+------------+-----------+-----------+
| ID                              |   PredCNV_Score |   rank | SV chrom   |   SV start |    SV end | SV type   |
|---------------------------------+-----------------+--------+------------+------------+-----------+-----------|
| nssv16214430                    |     1.00000e+00 |      1 | 1          |    1020163 |   2306775 | DEL       |
| DEL_pindel_34282                |     1.59618e-08 |      2 | 11         |   92050773 |  92050934 | DEL       |
| BI_GS_CNV_2_179516678_179528763 |     2.62666e-15 |      3 | 2          |  178651951 | 178664036 | DEL       |
| ALU_umary_ALU_666               |     6.02897e-25 |      4 | 1          |  188579454 | 188579455 | ALU       |
| ALU_umary_ALU_6671              |     4.44479e-25 |      5 | 8          |   34341082 |  34341083 | ALU       |
| ALU_umary_ALU_9350              |     9.27781e-26 |      6 | 12         |   76881950 |  76881951 | ALU       |
| ALU_umary_ALU_948               |     2.35165e-26 |      7 | 1          |  247589325 | 247589326 | ALU       |
| ALU_umary_ALU_12787             |     8.06472e-27 |      8 | X          |  127645122 | 127645123 | ALU       |
| ALU_umary_ALU_778               |     8.06472e-27 |      9 | 1          |  214046225 | 214046226 | ALU       |
+---------------------------------+-----------------+--------+------------+------------+-----------+-----------+

```

## Note
For any questions or comments please contact: azza.althagafi@kaust.edu.sa
