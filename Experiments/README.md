# DeepSVP: Evaluation
                                                   
## Dataset
The created synthetic genomes for eveluating each method can be downloaded from [here](https://bio2vec.cbrc.kaust.edu.sa/data/DeepSVP/experiments.zip).

## Install and run different methods:

- Install the following methods using the instructions they provided and run all using the default setting:
  - [StrVCTVRE](https://github.com/andrewSharo/StrVCTVRE)
  - [CADD-SV](https://cadd-sv.bihealth.org/score)
  - [AnnotSV](https://github.com/lgmgeo/AnnotSV)
  
  
## Scripts 
- **evaluation.py**: This script is used to compute metrics (recall@rank, AUC, AUPR); ties are broken randomly.

  ``` 
  Usage: evaluation.py [OPTIONS]

  Options:
  -m, --method TEXT      Method to evaluate
  -dir, --data-dir TEXT  Predictions directories
  -n, --num-p INTEGER    Number of patients
  --help                 Show this message and exit.
  
  example:
  python evaluation.py -m CADD-SV -n 1503 -dir './data/CADDSV/*.tsv'
  python evaluation.py -m AnnotSV -n 1503 -dir './data/AnnotSV/*.tsv'
  python evaluation.py -m StrVCTVRE -n 1503 -dir './data/StrVCTVRE/*.tsv'
  python evaluation.py -m DeepSVP_hp_max -n 1503 -dir './data/DeepSVP_hp_max/*full.tsv'
  ```
  
- **deepsvp_prediction.py**:  This script is used to run different DeepSVP models, and generate the input files to run using `evaluation.py`. 
 
   ```
   Usage: deepsvp_prediction.py [OPTIONS]

   Options:
   -m, --onto TEXT     Ontology model, one of the following (go , mp , hp, cl, uberon, union)
   -ag, --ope TEXT   Aggregation operation for the genes within CNV (max , mean)
   -sub, --subsets TEXT   Phenotypes subsets (full)

   example:
   python deepsvp_prediction.py -m cl -ag max -sub full
   ```

## Final notes
For any questions or comments please contact: azza.althagafi@kaust.edu.sa


