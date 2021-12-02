# DeepSVP: Evaluation
                                                   
## Dataset
The created synthetic genome for eveluating each method can be downloaded from [data]().

## Install and run different methods:

- Install the following methods by follow the insterction they provided and run all using the default setting:
  - [StrVCTVRE](https://github.com/andrewSharo/StrVCTVRE)
  - [CADD-SV](https://cadd-sv.bihealth.org/score)
  - [AnnotSV](https://github.com/lgmgeo/AnnotSV)
  
  
## Scripts 
- **evaluation.py**: This script is used to break the ties randommly and report different perfoermance measure.

  ``` 
  Usage: evaluation.py [OPTIONS]

  Options:
  -m, --method TEXT      Method to run
  -dir, --data-dir TEXT  Predictions directories
  -n, --num-p INTEGER    Number of patients
  --help                 Show this message and exit.
  
  example:
  python evaluation.py -m CADD-SV -n 1503  -dir user/data/ 
  ```
  
- **deepsvp_prediction.py**:  This script is used to run different DeepSVP models, and generate the input files to run using `evaluation.py`. 
 
   ```
   Usage: deepsvp_prediction.py [OPTIONS]

   Options:
   -m, --model_type TEXT     Ontology model, one of the following (go , mp , hp, cl, uberon, union)
   -ag, --aggregation TEXT   Aggregation method for the genes within CNV (max , mean)
   -sub, --subsets TEXT   Phenotypes subsets (full)

   example:
   python deepsvp_prediction.py -m cl -ag max -sub full
   ```

## Final notes
For any questions or comments please contact: azza.althagafi@kaust.edu.sa

