from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso
import matplotlib.pyplot as plt
import time
from matplotlib import offsetbox
import matplotlib.patheffects as PathEffects
import pathlib
import pandas as pd
from pandas.plotting import scatter_matrix
import pandas.plotting
from pandas.plotting import parallel_coordinates
import numpy as np
import seaborn as sns  # for intractve graphs
from sklearn.utils import class_weight
from sklearn.preprocessing import StandardScaler
import pandas as pd
import gensim
import numpy as np
np.random.seed(0) 
import re
import pickle as pkl
from numpy import array
from numpy import argmax
import mygene
import gzip
import sys
import warnings
warnings.simplefilter(action='ignore')
import random
from sklearn.model_selection import train_test_split
import os 


# Constants and helper functions
def tidy_split(df, column, sep=',', keep=False):
    """
    Split the values of a column and expand so the new DataFrame has one split
    value per row. Filters rows where the column is missing.

    Params
    ------
    df : pandas.DataFrame
        dataframe with the column to split and expand
    column : str
        the column to split and expand
    sep : str
        the string used to split the column's values
    keep : bool
        whether to retain the presplit value as it's own row

    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the same columns as `df`.
    """
    indexes = list()
    new_values = list()
    for i, presplit in enumerate(df[column].astype(str)):
        values = presplit.split(sep)
        if keep and len(values) > 1:
            indexes.append(i)
            new_values.append(presplit)
        for value in values:
            indexes.append(i)
            new_values.append(value)
    new_df = df.iloc[indexes, :].copy()
    new_df[column] = new_values
    return new_df

def get_HPO(disease):
    """
    generate the human disease and hp association
    
    Params
    ------
    disease
    
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the OMIM and HP annotation
    """
    all_dis_phe=dict()
    my_dis={}
    with open("data/phenotype_annotation.tab","r",encoding='utf-8') as f:
        for line in f.readlines():
            data=line.split("\t")
            if (data[0]=="OMIM")&(data[5][:4]=="OMIM"):
                d = '<http://purl.obolibrary.org/obo/'+str(data[4].strip().replace(":","_"))+'>'
                if data[5].strip() not in all_dis_phe:
                    all_dis_phe[data[5].strip()]=[d]
                else:
                    all_dis_phe[data[5].strip()].append(d)
 
    with open("data/phenotype_to_genes.txt","r",encoding='utf-8') as f:
        for line in f.readlines():
            data=line.split("\t")
            if (data[5]=="mim2gene"):
                d = '<http://purl.obolibrary.org/obo/'+str(data[0].strip().replace(":","_"))+'>'
                if data[6].strip() not in all_dis_phe:
                    all_dis_phe[data[6].strip()]=[d]
                else:
                    all_dis_phe[data[6].strip()].append(d)
          
    for d in disease:
        if d in all_dis_phe:
            my_dis[d] = all_dis_phe[d]
            
    return my_dis, all_dis_phe

def GCcontent(data, types):
    data[types] = data[types].astype(str)
    data[types] = data[types].str.split('(').str.get(0).astype(float)
    return data

def preprocessing_numeric(features,cols,mean_std_train):
    for column in cols :
        if not column.startswith('HI_CGscore') and not column.startswith('TriS_CGscore'):
            features[column+"_indicator"] = np.where(features[column].isnull(), 1, 0)
            features[column] = (features[column]-mean_std_train[column]['mean'])/mean_std_train[column]['std']
            features[column].fillna(0, inplace=True)    
    return features

def preprocessing_categorical(features,cat):  
    features[cat] = features[cat].fillna('undefined') 
    values= [3.0,2.0,1.0,0.0,40.0,30.0,'undefined']
    features[cat] = features[cat].astype(pd.CategoricalDtype(values))
    features = pd.concat([features,pd.get_dummies(features[cat], prefix=cat)],axis=1)
    features.drop([cat],axis=1, inplace=True)
    return features

def merge_fetaures(cnv,genes,onto,operation):
    global mean_std_train
    features = cnv.merge(genes, on=['ID']) 
    
    categorical=['HI_CGscore','TriS_CGscore']
    for i in categorical:
        features = preprocessing_categorical(features,i)
        
    # split the data into training and testing 
    train = features[features['OMIM'].isin(dis_train)]
    test = features[features['OMIM'].isin(dis_test)]
    
    # save the mean and std to use it for the testing  
    col = list(train.loc[:,'SV length':].columns)
    for c in col:
        if c not in mean_std_train and not c.startswith('HI_CGscore') and not c.startswith('TriS_CGscore'):
            mean_std_train[c]={'mean':train[train[c].notna()][c].astype(float).mean(), 'std':train[train[c].notna()][c].astype(float).std()}   

    train = preprocessing_numeric(train,col,mean_std_train)
    test = preprocessing_numeric(test,col,mean_std_train)
            
    train.to_csv("data/"+onto+operation+"training.tsv", sep='\t', index=False)
    test.to_csv("data/"+onto+operation+"testing.tsv", sep='\t', index=False)
    
    return train, test


def get_scores(dl2vec_scores,dis, gene):
    g=str(gene).strip('.0')
    if dis in dl2vec_scores:
        if g in dl2vec_scores[dis]:
            results=dl2vec_scores[dis][g]
        else:
            results=np.nan        
        return results


def collect_features(data, onto):    
    # 1. Gene Features 
    #---------------------
    # get DL2vec scores
    Gene_features = data[data['AnnotSV type']=='split']
    if (onto == 'all'):
        onto_types = ["go", "mp", "uberon", "hp", "cl"]
        for i in onto_types:     
            with open("data/"+i+"_quantile_dl2vec_ranks.pkl","rb") as f: 
                dl2vec=pkl.load(f)
            Gene_features[i+'_score'] = Gene_features.apply(lambda x: get_scores(dl2vec,x['OMIM'], x['entrezgene']), axis=1)
        GF = ['ID','HI_DDDpercent','go_score','mp_score','cl_score','hp_score','uberon_score','CDS length','tx length','synZ_ExAC','misZ_ExAC','pLI_ExAC', 'delZ_ExAC', 'dupZ_ExAC', 'cnvZ_ExAC'] #'

    else:
        with open("data/"+onto+"_quantile_dl2vec_ranks.pkl","rb") as f: 
            dl2vec=pkl.load(f)
        Gene_features[onto+'_score'] = Gene_features.apply(lambda x:  get_scores(dl2vec,x['OMIM'], x['entrezgene']), axis=1)
        GF = ['ID','HI_DDDpercent',onto+'_score','CDS length','tx length','synZ_ExAC','misZ_ExAC', 'pLI_ExAC', 'delZ_ExAC', 'dupZ_ExAC', 'cnvZ_ExAC'] 
    Gene_features = Gene_features[GF].copy()
    Gene_features[GF[1:]]=Gene_features[GF[1:]].astype(float)
    
    # Gene operations: max and mean
    gene_features_max = Gene_features.groupby('ID')[GF[1:]].max().reset_index()
    gene_features_mean = Gene_features.groupby('ID')[GF[1:]].mean().reset_index()
   
    #2. CNV_features
    #---------------------
    CF=['SV chrom', 'SV start', 'SV end','label','OMIM','ID','SV type','SV length','GCcontent_right', 'GCcontent_left','NumGenes','NumPromoters','HI_CGscore','TriS_CGscore']#'StrVCTVRE'] #'SV type',','phastConsElements100way','genomicSuperDups']
    cnv_features = data[data['AnnotSV type']=='full']
    cnv_features['NumGenes'] = cnv_features['Gene name'].apply(lambda x: len(str(x).split('/'))).astype(int)
    cnv_features['NumPromoters'] = cnv_features['promoters'].apply(lambda x: len(str(x).split('/'))).astype(int)
    cnv_features = GCcontent(cnv_features, 'GCcontent_right')
    cnv_features = GCcontent(cnv_features, 'GCcontent_left')
    cnv_features['TriS_CGscore'] = cnv_features['TriS_CGscore'].replace('Not yet evaluated',np.NaN)
    cnv_features['SV length'] = cnv_features['SV length'].abs()
    cnv_features['SV type'] = pd.Categorical(cnv_features['SV type'])
    cnv_features['SV type'] =  cnv_features['SV type'].cat.codes
    cnv_features.drop_duplicates(inplace=True)
    cnv_features = cnv_features[CF].copy()
    cnv_features[CF[8:]]=cnv_features[CF[8:]].astype(float)

    # Merge CNV + Gene features and preprocessing
    #--------------------------------------------
    features_max_train, features_max_test = merge_fetaures(cnv_features,gene_features_max,onto,'max')
    features_avg_train, features_avg_test = merge_fetaures(cnv_features,gene_features_mean,onto,'mean')
     
    return features_max_train, features_max_test, features_avg_train, features_avg_test

def get_features(onto):
    print("Generate a training data for: ", onto)
    features = all_data.copy()
    train_max, test_max, train_avg, test_avg = collect_features(features,onto)    
    print('Done.')
    
if __name__ == '__main__':
    
    dbvar_benign='data/GRCh38.variant_call.clinical.benign_or_likely_benign.vcf'
    dbvar_pathogenic='data/GRCh38.variant_call.clinical.pathogenic_or_likely_pathogenic.vcf'
    dbvar_benign_annotsv='data/dbvar_benign.tsv'
    dbvar_pathogenic_annotsv='data/pathogenic_annotsv.tsv'
    mean_std_train={}

    Benign = pd.read_csv(dbvar_benign, skiprows=34, sep='\t',low_memory=False)
    Pathogenic = pd.read_csv(dbvar_pathogenic, skiprows=34, sep='\t',low_memory=False)
    Benign = Benign[Benign['ALT'].isin(['<DEL>','<DUP>'])] 
    Pathogenic = Pathogenic[Pathogenic['ALT'].isin(['<DEL>','<DUP>'])] 
    print('Number of Benign and Pathogenic varaints, Pathogenic: {}, Benign: {}'.format(Pathogenic.shape[0], Benign.shape[0]))
    Pathogenic['OMIM'] = Pathogenic['INFO'].str.findall(r'OMIM:\d{6}')

    # Pathogenic not caustive
    Pathogenic['OMIM'] = Pathogenic['OMIM'].apply(lambda y: np.nan if len(y)==0 else y)
    Pathogenic['label'] = np.where((Pathogenic['OMIM'].notna()), 1, 0)
    Benign['label'] = 0
    Benign['OMIM'] = np.nan

    #split the varaints with multiple OMIM to a separate variant
    Pathogenic = tidy_split(Pathogenic, "OMIM", ',', keep=False)
    Pathogenic['OMIM'] = Pathogenic['OMIM'].str.replace(r"[^:0-9a-zA-Z ]+", "").str.strip().replace('nan',np.nan)
    Pathogenic.drop_duplicates(inplace=True)

    # Map all the OMIM to HPO 
    dis = list(Pathogenic[Pathogenic['OMIM'].notna()]['OMIM'].drop_duplicates()) 
    dis_HPO , all_dis_phe = get_HPO(dis)
    print('Number of varaints with no phenotypes = {}'.format(Pathogenic[Pathogenic['OMIM'].isna()].shape[0])) 
    print('Number of OMIM with phenotypes: {}'.format(len(dis_HPO)))
            
    dis_HPO_list = list(dis_HPO.keys()) 

    all_data = pd.concat([Pathogenic,Benign])

    # pathoginc but not causing the diseases
    nopheno = all_data[(all_data['label']==0)]
    np.random.seed(2)
    random_pheno = np.random.choice(dis_HPO_list,int(nopheno.shape[0])) 
    all_data.loc[all_data['label']==0,'OMIM'] = random_pheno

    allD = list(all_data['OMIM'].drop_duplicates())
    dis_train, dis_test = train_test_split(allD,test_size=0.15) 
    
    # Collect the features for the genes and varaints: annotation.sh 
    all_data = all_data[['ID','OMIM','label']].copy()
    AnnotSV_benign = pd.read_csv(dbvar_benign_annotsv, sep='\t',low_memory=False)
    AnnotSV_path = pd.read_csv(dbvar_pathogenic_annotsv, sep='\t',low_memory=False)
    annotation = pd.concat([AnnotSV_benign,AnnotSV_path])

    all_data = all_data.merge(annotation, on=['ID'])
    
    # visulaize the class features 
    l = all_data[all_data['AnnotSV type']=='full']['label']
    print ("Number of classes:\n", l.value_counts())
    carrier_count = l.value_counts()
    sns.set(style="darkgrid")
    splot = sns.barplot(carrier_count.index, carrier_count.values, alpha=0.9)
    plt.title('Frequency Distribution of class')
    plt.ylabel('Number of Occurrences', fontsize=12)
    plt.xlabel('class', fontsize=12)
    for p in splot.patches:
        splot.annotate(format(p.get_height(), ''), (p.get_x() + p.get_width() / 2., p.get_height()), 
                    ha = 'center', va = 'center', xytext = (0, 6), textcoords = 'offset points')
    plt.show()
    #Number of classes: 0    36041, 1     5907
    
    # Find the enterz gene id 
    g = tidy_split(all_data, "Gene name", '/', keep=False)
    gene = pd.read_csv("data/homo_entrez_gene.txt", sep=' ',low_memory=False, header=None)
    gene.columns=['id','entrezgene','Gene name']
    all_data = all_data.merge(gene, on=['Gene name'], how='left')

    onto_types = ["go", "mp", "uberon", "hp","cl","all"]
    for i in onto_types:
        get_features(i)
    
    #save the mean and std to use for the testing 
    with open('data/mean_std_train.pkl', 'wb') as handle:
        pkl.dump(mean_std_train, handle)
        
