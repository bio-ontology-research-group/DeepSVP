import os
import mygene
import numpy
import pandas as pd
import numpy as np 
import tensorflow as tf
import pickle
import gensim
import random
import pickle as pkl
import sys
from pandas.api.types import CategoricalDtype 
from statistics import mean 
import sklearn.metrics as metrics
import click as ck
pd.options.mode.chained_assignment = None  


@ck.command()
@ck.option(
    '--onto', '-m', default=f'mp',
    help='Ontology model, one of the following (go , mp , hp, cl, uberon, union')
@ck.option(
    '--ope', '-ag', default=f'max',
    help='Aggregation method for the genes within CNV (max , mean)')
@ck.option(
    '--subsets', '-sub', default='full',
    help='Phenotypes subsets (full)')
         
def main(onto,ope,subsets):
    allrank = []
    samples=[]
    tophits_auc= pd.DataFrame()
    count=0
    aucs=[]
    id_score={}              
    annotsv_path='data/benchmark/Annootsv_postive.tsv'
    annotsv_ben='data/benchmark/Annotsv_negative.tsv'        
    annotsv_p = pd.read_csv(annotsv_path, sep='\t', low_memory=False)  
    annotsv_b = pd.read_csv(annotsv_ben, sep='\t', low_memory=False)
    annotsv_b.drop(annotsv_b.iloc[:, 12:2517], inplace = True, axis = 1) 
    annotsv_b['ID']=annotsv_b['ID'].astype(str)    

    annotsv_p = tidy_split(annotsv_p, "Gene name", '/', keep=False)
    gene = pd.read_csv("data/benchmark/homo_entrez_gene.txt", sep=' ',low_memory=False, header=None)
    gene.columns=['id','entrezgene','Gene name']
    final_path = annotsv_p.merge(gene, on=['Gene name'],how='left')
    final_path['class']=1

    annotsv_b = tidy_split(annotsv_b, "Gene name", '/', keep=False)
    gene = pd.read_csv("data/benchmark/homo_entrez_gene.txt", sep=' ',low_memory=False, header=None)
    gene.columns=['id','entrezgene','Gene name']
    final_ben =  annotsv_b.merge(gene, on=['Gene name'],how='left')
    final_ben['class']=0

    with open('data/benchmark/patients', "r") as f: 
        for line in f.readlines():
            samples.append(line.split()[0])
            
    comb_data= pd.concat([final_ben,final_path])
    all_genes=list(comb_data['entrezgene'].drop_duplicates())
    comb_data.drop_duplicates(inplace=True)

    if onto != 'union':
        dlres=get_allscore(onto,samples,all_genes,subsets)
        res_go={}
        res_mp={}
        res_cl={}
        res_ub={}
        res_hp={}
    else:
        dlres={}
        onto_types = ["go", "mp", "hp", "cl","uberon"]
        res_go=get_allscore(onto_types[0],samples,all_genes,subsets)
        res_mp=get_allscore(onto_types[1],samples,all_genes,subsets)
        res_hp=get_allscore(onto_types[2],samples,all_genes,subsets)
        res_cl=get_allscore(onto_types[3],samples,all_genes,subsets)
        res_ub=get_allscore(onto_types[4],samples,all_genes,subsets)          
    j=0
    
    for ID in samples:
        allData = final_path[final_path['ID']==ID]
        allData = final_ben.append(allData)
        allData['num'] = allData['Gene name'].apply(lambda x: len(str(x).split('/'))).astype(int)
        allData = pd.concat([func(row) for _, row in allData.iterrows()], ignore_index=True, axis=1).T
        allData['ID'] = allData['AnnotSV ID']    
        allData['entrezgene']= allData['entrezgene'].astype(str)
        allData['entrezgene']= allData.entrezgene.apply(lambda x: x[0:-2])
        features = collect_features(allData, onto, ope, ID, dlres,  res_go,res_mp,res_hp,res_cl,res_ub)
        features.drop_duplicates(inplace=True) 
        res = test_model_bench(onto, features, ope) 
        f='data/DeepSVP_'+onto+'_'+ope+'/'+str(j)+'_'+subsets+'.tsv'
        res.to_csv(f, index=False, sep='\t')       
        j+=1    

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

def GCcontent(data, types):
    data[types] = data[types].astype(str)
    data[types] = data[types].str.split('(').str.get(0).astype(float)
    return data

def preprocessing_numeric(features,cols,onto, ope):
    with open('data/benchmark/'+onto+ope+'_full_mean_std_scaler.pkl', 'rb') as handle:   
        mean_std_train = pickle.load(handle)
    for column in cols:
        if not column.startswith('HI_CGscore') and not column.startswith('TriS_CGscore'):
            features[column+"_indicator"] = np.where(features[column].isnull(), 1, 0)
            features[column] = (features[column]-mean_std_train[column]['mean'])/mean_std_train[column]['std']
            features[column] = features[column].fillna(0)                
    return features

def preprocessing_categorical(features,cat):  
        features[cat] = features[cat].fillna('undefined') 
        values= [3.0,2.0,1.0,0.0,40.0,30.0,'undefined']
        features[cat] = features[cat].astype(pd.CategoricalDtype(values))
        features = pd.concat([features,pd.get_dummies(features[cat], prefix=cat)],axis=1)
        features = features.drop([cat],axis=1)
        return features

def get_scores(dl2vec_scores,dis, gene):
    g=str(gene).strip('.0')
    if dis in dl2vec_scores:
        if g in dl2vec_scores[dis]:
            results=dl2vec_scores[dis][g]
        else:
            results=np.nan
    else:
        results=np.nan
        
    return results
    
    
def test_model_bench(ontology, test_data, op):       
    cleaned_df = test_data.iloc[:, 2:test_data.shape[1]].values
    tf.keras.backend.clear_session()
    model = tf.keras.models.load_model("data/models_final/final/"+ontology+"_"+op+".h5") 
    test_pred = model.predict(cleaned_df)
    test_result = pd.DataFrame()
    test_result['ID'] = test_data['ID']
    test_result['score'] = test_pred
    test_result['label'] =  test_data['class']
    test_result['score'] = test_result['score'].astype(float) 
    return test_result

    
def collect_features(data, onto, operation, ids, dl2vec, res_go,res_mp,res_hp,res_cl,res_ub):
    # 1. Gene Features 
    #---------------------
    # get DL2vec scores
    Gene_features = data[data['AnnotSV type']=='split']
    if (onto == 'union'):
        Gene_features['go_score'] = Gene_features.apply(lambda x: get_scores(res_go,x['ID'], x['entrezgene']), axis=1)
        Gene_features['mp_score'] = Gene_features.apply(lambda x: get_scores(res_mp,x['ID'], x['entrezgene']), axis=1)
        Gene_features['cl_score'] = Gene_features.apply(lambda x: get_scores(res_cl,x['ID'], x['entrezgene']), axis=1)
        Gene_features['hp_score'] = Gene_features.apply(lambda x: get_scores(res_hp,x['ID'], x['entrezgene']), axis=1)
        Gene_features['uberon_score'] = Gene_features.apply(lambda x: get_scores(res_ub,ids, x['entrezgene']), axis=1)                                                                
        GF = ['ID','HI_DDDpercent','go_score','mp_score','cl_score','hp_score','uberon_score','CDS length','tx length','synZ_ExAC','misZ_ExAC', 'pLI_ExAC', 'delZ_ExAC', 'dupZ_ExAC', 'cnvZ_ExAC']
    else:
        Gene_features[onto+'_score'] = Gene_features.apply(lambda x: get_scores(dl2vec,ids, x['entrezgene']), axis=1)
        GF = ['ID','HI_DDDpercent',onto+'_score','CDS length','tx length','synZ_ExAC','misZ_ExAC', 'pLI_ExAC', 'delZ_ExAC', 'dupZ_ExAC', 'cnvZ_ExAC']
            
    #1. CNV_features
    #---------------------
    CF=['ID','class','SV type','SV length','GCcontent_right', 'GCcontent_left','NumGenes','NumPromoters','HI_CGscore','TriS_CGscore']
    cnv_features = data[data['AnnotSV type']=='full']
    cnv_features['NumGenes'] = cnv_features['Gene name'].apply(lambda x: len(str(x).split('/'))).astype(int)
    cnv_features['NumPromoters'] = cnv_features['promoters'].apply(lambda x: len(str(x).split('/'))).astype(int)
    cnv_features['TriS_CGscore'] = cnv_features['TriS_CGscore'].replace('Not yet evaluated',np.NaN)
    
    cnv_features = GCcontent(cnv_features, 'GCcontent_right')
    cnv_features = GCcontent(cnv_features, 'GCcontent_left')
    cnv_features['SV length'] = cnv_features['SV length'].abs()
    cnv_features['SV type'] = pd.Categorical(cnv_features['SV type'])
    cnv_features['SV type'] =  cnv_features['SV type'].cat.codes
    cnv_features = cnv_features[CF].copy()
    cnv_features[CF[3:]]=cnv_features[CF[3:]].astype(float)    

    #2. Gene Features 
    Gene_features = Gene_features[GF].copy()
    Gene_features[GF[1:]]=Gene_features[GF[1:]].astype(float)

    if (operation == "max"):       
       gene_features = Gene_features.groupby('ID')[GF[1:]].max().reset_index()
                 
    elif (operation == "avg"):
        gene_features = Gene_features.groupby('ID')[GF[1:]].mean().reset_index()
    else:
        print('error! choose avg, max')
        exit()
    # merge the features     
    features = cnv_features.merge(gene_features, on=['ID'], how='left') 
    col = list(features.loc[:,'SV length':].columns)
    categorical=['HI_CGscore','TriS_CGscore']
    for i in categorical:
        features = preprocessing_categorical(features,i)
            
    features = preprocessing_numeric(features,col,onto, operation)    
    return features

def get_allscore(onto,samples,all_genes,subsets):
    with open("data/benchmark/"+onto+"_"+subsets+"_dl2vec.pkl","rb") as f:
        dl2vec=pkl.load(f)  
    res={}
    l=[]
    for i in samples:
        if i not in dl2vec:
            l.append(i)
            continue
        for g in all_genes:
            gg=str(g).strip('.0')
            if gg in dl2vec[i]:
                if i not in res:
                    res[i]={}
                res[i][gg] =  dl2vec[i][gg]
    return res

def func(row):
   if row['num'] == 1:
        row2 = row.copy()
        row2['AnnotSV type']='split'
        return pd.concat([row, row2], axis=1)       
   return row
                         
if __name__ == '__main__':
    main() 
    
