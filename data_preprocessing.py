#!/usr/bin/env python
# coding: utf-8

from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso
import matplotlib.pyplot as plt
import time
from matplotlib import offsetbox
import matplotlib.patheffects as PathEffects
get_ipython().magic(u'matplotlib inline')
import pathlib
import pandas as pd
from pandas.plotting import scatter_matrix
import pandas.plotting
from pandas.plotting import parallel_coordinates
import numpy as np
import seaborn as sns  # for intractve graphs
from sklearn.utils import class_weight
from sklearn.preprocessing import StandardScaler
seed = 42
np.random.seed(seed)
import pandas as pd
import gensim
import numpy as np
import re
import pickle as pkl
from numpy import array
from numpy import argmax
import utils
import mygene


# helper functions
def GCcontent(data, types):
    data[types] = data[types].astype(str)
    data[types] = data[types].str.split('(').str.get(0).astype(float)
    return data

def count(data, types, y):
    #data.insert(0,x,0)
    data[y] = data[types].apply(lambda x: len(str(x).split('/')))
    data[y] = data[y].astype(int)
    f_max = data.groupby('ID')[y].max().reset_index()
    del data[y]
    data = data.merge(f_max, on=['ID'])
    return data
    

def normalize_series(s):
    #return (s - s.mean())/s.std()
    return (s - s.min()) / (s.max() - s.min())


def preprocessing(features, columnNames, categorical):
    for col in columnNames:
        if col not in categorical:
            features[col + "_indicator"] = np.where(np.isnan(features[col]), 1, 0)
            features[col] = features[col].fillna(0)  #features[col].mean()
            #min max normlization
            features[col] = features[col].astype(float)
            features[col] = normalize_series(features[col])

    for cat in categorical:
        features[cat] = features[cat].fillna("undefined")
        features[cat] = features[cat].astype('category')
        features[cat] = features[cat].cat.codes
        features[cat] = features[cat].astype(float)
        features[cat] = normalize_series(features[cat])

    return features


def DL2vec(data, onto):
    #----------------------------------
    with open("data/" + onto + "_cnv_ranking_score_training.pkl", "rb") as f:
        prediction = pkl.load(f)
    prediction = pd.DataFrame(prediction.items(),
                              columns=['OMIM', 'DL2vec_score'])
    prediction['OMIM'] = prediction['OMIM'].astype(str)
    prediction['entrezgene'] = prediction['OMIM'].str.split('_').str.get(1)
    prediction['OMIM'] = prediction['OMIM'].str.split('_').str.get(0)
    prediction['DL2vec_score'] = prediction['DL2vec_score'].astype(float)
    prediction['entrezgene'] = prediction['entrezgene'].astype(str)
    #----------------------------------
    # model score
    prediction['entrezgene'] = prediction['entrezgene'].astype(str)
    data = data.merge(prediction, on=['entrezgene', 'OMIM'], how='left')
    return data


def collect_features(data, onto):
    # DL2vec score
    features = DL2vec(data, onto)

    # collect all the related scores
    columnNames = [
        'SV length', 'CDS length', 'tx length', 'HI_CGscore', 'TriS_CGscore',
        'synZ_ExAC', 'misZ_ExAC', 'pLI_ExAC', 'delZ_ExAC', 'dupZ_ExAC',
        'cnvZ_ExAC', 'HI_DDDpercent', 'GCcontent_left', 'GCcontent_right',
        'phastConsElements100way', 'genomicSuperDups', 'DL2vec_score',
        'StrVCTVRE', 'GADO_zscore'
    ]

    categorical = ['HI_CGscore', 'TriS_CGscore', 'SVtype']
    features[columnNames] = features[columnNames].astype(float).abs()

    features_mean = features.groupby('ID')[columnNames].mean().reset_index()
    features_max = features.groupby('ID')[columnNames].max().reset_index()

    features_match = features[[
                    'ID', 'class', 'SVtype', 'NumGenes', 'NumPromoters'
                    ]].copy()
    features_match.drop_duplicates(inplace=True)

    features_max = features_max.merge(features_match, on=['ID'])
    features_mean = features_mean.merge(features_match, on=['ID'])

    features_max = preprocessing(features_max, columnNames, categorical)
    features_mean = preprocessing(features_mean, columnNames, categorical)

    features_mean.drop_duplicates(inplace=True)
    features_max.drop_duplicates(inplace=True)

    return features_max, features_mean


def save_data(f_data1, f_data2, onto, op1, op2):
    f_data1.to_csv("../data/" + onto + "_" + op1 + "_" +
                   "training_predcnv.tsv",
                   sep='\t',
                   index=False)
    f_data2.to_csv("../data/" + onto + "_" + op2 + "_" +
                   "training_predcnv.tsv",
                   sep='\t',
                   index=False)



def feature_selection(onto, operation):
    f_data = pd.read_csv("../data/" + onto + "_" + operation + "_" +
                         "training_predcnv.tsv",
                         sep='\t')
    #Loading the dataset
    X = f_data.drop("class", 1).iloc[:, 1:f_data1.shape[1]]  #Feature Matrix
    y = f_data["class"]  #Target Variable
    reg = LassoCV(eps=0.001, tol=0.0001, normalize=True, max_iter=5000, cv=10)
    reg.fit(X, y)
    print("Best alpha using built-in LassoCV: ", reg.alpha_)
    print("Best score using built-in LassoCV: ", reg.score(X, y))
    coef = pd.Series(reg.coef_, index=X.columns)
    print("Lasso picked " + str(sum(coef != 0)) +
          " variables and eliminated the other " + str(sum(coef == 0)) +
          " variables")
    selected_features = (coef != 0)
    selected_features = list(
        selected_features[selected_features == True].index)
    global models_features
    k = onto + "_" + operation
    models_features[k] = selected_features
    





         
# # Dataset
# Annovar annotation
benign = pd.read_csv("data/benign.hg38_multianno.txt",
                     sep='\t',
                     low_memory=False)
                     
pathog = pd.read_csv("data/pathogenic.hg38_multianno.txt",
                    sep='\t',
                    low_memory=False)
    
annovar = pd.concat([benign, pathog])

annovar_feature = annovar[[
    'phastConsElements100way', 'genomicSuperDups', 'Otherinfo6']].copy()

annovar_feature['phastConsElements100way'] = annovar_feature['phastConsElements100way'].str.split('Score=').str.get(1).str.split(';').str.get(0)
    
annovar_feature['genomicSuperDups'] = annovar_feature['genomicSuperDups'].str.split('Score=').str.get(1).str.split(';').str.get(0)
    
annovar_feature.rename(columns={'Otherinfo6': "ID"}, inplace=True)

# AnnotSV annotation
benign = pd.read_csv("data/dbvar_benign.tsv",
                     sep='\t',
                     low_memory=False)  
                     
pathog = pd.read_csv("data/pathogenic_annotsv.tsv",
                     sep='\t',
                     low_memory=False)
                     
benign['class'] = 0
pathog['class'] = 1
data = [benign, pathog]
alldata = pd.concat(data)
alldata = alldata.merge(annovar_feature, on=['ID'], how='left')
print(alldata.shape)

# Find the enterz gene id
mg = mygene.MyGeneInfo()
gene = alldata[['Gene name']].copy()
gene = gene[~gene['Gene name'].str.contains('/', na=False)]
gene['Gene name'].drop_duplicates(inplace=True)
gene = utils.get_geneID(gene['Gene name'], mg)
alldata = alldata.merge(gene, on=['Gene name'], how='left')

# Ensamble gene
mg = mygene.MyGeneInfo()
gene = alldata[['entrezgene']].copy()
gene = gene[~gene['entrezgene'].str.contains('/', na=False)]
gene['entrezgene'].drop_duplicates(inplace=True)
gene = utils.get_geneEnsID(gene['entrezgene'], mg)
alldata = alldata.merge(gene, on=['entrezgene'], how='left')

# Process the pathoginc and bengin variants
allData = alldata.copy()
allData['SVtype'] = allData['INFO'].str.split('SVTYPE=').str.get(1).str.split(';').str.get(0)

allData = allData.loc[~allData['SVtype'].str.contains('IN')]
allData['OMIM'] = allData['INFO'].str.findall(r'OMIM:\d{6}')
allData = utils.tidy_split(allData, "OMIM", ',', keep=False)
allData['OMIM'] = allData['OMIM'].astype(str)
allData.loc[allData["OMIM"] == "[]", "OMIM"] = np.nan
allData['OMIM'] = allData['OMIM'].str.replace(r"\['", "").str.replace(r"\]", "").str.replace(r"\'", "").str.strip()
allData.drop_duplicates(inplace=True)

# pathoginc but not causing the diseases, if class == 1 and no pheno, replace it with 0
allData.loc[allData['OMIM'].isnull(), 'class2'] = 0
nopheno = allData[(allData['OMIM'].isnull())]
random_pheno = np.random.choice(allData.OMIM.dropna(), int(y.OMIM.isna().sum()))

allData.loc[allData.OMIM.isnull(), 'OMIM'] = random_pheno
allData.loc[(allData['class'] == 1.0) & (allData['class2'] == 0.0),'class'] = 0

# visulaize the class features
get_ipython().magic(u'matplotlib inline')
import seaborn as sns
import matplotlib.pyplot as plt
carrier_count = allData['class'].value_counts()
sns.set(style="darkgrid")
splot = sns.barplot(carrier_count.index, carrier_count.values, alpha=0.9)
plt.title('Frequency Distribution of SVtype')
plt.ylabel('Number of Occurrences', fontsize=12)
plt.xlabel('class', fontsize=12)
for p in splot.patches:
    splot.annotate(format(p.get_height(), ''),
                   (p.get_x() + p.get_width() / 2., p.get_height()),
                   ha='center',
                   va='center',
                   xytext=(0, 6),
                   textcoords='offset points')
print(allData['class'].value_counts())
plt.show()




# collect GADO score
features = utils.get_GADO(allData)

# collect StrVCTVRE score
benign_StrVCTVRE = pd.read_csv("data/benign_StrVCTVRE.vcf",
                               comment='#',
                               header=None,
                               sep='\t')
path_StrVCTVRE = pd.read_csv("data/path_StrVCTVRE.vcf",
                             comment='#',
                             header=None,
                             sep='\t')
StrVCTVRE = pd.concat([benign_StrVCTVRE, path_StrVCTVRE])
StrVCTVRE['StrVCTVRE'] = StrVCTVRE[7].str.split(';StrVCTVRE=').str.get(1)
StrVCTVRE = StrVCTVRE.rename(columns={2: "ID"})
StrVCTVRE['StrVCTVRE'] = StrVCTVRE['StrVCTVRE'].replace({
    'not_exonic':
    np.nan,
    'less_than_50bp':
    np.nan,
    'not_dup_or_del':
    np.nan,
    'missing_END_or_SVTYPE':
    np.nan,
    'not_valid_chrom':
    np.nan
})
StrVCTVRE = StrVCTVRE[['ID', 'StrVCTVRE']].copy()
features = features.merge(StrVCTVRE, on=['ID'])
features['TriS_CGscore'] = features['TriS_CGscore'].replace(
    'Not yet evaluated', np.NaN)

# other features: Number of Genes and promoters, and the GC content
features = count(features, 'Gene name', 'NumGenes')
features = count(features, 'promoters', 'NumPromoters')
features = GCcontent(features, 'GCcontent_right')
features = GCcontent(features, 'GCcontent_left')

onto_types = [
    "go", "mp", "uberon", "hp", "go_ppi", "mp_ppi", "uberon_ppi", "hp_ppi"
]
for onto in onto_types:
    print("Generate a training data for: ", onto)
    max_score, mean_score = collect_features(features, onto)
    save_data(max_score, mean_score, onto, 'max', 'mean')
    print('Done.')

# # Feature Selection Embedded Method
models_features = {}
for onto in onto_types:
    print("Feature Selection Embedded Method: ", onto)
    feature_selection(onto, 'max')
    feature_selection(onto, 'mean')
    print("------------------------------------------")

with open('data/features.pickle', 'wb') as handle:
    pkl.dump(models_features, handle)


imp_coef = coef.sort_values()
import matplotlib
matplotlib.rcParams['figure.figsize'] = (10.0, 10.0)
locs, labels = plt.xticks()
plt.xticks(np.arange(-9, 9, step=0.01))
imp_coef.plot(kind="barh", linewidth=0.1)
plt.title("Feature importance using Lasso Model")
