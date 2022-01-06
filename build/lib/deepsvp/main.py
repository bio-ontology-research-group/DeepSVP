#!/usr/bin/env python
import click as ck
import numpy as np
import pandas as pd
import gzip
import os
import sys
import logging
import numpy
import tensorflow as tf
import gensim
import random
import pickle as pkl
from pandas.api.types import CategoricalDtype
from statistics import mean
import sklearn.metrics as metrics
from progress.bar import Bar
import time
import subprocess
import argparse
import os.path
import networkx as nx
import pickle as pkl
import json
import tempfile
import shutil
import pdb
from os import path
import networkx as nx
import json
import multiprocessing as mp
from threading import Lock
import deepsvp.model
from deepsvp.model import Rank_model
import torch
import shutil
from torch._C import *
import torch.optim as optim
from scipy.stats import rankdata
from argparse import ArgumentParser
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset
from torch.autograd import Variable
from networkx.readwrite import json_graph
from gensim.models import Word2Vec, Phrases, phrases, KeyedVectors
from random import seed
logging.getLogger("urllib3").setLevel(logging.WARNING)
np.random.seed(42)
from statistics import mean
from operator import itemgetter
import scipy.stats as ss
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None  # default='warn'
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"
tmp = tempfile.mkdtemp(prefix='DeepSVP', suffix='/')
bar = Bar(max=4, fill='=', suffix='%(percent)d%%')
lock = Lock()
WALK_LEN = 25
N_WALKS = 100

logger = logging.getLogger('my-logger')
logger.propagate = False
logger.disabled = True

@ck.command()
@ck.option('--data-root', '-d',
           default='data/',
           help='Data root folder',
           required=True)
@ck.option('--in-file', '-i', help='Annotated Input file', required=True)
@ck.option('--hpo','-p',
           help='List of phenotype ids separated by commas',
           required=True)
@ck.option('--maf_filter','-maf',
           help='Allele frequency filter using gnomAD and 1000G default<=0.01',
           default=0.01)
@ck.option('--model_type','-m',
            default='mp',
            help='Ontology model, one of the following (go , mp , hp, cl, uberon, union), default=mp')
@ck.option('--aggregation','-ag',
            help='Aggregation method for the genes within CNV (max or mean) default=max',
            default='max')
@ck.option('--outfile','-o',
            default='cnv_results.tsv',
            help='Output result file')

def main(data_root, in_file, hpo, maf_filter, model_type, aggregation, outfile):
    # Check data folder and required files
    """DeepSVP: A phenotype-based tool to prioritize caustive CNV using WGS data and Phenotype/Gene Functional Similarity"""
    try:
        if os.path.exists(data_root):
            in_file = os.path.join(data_root, in_file)
            model_file = model_type + '_' + aggregation + '.h5'
            model_file = os.path.join(data_root, model_file)
            if not os.path.exists(in_file):
                raise Exception(
                    f'Annotated Input file ({in_file}) is missing!')
            if not os.path.exists(model_file):
                raise Exception(f'Model file ({model_file}) is missing!')
        else:
            raise Exception(f'Data folder {data_root} does not exist!')
    except Exception as e:
        logging.error(e)
        sys.exit(1)

    state = 'START'
    while state != 'DONE':
        # Read input data, Load and Run pheno model
        bar.next()
        load_pheno_model(in_file, hpo, model_type, data_root)

        bar.next()
        print(" Phenotype prediction... ")

        # Load and Run cnv model
        bar.next()
        output = load_cnv_model(in_file, model_type, aggregation, data_root, maf_filter)
        print(" CNV Prediction... ")

        # Write the results to a file
        bar.next()
        print_results(output, outfile)
        print(' DONE! You can find the prediction results in the output file:',
              outfile)

        state = 'DONE'


def load_pheno_model(in_file, hpo, model_type, data_root):
    """
    This function load input data and formats it and run phenotype model
    """
    with open(tmp + "pheno.txt", "w") as fp:
        for pheno in hpo.split(','):
            fp.write('patient' + " " + "<http://purl.obolibrary.org/obo/" +pheno.replace(':', '_') + ">")
            fp.write("\n")
        fp.close()
    print(" Reading the input phenotypes...")
    pheno = tmp + "pheno.txt"
    #load_pheno_model
    if (model_type == 'union'):
        for i in ['go', 'mp', 'hp', 'cl', 'uberon']:
            pheno_model(in_file, i, pheno, tmp, data_root)
    else:
        pheno_model(in_file, model_type, pheno, tmp, data_root)


def print_results(results, out_file):
    """
    Write results to a file
    """
    results.to_csv(out_file, sep='\t', index=False)


#Pheno model Functions:-------------------------------------------------------
def run_random_walks(G, nodes, data_type, num_walks=N_WALKS):
    """
    This function to run the random walks
    """
    pairs = []
    seed(1)
    for count, node in enumerate(nodes):
        if G.degree(node) == 0:
            continue
        for i in range(num_walks):
            curr_node = node
            walk_accumulate = []
            for j in range(WALK_LEN):
                next_node = random.choice(list(G.neighbors(curr_node)))
                type_nodes = G.edges[curr_node, next_node]["type"]
                if curr_node == node:
                    walk_accumulate.append(curr_node)
                walk_accumulate.append(type_nodes)
                walk_accumulate.append(next_node)
                curr_node = next_node

            pairs.append(walk_accumulate)
    write_file(pairs, data_type)


def run_walk(nodes, G, data_type):
    """
    This function to run the walks in parallel
    """
    number = 30
    length = len(nodes) // number
    processes = [
        mp.Process(target=run_random_walks,
                   args=(G, nodes[(index) * length:(index + 1) * length],
                         data_type)) for index in range(number - 1)
    ]
    processes.append(
        mp.Process(target=run_random_walks,
                   args=(G, nodes[(number - 1) * length:len(nodes) - 1],
                         data_type)))
    for p in processes:
        p.start()
    for p in processes:
        p.join()


def write_file(pair, data_type):
    """
    This function to write the walk to new file
    """
    with lock:
        with open(tmp + data_type + "_walks.txt", "a") as fp:
            for p in pair:
                for sub_p in p:
                    fp.write(str(sub_p) + " ")
                fp.write("\n")

def hash(astring):
    return ord(astring[0])

def gene_node_vector(data_folder, graph, entity_list, data_type):
    nodes_set = set()
    with open(entity_list, "r") as f:
        for line in f.readlines():
            data = line.strip().split()
            for da in data:
                nodes_set.add(da)

    nodes = [n for n in nodes_set]
    nodes.sort(reverse=True)
    print("Update the graph statring from patient node.")
    ego_graph = nx.ego_graph(graph,
                             nodes[0],
                             undirected=True,
                             center=True,
                             radius=100)
    run_walk(nodes, ego_graph, data_type)
    print("start to train the word2vec models")
    mymodel = KeyedVectors.load(data_folder + data_type + "_embedding")
    mymodel.wv.min_count = 0
    mymodel.workers = 10
    mymodel.hashfix = hash
    sentences = gensim.models.word2vec.LineSentence(tmp + data_type +"_walks.txt")
    mymodel.build_vocab(sentences, update=True, keep_raw_vocab=True)
    mymodel.train(sentences, total_examples=len(mymodel.wv.vocab), epochs=10)
    return mymodel


def update_graph(filename, typee, G):
    """
    This function to add the patient phenotypes data
    """
    with open(filename, "r") as f:
        for line in f.readlines():
            entities = line.split(' ')
            G.add_edge(entities[0].strip(), entities[1].strip())
            G.edges[entities[0].strip(), entities[1].strip()]["type"] = "HasAssociation"
    return G

def convert_to_torch_format(data):
    data = np.array(data, dtype="float64")
    return torch.from_numpy(data).double()

def look_up_embed(embed_dic_, entity):
    entity_emb = embed_dic_[entity]
    return entity_emb

def predict(model_, disease, embed_dic_, ranking_score, genes):
    """
    This function to predict phenotype score
    """
    model_.eval()
    testing_data = []
    if disease in embed_dic_:
        disease_vec = look_up_embed(embed_dic_, disease)
        for gene in genes:
            if gene in embed_dic_:
                gene_vec = look_up_embed(embed_dic_, gene)
                testing_data = [gene_vec, disease_vec]
                ids = disease + "_" + gene
                testing_data = convert_to_torch_format(testing_data)
                testing_data = Variable(testing_data)
                test_result = model_.predict(testing_data).data
                ranking_score[ids] = test_result.item()
    return ranking_score

# decorater used to block function printing to the console
def blockPrinting(func):
    def func_wrapper(*args, **kwargs):
        # block all printing to the console
        sys.stdout = open(os.devnull, 'w')
        # call the method in question
        value = func(*args, **kwargs)
        # enable all printing to the console
        sys.stdout = sys.__stdout__
        # pass the return value of the method back
        return value

    return func_wrapper

@blockPrinting
def pheno_model(in_file, typee, filename, tmp, data_folder):
    print("loading data....")
    G = json_graph.node_link_graph(json.load(open(data_folder + typee + "_graph.json")))
    #Rank all the genes
    with open(data_folder + typee + "_gene_set.pkl","rb") as f:
        genes=pkl.load(f)
        
    entities = set()
    entities.add('patient')
    for g in genes:
        entities.add(g)

    updateG = update_graph(filename, typee, G)
    word2vec_model = gene_node_vector(data_folder, updateG, filename, typee)
    ranking_score = {}
    dic = dict()
    for entity in entities:
        if entity in word2vec_model.wv.vocab:
            dic[entity] = word2vec_model.wv[entity]

    ranking_score = dict()
    bestmodel = data_folder + typee + '_dl2vec_best_performance.pt'
    device = torch.device("cpu")
    model = Rank_model(num_feature=100).double()
    model.load_state_dict(torch.load(bestmodel))
    model.to(device)

    ranking_score = predict(model, 'patient', dic, ranking_score, genes)
    
    dis={}
    for k, v in ranking_score.items():
        d = k.split('_')[0]
        g = k.split('_')[1]
        if d not in dis:
            dis[d]={}
            dis[d]['gene']=[]
            dis[d]['score']=[]
        dis[d]['gene'].append(g)
        dis[d]['score'].append(v)
    
    final={}
    for k,v in dis.items():
        u = list(zip(v['gene'],v['score']))
        g = list(map(itemgetter(0), u))
        s = list(map(itemgetter(1), u))
        ranks = ss.rankdata(s)
        norm = [(float(i)-min(ranks))/(max(ranks)-min(ranks)) for i in ranks]
        g_s = dict(zip(g, norm))
        final[k] = g_s

    with open(tmp + typee + "_dl2vec.pkl", "wb") as f:
        pkl.dump(final, f)

#CNV Functions:-------------------------------------------------------
def cnv_model(model_type, features, operation, data_root):
    df = features.iloc[:, 1:features.shape[1]].values
    tf.keras.backend.clear_session()
    model = tf.keras.models.load_model(data_root + model_type + "_" + operation + ".h5")
    pred = model.predict(df)
    result = pd.DataFrame()
    result['ID'] = features['ID']
    result['DeepSVP_Score'] = pred
    result['DeepSVP_Score'] = result['DeepSVP_Score'].astype(float)
    return result

def GCcontent(data, types):
    data.loc[data.index, types] = data[types].astype(str)
    data.loc[data.index, types] = data[types].str.split('(').str.get(0).astype(float)
    return data

def preprocessing_numeric(features,cols,mean_std_train):
    for column in cols:
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

def get_scores(dl2vec_scores, gene):
    k='patient'
    g=str(gene).strip('.0')
    if g in dl2vec_scores[k]:
        results=dl2vec_scores[k][g]
    else:
        results=np.nan
    return results

def collect_features(data, onto, operation, data_root):
    # Gene Features
    #---------------------
    Gene_features = data[data['AnnotSV type'] == 'split']
    if (onto == 'union'):
        onto_types = ["go", "mp", "uberon", "hp", "cl"]
        for i in onto_types:
            with open(tmp + i +"_dl2vec.pkl","rb") as f: 
                dl2vec=pkl.load(f)
            Gene_features.loc[Gene_features.index,i+'_score'] = Gene_features.apply(lambda x: get_scores(dl2vec, x['entrezgene']), axis=1)            
        GF = ['ID','HI_DDDpercent','go_score','mp_score','cl_score','hp_score','uberon_score',
              'CDS length','tx length','synZ_ExAC','misZ_ExAC', 'pLI_ExAC', 'delZ_ExAC', 'dupZ_ExAC', 'cnvZ_ExAC']
    else:
        with open(tmp + onto+ "_dl2vec.pkl","rb") as f: 
            dl2vec=pkl.load(f)
        Gene_features.loc[Gene_features.index,onto+'_score'] = Gene_features.apply(lambda x: get_scores(dl2vec, x['entrezgene']), axis=1) 
        GF = ['ID','HI_DDDpercent',onto+'_score','CDS length','tx length','synZ_ExAC','misZ_ExAC', 
              'pLI_ExAC', 'delZ_ExAC', 'dupZ_ExAC', 'cnvZ_ExAC']

    #1. CNV_features
    #---------------------
    CF=['ID','SV type','SV length','GCcontent_right', 'GCcontent_left','NumGenes','NumPromoters','HI_CGscore','TriS_CGscore']
    cnv_features = data[data['AnnotSV type']=='full']
    cnv_features.loc[cnv_features.index,'NumGenes'] = cnv_features['Gene name'].apply(lambda x: len(str(x).split('/'))).astype(int)
    cnv_features.loc[cnv_features.index,'NumPromoters'] = cnv_features['promoters'].apply(lambda x: len(str(x).split('/'))).astype(int)
    cnv_features.loc[cnv_features.index,'TriS_CGscore'] = cnv_features['TriS_CGscore'].replace('Not yet evaluated',np.NaN)
    cnv_features = GCcontent(cnv_features, 'GCcontent_right')
    cnv_features = GCcontent(cnv_features, 'GCcontent_left')
    cnv_features.loc[cnv_features.index,'SV length'] = cnv_features['SV length'].abs()
    cnv_features = cnv_features[(cnv_features['SV type'] == 'DUP') | (cnv_features['SV type'] == 'DEL')]
    cnv_features.loc[cnv_features.index,'SV type'] = pd.Categorical(cnv_features['SV type'])
    cnv_features['SV type'] =  cnv_features['SV type'].cat.codes
    cnv_features = cnv_features[CF].copy()
    cnv_features[CF[3:]]=cnv_features[CF[3:]].astype(float) 
    
    #2. Gene Features 
    Gene_features = Gene_features[GF].copy()
    Gene_features[GF[1:]]=Gene_features[GF[1:]].astype(float)
    if (operation == "max"):
       gene_features = Gene_features.groupby('ID')[GF[1:]].max().reset_index()              
    elif (operation == "mean"):
        gene_features = Gene_features.groupby('ID')[GF[1:]].mean().reset_index()
    else:
        print('ERROR! choose mean, max operation')
        exit()

    # merge the features     
    features = cnv_features.merge(gene_features, on=['ID'], how='left') 
    categorical=['HI_CGscore','TriS_CGscore']
    for i in categorical:
        features = preprocessing_categorical(features,i)    
    col = list(features.loc[:,'SV length':].columns)
    
    with open(data_root+'mean_std_train.pkl', 'rb') as handle:
        mean_std_train = pkl.load(handle)
    
    features = preprocessing_numeric(features,col,mean_std_train)
    
    return features

@blockPrinting
def load_cnv_model(in_file, model_type, aggregation, data_root, maf):
    """
    This function to run the final cnv model
    """
    data = pd.read_csv(in_file, sep='\t', low_memory=False)
    data['1000g_AF'].fillna(0, inplace=True)
    data['GD_AF'].fillna(0, inplace=True)
    data = data[(data['1000g_AF'] <= maf) & (data['GD_AF'] <= maf)]
    genes = pd.read_csv(data_root + "homo_entrez_gene.txt", sep=' ',low_memory=False, header=None)
    genes.columns=['id','entrezgene','Gene name']
    data = data.merge(genes, on=['Gene name'], how='left')
    features = collect_features(data, model_type, aggregation, data_root)
    features.drop_duplicates(inplace=True)
    results = cnv_model(model_type, features, aggregation, data_root)
    data = data[data['AnnotSV type'] == 'full']
    data = data[['ID', 'SV chrom', 'SV start', 'SV end', 'SV type', 'REF', 'ALT','Gene name']].copy()
    results = data.merge(results, on=['ID']).drop_duplicates()
    results.fillna("-", inplace=True)
    results = results.sort_values(by='DeepSVP_Score', ascending=0)
    results['rank'] = range(1, 1 + len(results))

    shutil.rmtree(tmp)
    
    return results

if __name__ == '__main__':
    main() 
