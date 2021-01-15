#!/usr/bin/python
from __future__ import print_function
import logging
import networkx as nx
import numpy as np
import gensim
import sys
import random
import json
import mygene
import os
import multiprocessing as mp
from threading import Lock
import pickle as pkl
from model import Rank_model
import torch
from torch._C import *
import torch.optim as optim
from scipy.stats import rankdata
from argparse import ArgumentParser
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset
from torch.autograd import Variable
from networkx.readwrite import json_graph
from gensim.models import Word2Vec, Phrases, phrases, KeyedVectors
import pandas as pd
from random import seed
seed(1)
from statistics import mean 
from collections import OrderedDict
import scipy.stats as ss
from operator import itemgetter
np.random.seed(42)
lock = Lock()

onto_type=sys.argv[1]
tmp='/tmp/'

path='./data/'
WALK_LEN=10 
N_WALKS=100 
global data_pairs
data_pairs = []

logger = logging.getLogger('my-logger')
logger.propagate = False

def run_random_walks(G, nodes,data_type, num_walks=N_WALKS):
    pairs = []
    seed(1)
    for count, node in enumerate(nodes):
        if G.degree(node) == 0:
            continue
        for i in range(num_walks):
            curr_node = node
            walk_accumulate=[]
            for j in range(WALK_LEN):
                next_node = random.choice(list(G.neighbors(curr_node)))
                type_nodes = G.edges[curr_node, next_node]["type"]
                if curr_node ==node:
                    walk_accumulate.append(curr_node)
                walk_accumulate.append(type_nodes)
                walk_accumulate.append(next_node)
                curr_node = next_node

            pairs.append(walk_accumulate)
    write_file(pairs,data_type)

def run_walk(nodes,G,data_type):
        global data_pairs
        #Start random walk
        number=30
        length = len(nodes) // number
    
        processes = [mp.Process(target=run_random_walks, args=(G, nodes[(index) * length:(index + 1) * length],data_type )) for index
                     in range(number-1)]
        processes.append(mp.Process(target=run_random_walks, args=(G, nodes[(number-1) * length:len(nodes) - 1], data_type)))
    
        for p in processes:
            p.start()
        for p in processes:
            p.join()
    
def write_file(pair,data_type):
        with lock:
            with open(tmp+data_type+"_walks.txt", "w") as fp:
                for p in pair:
                    for sub_p in p:
                        fp.write(str(sub_p)+" ")
                    fp.write("\n")

def rule(word, count, min_count):
    if word.startswith("p"):
        return gensim.utils.RULE_KEEP
    else:
        return gensim.utils.RULE_DEFAULT 
        
def hash(astring):
    return ord(astring[0])
                          
def gene_node_vector(graph, k,v,data_type):
    nodes_set=set()
    nodes_set.add(k)
    for da in v:  
        nodes_set.add(da)

    nodes= [n for n in nodes_set]
    nodes.sort(reverse=True)
    
    ego_graph = nx.ego_graph(graph, nodes[0], undirected=True, center=True, radius=100) 
    run_walk(nodes,ego_graph,data_type)#graph

    #start to train the word2vec models
    mymodel = KeyedVectors.load(path+onto_type+"_embedding")
    mymodel.wv.min_count=1
    mymodel.workers=20
    mymodel.hashfix=hash
    sentences =gensim.models.word2vec.LineSentence(tmp+data_type+"_walks.txt")
    mymodel.build_vocab(sentences, update=True, keep_raw_vocab=True,trim_rule=rule) #trim_rule=rule
    mymodel.train(sentences, total_examples=len(mymodel.wv.vocab) , epochs = 10)    
    return mymodel
    
def update_graph(typee,G,k,v):
    for i in v:
        G.add_edge(k.strip(), i.strip())
        G.edges[k.strip(), i.strip()]["type"] = "HasAssociation"
    return G
    
def convert_to_torch_format(data):
    data = np.array(data, dtype="float64")
    return torch.from_numpy(data).double()

def look_up_embed(embed_dic_, entity):
    entity_emb = embed_dic_[entity]
    return entity_emb

def predict(model_, disease, embed_dic_, ranking_score,genes):
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

if __name__ == '__main__':

    with open(path+onto_type+"_gene_set.pkl","rb") as f:
        genes=pkl.load(f)

    entities = set()
    for g in genes:
        entities.add(g) 
        
    with open(path+onto_type+"_dis_HPO.pkl","rb") as f:
        dis_HPO = pkl.load(f)

    bestmodel = path+onto_type+'_dl2vec_best_performance.pt'
    device = torch.device("cpu")
    #torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = Rank_model(num_feature=100).double()
    model.load_state_dict(torch.load(bestmodel))
    model.to(device)

    ranking_score={}
    for omim, hpo in dis_HPO.items():    
        print("loading data....") 
        G = json_graph.node_link_graph(json.load(open(path+onto_type+"_graph.json")))                              
        updateG = update_graph(onto_type,G,omim,hpo)
        word2vec_model = gene_node_vector(updateG, omim, hpo, onto_type)
        entities.add(omim)
        dic = dict()    
        for entity in entities:
            if entity in word2vec_model.wv.vocab:
                dic[entity] = word2vec_model.wv[entity]

        ranking_score = predict(model, omim , dic, ranking_score, genes)
        entities.remove(omim)
 
    dis={}
    for k, v in ranking_score.items():
        d = k.split('_')[0]
        g = k.split('_')[1]
        if d not in dis:
            dis[d]={}
            dis[d]['gene']=[]
            dis[d]['score']=[]
        else:
            dis[d]['gene'].append(g)
            dis[d]['score'].append(v)
    
    for k, v in dis.items():
        print(len(v['gene']))
        break

    #Quantiles ranking 
    final={}
    for k,v in dis.items():
        u = list(zip(v['gene'],v['score']))
        g = list(map(itemgetter(0), u))
        s = list(map(itemgetter(1), u))
        ranks = ss.rankdata(s)
        norm = [(float(i)-min(ranks))/(max(ranks)-min(ranks)) for i in ranks]
        g_s = dict(zip(g, norm))
        final[k] = g_s
        
    with open(path + typee + "_dl2vec.pkl","wb") as f:
        pkl.dump(final, f)
            
   