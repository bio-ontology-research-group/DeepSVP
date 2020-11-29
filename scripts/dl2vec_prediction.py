from __future__ import print_function
import json
import numpy as np
import pickle as pkl
import gensim
from networkx.readwrite import json_graph
from argparse import ArgumentParser
from model import Rank_model
import torch
from torch._C import *
import random
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset
from torch.autograd import Variable
import torch.optim as optim
from scipy.stats import rankdata
import sys
import random
import argparse
import sys
import os

data_type = sys.argv[1]
bestmodel = sys.argv[2]


def convert_to_torch_format(data):
    data = np.array(data, dtype="float64")
    return torch.from_numpy(data).double()


def look_up_embed(embed_dic_, entity):
    entity_emb = embed_dic_[entity]
    return entity_emb


def predict(model_, disease, embed_dic_, ranking_score):
    model_.eval()
    testing_data = []
    positive_genes = disease_genes[disease]
    if disease in embed_dic_:
        disease_vec = look_up_embed(embed_dic_, disease)

        for gene in positive_genes:
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

    with open("training_pvp/dis_gene_training_cnv.pkl", "rb") as f:
        disease_genes = pkl.load(f)

    entities = set()
    for disease in disease_genes.keys():
        entities.add(disease)
        for gene in disease_genes[disease]:
            entities.add(gene)

    word2vec_model = gensim.models.Word2Vec.load("./DL2vec/output/" +
                                                 data_type + "_embedding")

    ranking_score = {}
    dic = dict()
    for entity in entities:
        if entity in word2vec_model.wv.vocab:
            dic[entity] = word2vec_model[entity]

    print(len(dic))

    diseases = []
    for disease in disease_genes.keys():
        diseases.append(disease)

    ranking_score = dict()
    device = torch.device("cpu")
    model = Rank_model(num_feature=100).double()
    model.load_state_dict(torch.load(bestmodel))
    model.to(device)

    for test_data in diseases:
        ranking_score = predict(model, test_data, dic, ranking_score)

    with open("./data/" + data_type + "_cnv_ranking_score_training.pkl",
              "wb") as f:
        pkl.dump(ranking_score, f)
