"""
Frequently used utility functions.
"""
import pandas as pd
import gensim
import pickle as pkl
import numpy


#------Preprocessing Functions--------------#
def tidy_split(df, column, sep, keep=False):
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
    df = df.dropna(subset=[column])
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


def get_geneID(gene, mg):
    """
    map gene name to enterez id
    
    Params
    ------
    gene : gene name as a panda column
    
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the gene Gene name and entrezgene
    """

    #mg = mygene.MyGeneInfo()
    out = mg.querymany(gene,
                       scopes='symbol',
                       fields='entrezgene',
                       species='human',
                       returnall=True,
                       as_dataframe=True)
    geneid = pd.DataFrame.from_dict(out['out']).reset_index()
    geneid = geneid.dropna(subset=['entrezgene'])
    geneid.rename(columns={"query": "Gene name"}, inplace=True)

    return geneid


def get_geneEnsID(gene, mg):
    """
    map gene name to ensamble id
    
    Params
    ------
    gene : gene name as a panda column
    
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the gene Gene name and entrezgene
    """

    #mg = mygene.MyGeneInfo()
    out = mg.querymany(gene,
                       scopes='entrezgene',
                       fields='ensembl.gene',
                       species='human',
                       returnall=True,
                       as_dataframe=True)
    geneid = pd.DataFrame.from_dict(out['out']).reset_index()
    geneid = geneid.dropna(subset=['ensembl.gene'])
    print(geneid.head())
    del geneid['_id'], geneid['_score'], geneid['ensembl']
    geneid.rename(columns={"query": "entrezgene"}, inplace=True)

    return geneid


class Stack(object):
    def __init__(self):
        self.items = []

    def is_empty(self):
        return self.items == []

    def peak(self):
        return self.items[len(self.items) - 1]

    def size(self):
        return len(self.items)

    def push(self, item):
        self.items.append(item)

    def pop(self):
        return self.items.pop()
