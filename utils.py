"""
Frequently used utility functions.
"""
import pandas as pd
import gensim
import pickle as pkl

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
    out = mg.querymany(gene, scopes='symbol', fields='entrezgene', species='human',returnall=True,as_dataframe=True)
    geneid = pd.DataFrame.from_dict(out['out']).reset_index()
    geneid = geneid.dropna(subset=['entrezgene'])
    #del geneid['_id'],geneid['_score'],geneid['notfound'],
    geneid.rename(columns={"query": "Gene name"} , inplace=True)

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
    out = mg.querymany(gene, scopes='entrezgene', fields='ensembl.gene', species='human',returnall=True,as_dataframe=True)
    geneid = pd.DataFrame.from_dict(out['out']).reset_index()
    geneid = geneid.dropna(subset=['ensembl.gene'])
    print(geneid.head())
    #del geneid['_id'],geneid['_score'], geneid['ensembl']
    geneid.rename(columns={"query": "entrezgene"} , inplace=True)

    return geneid


def get_HPO():
    """
    generate the human disease and hp association
    
    Params
    ------
    none
    
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the OMIM and HP annotation
    """
    
    dis_phe=dict()
    with open("../../EmbedPVP/phenotype_annotation.tab","r",encoding='utf-8') as f:
        for line in f.readlines():
            data=line.split("\t")
            if (data[0]=="OMIM")&(data[5][:4]=="OMIM"):
                try:
                    dis_phe[data[5].strip()].append(data[4].strip())
                except:
                    dis_phe[data[5].strip()]=[data[4].strip()]

    hp = pd.DataFrame.from_dict(dis_phe.items(),orient='columns')
    hp.rename(columns={0:'OMIM',1:'HPO'}, inplace=True)
    
    
    return hp

def get_GADO_test(data,ids): 
    """
    collect GADO score for the data set  
    
    Params
    ------
    data : pandas.DataFrame
        dataframe with OMIM, and Gene column 
    
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the a new column for GADO
    """
        
    gado = pd.read_csv("../data/annotsv_res/testing/"+ids+".txt")
    gado.drop_duplicates(inplace=True)
    gado = gado[['Ensg','Zscore']].copy()
    gado = gado.rename(columns={"Ensg": "ensembl.gene"})
    allData = data.merge(gado, on=['ensembl.gene'], how='left')

    return allData

    
def get_GADO(data): 
    """
    collect GADO score for the data set  
    
    Params
    ------
    data : pandas.DataFrame
        dataframe with OMIM, and Gene column 
    
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the a new column for GADO
    """
    with open('../annotation/GadoCommandline-1.0.1/gado_omim_genes.pickle', 'rb') as f:
        gado_omim_genes = pkl.load(f)
        
    omim = data[data['ensembl.gene'].notna()].copy()
    omim.drop_duplicates(inplace=True)
    omim = omim[['OMIM','ensembl.gene']].copy()
    omim = omim.groupby('OMIM')['ensembl.gene'].agg(lambda column: ",".join(column)).reset_index()
    omim = omim.set_index(['OMIM']).T.to_dict('list')

    res={}
    for k, v in omim.items():
        if k in gado_omim_genes:
            for g in v[0].split(','):
                if g in gado_omim_genes[k]:
                    kg=k+'_'+g
                    res[kg]=gado_omim_genes[k][g]
                    
    #convert to panda and merge 
    gado = pd.DataFrame(res.items(), columns=['OMIM_ensembl.gene', 'GADO_zscore'])
    gado['OMIM'] = gado['OMIM_ensembl.gene'].str.split('_').str.get(0)
    gado['ensembl.gene'] = gado['OMIM_ensembl.gene'].str.split('_').str.get(1)
    del gado['OMIM_ensembl.gene']
    allData = data.merge(gado, on=['OMIM','ensembl.gene'], how='left')

    return allData


def get_vectors(word2vec_model, data, d):
    """
    collect DL2vec vectors for the data set  
    
    Params
    ------
    word2vec_model : 
        word2vec_model with vectores for OMIM and gene 
    
    Returns
    -------
    pandas.DataFrame
        Returns a dataframe with the a new columns for the vectores
    """
    
    vec={}
    entities = data[d].tolist()
    model=gensim.models.Word2Vec.load(word2vec_model)
    for a in entities:
        try:
            ve =model[a]
            vec[a] = ve
        except:
            none=''
    vec_omim = pd.DataFrame.from_dict(vec, orient='index')
    vec_omim = vec_omim.reset_index()
    vec_omim1 = vec_omim.rename(columns={"index": d})
    #data = data[[d]].copy()
    allData = pd.merge(data, vec_omim1, on=[d], how='left')
    
    return allData

  #------Model Functions--------------#
  
      