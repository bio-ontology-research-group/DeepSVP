import pandas as pd
import numpy as np 
import sys 
import glob
from sklearn.utils import shuffle
import click as ck
pd.options.mode.chained_assignment = None  
import sklearn.metrics
import statistics as stat

@ck.command()
@ck.option(
    '--method', '-m', default=f'AnnotSV',
    help='Method to evaluate')
@ck.option(
    '--data-dir', '-dir', default=f'data/',
    help='Predictions directories')
@ck.option(
    '--num-p', '-n', default=1503,
    help='Number of patients')
         
def main (method, data_dir , num_p):
    prediction=glob.glob(data_dir)
    auc_ranks={}
    tp={}
    fp={}
    allrank=[]
    a=[]
    allscores = pd.DataFrame()
    total_num_var = 2392
    np.random.seed(123) 
    for pred in prediction:   
        df = pd.read_csv(pred, sep='\t')
        df = df.sort_values(by='score', ascending=False)
        df['rank'] = range(1, 1+len(df))
        df = df.groupby('score', as_index=False).agg(lambda x: list(x))
        df['label'] = df['label'].explode().astype(int).max(level=0)
        random_index = [np.random.randint(len(ranks)) for ranks in df['rank']]
        ranks = [ranks[idx] for idx, ranks in zip(random_index, df['rank'])]
        df = df.assign(rank=ranks)
        if method == 'AnnotSV':
            pr_auc = sklearn.metrics.average_precision_score(df['label'], df['rank'], average='micro')
        else:
            pr_auc = sklearn.metrics.average_precision_score(df['label'], df['score'], average='micro')
        a.append(pr_auc)
        index = df[df['label']==1]        
        index['ID'] = pred.split('/')[11].split('.')[0].split('_')[0]
        index['ID'] = index['ID'].astype(int)
        allscores = allscores.append(index)             
        rank = index['rank'].values[0]
        if rank not in auc_ranks:
            auc_ranks[rank] = 0
        auc_ranks[rank] += 1        
        allrank.append(rank)

    print("Benchmark evaluation for: ", method)
    top_k(allrank, 1, num_p)
    top_k(allrank, 10, num_p)
    top_k(allrank, 30, num_p)
    aucs = compute_rank_roc(auc_ranks,total_num_var)
    pr_auc = stat.mean(a)
    print(f'AUC: {aucs:.4f}\nPRAUC: {pr_auc:.4f}')

def compute_rank_roc(ranks, n_var):
    auc_x = list(ranks.keys())
    auc_x.sort()
    auc_y = []
    tpr = 0
    sum_rank = sum(ranks.values()) 
    for x in auc_x:
        tpr += ranks[x]
        auc_y.append(tpr / sum_rank)
    
    auc_x.append(n_var)
    auc_y.append(1)
    aucs = np.trapz(auc_y, auc_x) / n_var
    return aucs
        
def top_k(evaluation,k,total): 
     scores=[]
     top_k_result=[]
     for data in evaluation:
        if data <= k:
             top_k_result.append(data)
                         
     rank = len(top_k_result)/total
     print(f'Top@{k}: {len(top_k_result)} ({rank:.4f})')     
          
if __name__ == '__main__':
    main() 
    
