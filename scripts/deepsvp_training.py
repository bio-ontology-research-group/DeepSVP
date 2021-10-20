# -*- coding: utf-8 -*-
import sklearn
import os
import shutil
from sklearn.model_selection import KFold
import glob
from sklearn.utils import shuffle
import math
import pandas as pd
import time
import gensim
import numpy as np
seed = 1
np.random.seed(seed)
import random
from sklearn import metrics
import tensorflow as tf
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
from tensorflow.python.ops.resource_variable_ops import ResourceVariable
r = ResourceVariable([5.20])
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.layers import (Input,InputLayer, Dense, Embedding, Conv1D, Flatten, Concatenate, Dropout, Dot, LeakyReLU,Activation)
from tensorflow.keras.optimizers import Adam, RMSprop
from sklearn.model_selection import StratifiedKFold,StratifiedShuffleSplit, train_test_split
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint
from tensorflow.keras.models import Sequential, Model
from sklearn.metrics import (classification_report, confusion_matrix, precision_score,auc,precision_recall_curve,average_precision_score,
                            recall_score, f1_score,balanced_accuracy_score,accuracy_score,roc_auc_score,roc_curve)
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
from scipy import stats
print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))
import skopt
from skopt import gp_minimize, forest_minimize
from skopt.space import Real, Categorical, Integer
from skopt.plots import plot_convergence
from skopt.plots import plot_objective, plot_evaluations
from skopt.plots import plot_histogram, plot_objective_2D
from skopt.utils import use_named_argss
import sys
import pickle as pkl
from skopt import dump, load
import itertools

ontology=sys.argv[1]
operation=sys.argv[2]

# Hyperparamter 
dim_learning_rate = Real(low=1e-6, high=1e-2, prior='log-uniform', name='learning_rate')
dim_num_dense_layers = Integer(low=2, high=6, name='num_dense_layers')
dim_num_dense_nodes = Integer(low=50, high=512, name='num_dense_nodes')
dim_activation = Categorical(categories=['relu', 'sigmoid','selu'], name='activation')
dim_dropout = Real(low=0.1, high=0.5, name='dropout')
dimensions = [dim_learning_rate,
              dim_num_dense_nodes,
              dim_num_dense_layers,
              dim_activation,
              dim_dropout]
default_parameters = [1e-5, 100, 2, 'relu',0.2]
path_best_model = 'models/'+onto+'_'+ope+'.h5' 
W_path_best_model = 'models/weights_'+onto+'_'+ope+'.h5'    

#Helper-function for log-dir-name
def log_dir_name(learning_rate,num_dense_nodes,dim_num_dense_layers, activation,dropout):
    # The dir-name for the TensorBoard log-dir.
    s = "logs/"+ontology+operation+"_lr_{0:.0e}_layers_nodes_{1}_{2}_{3}_{4}/"
    # Insert all the hyper-parameters in the dir-name.
    log_dir = s.format(learning_rate, num_dense_nodes,dim_num_dense_layers,activation,dropout)
    return log_dir

def load_data(onto, type_): 
    X_train = pd.read_csv("data/"+onto+type_+"training.tsv", sep='\t',low_memory=False)    
    negatives = X_train[X_train['label']==0]
    postives = X_train[X_train['label']==1]
    return negatives, postives, X_train

negatives , postives , X_train  = load_data(ontology, operation)
data_shape = postives.shape[1]-6
    
def create_model(learning_rate,num_dense_nodes, num_dense_layers, activation,dropout):
    """
    Hyper-parameters:
    learning_rate:     Learning-rate for the optimizer.
    num_dense_nodes:   Number of nodes in each dense layer.
    activation:        Activation function for all layers.
    """  
    # Start construction of a Keras Sequential model.
    model = Sequential()

    # Add an input layer which is similar to a feed_dict in TensorFlow.
    model.add(InputLayer(input_shape=(data_shape,)))
    
    # Add fully-connected / dense layers. The number of nodes and the activation function.
    for i in range(num_dense_layers):
        name = 'layer_dense_{0}'.format(i+1)
        # Add the dense / fully-connected layer to the model.
        model.add(tf.keras.layers.Dropout(dropout))
        model.add(Dense(num_dense_nodes,activation=activation,name=name))    
        model.add(tf.keras.layers.Dropout(dropout))
                       
    # Last fully-connected
    model.add(Dense(1, activation='sigmoid')) 
    
    # Use the Adam method for training the network, find the best learning-rate for the Adam method.
    optimizer = Adam(lr=learning_rate)
    early_stop = EarlyStopping(monitor='val_loss', patience=10)
    
    # compile the model so it can be trained.
    model.compile(optimizer=optimizer,
                  loss='binary_crossentropy',
                  metrics=['accuracy'])   
                  
    return model
    
best_accuracy = 0.0
@use_named_args(dimensions=dimensions)
def fitness(learning_rate,num_dense_nodes, num_dense_layers ,activation,dropout):
    negatives , postives , X_train  = load_data(onto, ope)
    print(negatives.shape)
    """
    Hyper-parameters:
    learning_rate:     Learning-rate for the optimizer.
    num_dense_nodes:   Number of nodes in each dense layer.
    activation:        Activation function for all layers.
    num_dense_nodes    Number of dense nodes.
    dropout
    """
    # Print the hyper-parameters.
    print('learning rate: {0:.1e}'.format(learning_rate))
    print('num_dense_layers:', num_dense_layers)
    print('num_dense_nodes:', num_dense_nodes)
    print('activation:', activation)
    print('dropout:', dropout)
    
    # Create the neural network with these hyper-parameters.
    model = create_model(learning_rate=learning_rate,
                         num_dense_nodes=num_dense_nodes,
                         num_dense_layers=num_dense_layers,
                         activation=activation,
                         dropout=dropout)
    # Dir-name for the TensorBoard log-files.
    log_dir = log_dir_name(learning_rate,num_dense_nodes, num_dense_layers,activation, dropout)
    
    # callback-function for Keras which will be run after each epoch has ended during training
    my_callbacks = [
    EarlyStopping(monitor='val_loss', patience=10),
    TensorBoard(log_dir=log_dir,
                histogram_freq=0,
                write_graph=True,
                write_grads=False,
                write_images=False)
    ]
    folds = 5 
    cvscores=[]
    for i in range(folds):
        n = negatives.sample(n=postives.shape[0])       
        d = pd.concat([n,postives])
        d = shuffle(d)
        targetdata = np.asarray(d['label'])
        traindata = np.asarray(d.iloc[:, 6:d.shape[1]])
        
        history = model.fit(x=traindata,
                                y=targetdata,
                                epochs=100,
                                batch_size=128,
                                validation_split=0.15,
                                callbacks=my_callbacks)   
        accuracy = history.history['val_acc'][-1]
            
        cvscores.append(accuracy)         
        negatives = negatives.drop(n.index)        
        
    # Print the classification accuracy.
    print()
    print("%.2f%% (+/- %.2f%%)" % (np.mean(cvscores), np.std(cvscores)))
    accuracy = np.mean(cvscores)

    # Save the model if it improves on the best-found performance.
    global best_accuracy
    # If the classification accuracy of the saved model is improved
    if accuracy > best_accuracy:        
        model.save(path_best_model, include_optimizer=True)
        model.save_weights(W_path_best_model)
        print('Save model')
        
        # Update the classification accuracy.
        best_accuracy = accuracy

    # Delete the Keras model with these hyper-parameters from memory.
    del model
    print("Done 1")
    
    return -accuracy 

def plot_figures(search_result):
    plot_convergence(search_result)
    plt.savefig("figures/"+operation+"_plot_convergence.png",dpi=1200)
    space = search_result.space
    print(search_result.fun)
    fig, ax = plot_histogram(result=search_result, dimension_name='activation')
    plt.savefig("figures/"+operation+"_plot_histogram.png",dpi=1200)
    fig = plot_objective_2D(result=search_result,
                        dimension_name1='learning_rate',
                        dimension_name2='num_dense_layers',
                        levels=50)
    plt.savefig("figures/"+onto+"_plot_objective_2D.png",dpi=1200)
    dim_names = ['learning_rate', 'num_dense_nodes', 'num_dense_layers']
    fig, ax = plot_objective(result=search_result, dimension_names=dim_names)
    plt.savefig("figures/"+operation+"_plot_objective.png",dpi=1200)
    fig, ax = plot_evaluations(result=search_result, dimension_names=dim_names)
    plt.savefig("figures/"+operation+"_plot_evaluations.png",dpi=1200)

#ref: http://www.anaesthetist.com/mnm/stats/roc/Findex.htm
def standard_error(score,TP,TN) :
    q1 = score/(2-score)
    a2 = score**2
    q2 = (2*a2)/(1+score) 
    c1 = (score*(1-score))
    c2 = ((TP-1)*(q1-a2)) 
    c3 = ((TN -1)*(q2-a2))
    c4 = (TP*TN)
    ci =  math.sqrt(( c1 + c2 + c3 )/c4)
    return (ci)

def evaluate(path_best_model, X_test, space):
    original_stdout = sys.stdout # Save a reference to the original standard output
    allresults = open("results/"+ontology+operation+"_results.txt", "a")
    sys.stdout = allresults # Change the standard output to the file we created.        
    print(ontology+" : "+operation)
    print("The best parameteres: "+str(search_result.x))
    print(space.point_to_dict(search_result.x))
    print(search_result.fun)
    
    model = load_model(path_best_model)
    test = X_test.iloc[:, 6:X_test.shape[1]]
    result = model.evaluate(np.asarray(test), X_test['label'].values,verbose=0)
    for name, value in zip(model.metrics_names, result):
        print(name, value)
    print("{0}: {1:.2%}".format(model.metrics_names[1], result[1]))
    
    pred = model.predict(np.asarray(test))    
    fpr, tpr, thresholds = roc_curve(X_test['label'], pred)
    confusion = metrics.confusion_matrix(X_test['label'], pred.round())
    TP = confusion[1, 1]
    TN = confusion[0, 0]
    FP = confusion[0, 1]
    FN = confusion[1, 0]
    print("TP ", TP , 'TN', TN, 'FN', FN , 'FP', FP)

    average_precision = average_precision_score(X_test['label'], pred)
    # calculate precision and recall for each threshold
    lr_precision, lr_recall, _ = precision_recall_curve(X_test['label'], pred)

    # calculate scores
    merged_pred = list(itertools.chain(*pred))
    round_to_whole = [round(num) for num in merged_pred]
    lr_f1 = f1_score(X_test['label'].values,pred.round())
    se =  standard_error(lr_f1,TP,TN)  
    print('F1 score: {a:0.4f} and the Error {b:0.4f}'.format(a=lr_f1,b=se))
        
    # calculate the precision-recall auc
    precision, recall, _ = precision_recall_curve(X_test['label'].values, pred)
    auc_score = auc(recall, precision)
    se =  standard_error(auc_score,TP,TN)  
    print('PR AUC: {a:0.4f} and the Error {b:0.4f}'.format(a=auc_score,b=se))
        
    # AUC
    AUC = roc_auc_score(X_test['label'], pred, average='macro', max_fpr=None) 
    se =  standard_error(AUC,TP,TN)   
    print('AUC: {a:0.4f} and the Error {b:0.4f}'.format(a=AUC,b=se))

    #The diagnostic odds ratio is defined mathematically as: https://en.wikipedia.org/wiki/Diagnostic_odds_ratio
    DOR = (TP/FN)/(FP/TN)
    se = math.sqrt( (1/TP) + (1/FN) + (1/FP) + (1/TN) )
    print('DOR: {a:0.4f} and the error {b:0.4f}'.format(a=DOR,b=se))
    print('\n')
    print('End-------------------------------------------------')
    sys.stdout = original_stdout # Reset the standard output to its original value
    allresults.close()


if __name__ == '__main__':
                        
    search_result = gp_minimize(func=fitness,
                                dimensions=dimensions,
                                acq_func='EI', # Expected Improvement.
                                n_calls=15,
                                x0=default_parameters)                                
    space = search_result.space
    plot_figures(search_result)
    dump(search_result, 'search_result_opt.pkl')
    
    # Test the model using the test dataset and report the perforsmenace 
    test_data = pd.read_csv("data/"+ontology+operation+"testing.tsv", sep='\t', low_memory=False)
    evaluate(path_best_model, test_data, space)
    
