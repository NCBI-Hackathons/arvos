#!/usr/bin/python
"""
grid search
parameter output can be hardcoded into run.py and into the notebook
"""

import csv,time,sys,os
import numpy as np

from sklearn.model_selection import GridSearchCV,RandomizedSearchCV
from sklearn import preprocessing,tree,svm
from sklearn.ensemble import RandomForestClassifier

#from mylib import *

METRIC = 'f1_micro'

def main():
    """
    grid search for churn study
    """

    #######################################
    ## load data
    #######################################
    time_start = time.time()


    sys.path.append(os.path.join("..","..","arvos"))
    from Pipeline import Pipeline

    results_dir = os.path.join(".","results")
    pline = Pipeline(results_dir)

    deseq_file = os.path.join(results_dir,"deseq.csv")
    deseq_matrix_file = os.path.join(results_dir,"deseq-samples.csv")
    targets_file = os.path.join(".","data","targets.csv")

    X,y = pline.generate_features_and_targets(deseq_file,deseq_matrix_file,targets_file)

    features = pline.dfe_genes

    print(features[:10])

    #file_name = 'train.npz'    
    #if not os.path.exists(file_name):
    #    raise Exception("you must save the training data before grid_search")
    #print("...loading saved data")
    #file_name = 'train.npz'
    #npz = np.load(file_name)
    #X = npz['X_train']
    #y = npz['y_train']

    fid = open("estimated-params.log","w")
    log_file = csv.writer(fid)
    
    ###################################
    ## run random forest
    #######################################
    print("...running random forest grid search")
    rf = RandomForestClassifier(n_estimators=20)
    parameters = {"max_depth": [3, None],
                  "max_features": [1, 3, 10],
                  "min_samples_split": [2, 5, 10],
                  "min_samples_leaf": [1, 3, 10],
                  "bootstrap": [True, False],
                  "criterion": ["gini", "entropy"]}

    clf = GridSearchCV(rf,parameters,n_jobs=-1,verbose=1,scoring=METRIC)
    clf.fit(X, y)
    log_file.writerow(["rf",str(clf.best_params_)])

    #######################################
    ## run svm
    #######################################
    print("...running svm grid search (this will take some time")
    subsample = np.arange(y.size)
    np.random.shuffle(subsample)
    subsample = subsample[:15000]
    parameters = {'kernel':('linear', 'rbf'),
                  'gamma': [1e-02, 5e-02],
                  'C':[0.1,10.0,100.0,200.0,300.0]}
    svr = svm.SVC()
    clf = GridSearchCV(svr, parameters,n_jobs=-1,verbose=1,scoring=METRIC)
    clf.fit(X[subsample,:], y[subsample])
    log_file.writerow(["svm",str(clf.best_params_)])

    ## print run time
    run_time = time.strftime('%H:%M:%S', time.gmtime(time.time()-time_start))
    fid.close()
    print("\nTOTAL RUN TIME (HH:MM:SS):%s"%run_time)
    
if __name__ == "__main__":
    main()
