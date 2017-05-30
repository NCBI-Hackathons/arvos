#!/usr/bin/env python
"""

run the pipeline

   (1) Associative stats
   (2) Predictive stats

"""


## make imports
import os,sys,ast,csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn import metrics, svm
from sklearn.ensemble import RandomForestClassifier

sys.path.append(os.path.join("..","..","arvos"))
from Pipeline import Pipeline

results_dir = os.path.join(".","results")
pline = Pipeline(results_dir)

figures_dir = os.path.join(results_dir, "figures")

temp_figures_dir = os.path.join("..", "..", "web_app", "static", "img")
if not os.path.exists(figures_dir):
    print("creating a figures directory")
print(os.getcwd())
print(temp_figures_dir)

deseq_file = os.path.join(results_dir,"deseq.csv")
deseq_matrix_file = os.path.join(results_dir,"deseq-samples.csv")
targets_file = os.path.join(".","data","targets.csv")

X,y = pline.generate_features_and_targets(deseq_file,deseq_matrix_file,targets_file)


## use saved grid search parameters
log_file = "estimated-params.log"
if not os.path.exists(log_file):
    print("training data has been saved.  Run grid_search.py then run this function again")
    sys.exit()

print("loading saved parameters...")
with open(log_file, 'r') as fid:
    reader = csv.reader(fid)
    params = {key: ast.literal_eval(value) for (key,value) in reader}

#### Run the classifiers ####
print("...running cross validated model(s)")

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)


## run a svm alone
print("...running svm")
clf_svm = svm.SVC(probability=True,
                  kernel=params['svm']['kernel'],
                  C=params['svm']['C'],
                  gamma=params['svm']['gamma'])
y_score_svm = clf_svm.fit(X, y).decision_function(X)
y_pred_svm = clf_svm.predict(X)

## run a random forest
clf_random = RandomForestClassifier(n_estimators=20,max_depth=params["rf"]["max_depth"], max_features=params["rf"]["max_features"], min_samples_split=params["rf"]["min_samples_split"], min_samples_leaf=params["rf"]["min_samples_leaf"])
y_pred_random = clf_random.fit(X, y).predict(X)

#### Compute metrics for all estimators
## svm metrics
fpr_svm, tpr_svm, _ = metrics.roc_curve(y, y_score_svm)
svm_roc_auc = metrics.auc(fpr_svm, tpr_svm)
f1_svm = metrics.f1_score(y, y_pred_svm)

## random forest metrics
fpr_rand, tpr_rand, _ = metrics.roc_curve(y, y_pred_random)
random_roc_auc = metrics.auc(fpr_rand, tpr_rand)
f1_random = metrics.f1_score(y, y_pred_random)

#### plot roc curves

plt.title('Receiver Operating Characteristic for Random Forest and SVM')

plt.plot(fpr_svm,
         tpr_svm, 'b',
         color='cornflowerblue',
         label='SVM AUC: %0.2f'% svm_roc_auc)
plt.plot(fpr_rand,
         tpr_rand, 'b',
         color='navy',
         label='Random Forest AUC: %0.2f'% random_roc_auc)

plt.legend(loc='lower right')
plt.plot([0,1],[0,1],'r--')
plt.xlim([-0.1,1.2])
plt.ylim([-0.1,1.2])
plt.savefig(os.path.join(temp_figures_dir, "roc_curves.png"))
plt.show()

print("Done")
sys.exit()
