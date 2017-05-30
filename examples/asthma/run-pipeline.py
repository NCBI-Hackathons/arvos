#!/usr/bin/env python

print(__doc__)
import os,sys,csv
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import LogisticRegression
from sklearn.cross_validation import train_test_split
from sklearn.metrics import f1_score
from sklearn import metrics
from multiprocessing import Pool, cpu_count

outliers = True
invivo = False
runCombat = True

if outliers == True and invivo == True:
    tag = 'invivo-outlier'
elif outliers == True and invivo == False:
    tag = 'noninvivo-outlier'
elif outliers == False and invivo == True:
    tag = 'invivo-nonoutlier'
elif outliers == False and invivo == False:
    tag = 'noninvivo-nonoutlier'

compareFile = os.path.join("results","meta-analysis-%s.csv"%tag)

if outliers:
    expFile = os.path.join("data","merged-exp-o.csv")
else:
    expFile = os.path.join("data","merged-exp.csv")
    
if not os.path.exists(compareFile):
    raise Exception("cannot find compare file: use 'runMetaAnalysis.py' to produce \n %s"%compareFile)

if not os.path.exists(expFile):
    raise Exception("cannot find exp file: use 'createSummaryFiles.py' to produce \n %s"%expFile)

name2study = {"Christenson":"GSE67472",
              "Voraphani":"GSE43696",
              "Giovannini-Chami":"GSE19187",
              "Kicic":"GSE18965",
              "Yang":"GSE8190",
              "unpublished":"unpublished",
              "Wagener":"GSE51392",
              "yang16":"yang16"}

study2name = {v: k for k, v in name2study.iteritems()}

## specify studies
if invivo == False:
    _includedStudies = ['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'unpublished','Wagener', 'Christenson','yang16']
else:
    _includedStudies = ['Giovannini-Chami', 'Yang', 'Voraphani', 'unpublished', 'Christenson','yang16']

includedStudies = np.array([name2study[i] for i in _includedStudies])
compare = {}
fid = open(compareFile,"r")
reader = csv.reader(fid)
header = reader.next()
for h in header:
    compare[h] = []
for linja in reader:
    for i,l in enumerate(linja):
        compare[header[i]].append(l)
for key,item in compare.iteritems():

    if key not in ['gene', 'symbol']:
        compare[key] = np.array([float(i) for i in item])
    else:
        compare[key] = np.array(item)

compareInds = np.where(compare["Stouffer-adjusted-pvalue"] < 0.1)[0]
compareGenes = compare['gene'][compareInds]
print("there are %s sig indices in the compare file"%compareInds.size)
        
fid = open(expFile,'r')
reader = csv.reader(fid)
dx = np.array(reader.next()[2:])
studies = np.array(reader.next()[2:])
age = np.array(reader.next()[2:])
gender = np.array(reader.next()[2:])

geneSymbol,geneName = [],[]
for linja in reader:
    geneName.append(linja[0])
    geneSymbol.append(linja[1])
fid.close()
geneName = np.array(geneName)
geneSymbol = np.array(geneSymbol)

print("loading matrix")
mat = np.genfromtxt(expFile, delimiter=',',skip_header=4,usecols = range(2,len(linja)))
if mat.shape[1] != dx.size:
    raise Exception("dimension mismatch")
mat = mat.T

## use only includedStudies
includedInds = np.where(np.in1d(studies,includedStudies)==True)[0]

if includedInds.size != studies.size:
    print("subseting original data")
    print("\t %s x %s"%(mat.shape[0],mat.shape[1]))
    dx = dx[includedInds]
    studies = studies[includedInds]
    age = age[includedInds]
    gender = gender[includedInds]
    mat = mat[includedInds,:]
    print("\t %s x %s"%(mat.shape[0],mat.shape[1]))

## organize covariates
covs = {'gender':gender,
        'age':age,
        'dx':dx,
        'geo':studies,
        'study':np.array([study2name[s] for s in studies])}

## remove genes with nans
print("removing nans...")
toKeep = np.where(np.isnan(mat.mean(axis=0)) == False)[0]
print("\t %s x %s"%(mat.shape[0],mat.shape[1]))
mat = mat[:,toKeep]

print("\t %s x %s"%(mat.shape[0],mat.shape[1]))
geneName = geneName[toKeep]
geneSymbol = geneSymbol[toKeep]

if geneName.size != mat.shape[1] or geneSymbol.size != mat.shape[1]:
    raise Exception("invalid gene list size")

## global standardiztion
mat = preprocessing.scale(mat)

## run combat

combatFile = "mat_combat-%s.csv"%tag
if not os.path.exists(combatFile):
    runCombat = True

print("DEBUG")
runCombat = True
    
if runCombat:
    np.savetxt("mat.csv",mat,delimiter=",")

    fidout = open("covs.csv","w")
    writer = csv.writer(fidout)
    writer.writerow(["study","geo","gender","age","dx"])
    for i in range(covs['study'].size):
        writer.writerow([covs[c][i] for c in ["study","geo","gender","age","dx"]])
    fidout.close()
              
    fidout = open("genes.csv","w")
    writer = csv.writer(fidout)
    writer.writerow(["gene","symbol"])
    for i,gene in enumerate(geneName):
        writer.writerow([gene,geneSymbol[i]])
    fidout.close()

    print('running combat')
    os.system("Rscript runCombat.R")

    os.system("mv mat.csv mat-%s.csv"%tag)
    os.system("mv mat_combat.csv mat_combat-%s.csv"%tag)
    os.system("mv covs.csv covs-%s.csv"%tag)
    os.system("mv genes.csv genes-%s.csv"%tag)

    
## load combat matrix
cmat = np.genfromtxt(combatFile,delimiter=',')
cmat = cmat.T

if mat.shape[1] != cmat.shape[1] or mat.shape[0] != cmat.shape[0]:
    raise Exception("dimension mismatch")

for study in name2study.keys():
    print("%s - %s"%(study,np.where(studies==name2study[study])[0].size))


## debug
print(cmat.shape)
sys.exit()



    
############################################################################################
## split training test
print("good")
def split_train_test(_testStudy):
    testStudy = name2study[_testStudy]
    if testStudy not in studies:
        raise Exception("invalid test study")
    
    testInds = np.where(studies==testStudy)[0]
    trainInds = np.where(studies!=testStudy)[0]

    dx_test = dx[testInds]
    dx_train = dx[trainInds]

    studies_test = studies[testInds]
    studies_train = studies[trainInds]

    age_test = age[testInds]
    age_train = age[trainInds]

    gender_test = gender[testInds]
    gender_train = gender[trainInds]

    cmat_test = cmat[testInds,:]
    cmat_train = mat[trainInds,:]

    train = {'inds':trainInds,'studies':studies_train,'dx':dx_train,'age':age_train,'gender':gender_train,'cmat':cmat_train}
    test = {'inds':testInds,'studies':studies_test,'dx':dx_test,'age':age_test,'gender':gender_test,'cmat':cmat_test}
    
    return train,test

############################################################################################
## setup features and targets
"""
The parameter l1_ratio corresponds to alpha in the glmnet R package while alpha corresponds to the lambda parameter in glmnet.
"""

def run_enet(args):

    train,test,_alpha,_lambda,ax,color,label,legend,legKeys = args
    y_train = np.zeros(train['dx'].size)
    y_train[np.where(train['dx'] == 'asthma')[0]] = 1
    y_test = np.zeros(test['dx'].size)
    y_test[np.where(test['dx'] == 'asthma')[0]] = 1

    X_train = train['cmat']
    X_test = test['cmat']

    classifier = ElasticNet(alpha=_lambda, l1_ratio=_alpha,max_iter=5000)
    #probas_ = classifier.fit(X_train, y_train).predict(X_test)
    y_pred = classifier.fit(X_train, y_train).predict(X_test)
    
    fpr, tpr, thresholds = metrics.roc_curve(y_test,y_pred,pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    print('lambda', _lambda, 'alpha', _alpha, 'auc', roc_auc)

    p = ax.plot(fpr,tpr, label='ROC curve (area = %0.2f)' % roc_auc)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('Receiver operating characteristic example')
     
    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC curve')
    legKeys.append((p[0],'%s - (%0.2f)' %(label,roc_auc)))
    if legend:
        l = ax.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')
        legKeys.append((l[0],'Random'))
        l1, l2 = zip(*legKeys)
        leg = ax.legend(l1,l2,loc="lower right")

        ltext = leg.get_texts()
        for i in range(len(ltext)):
            plt.setp(ltext[i],fontsize=10) 

    return classifier.coef_
        
    #fscore = f1_score(y_test, y_pred, average='binary')
    #print fscore

####################################################################################
## specify parameters determined from cross-validation
_alpha = 0.9
_lambda = 0.05

fig = plt.figure()
ax = fig.add_subplot(111)
legKeys = []
tstudy = 'unpublished'
train, test = split_train_test(tstudy)
args = (train,test,_alpha,_lambda,ax,'k','initial',True,legKeys)
coefs = run_enet(args)
plt.savefig("roc.png",dpi=400)
    
print(coefs.shape)
geneInds = np.where(np.abs(coefs) != 0.0)[0]
#sys.exit()

fid = open("enet-genes-coef.csv","w")
writer = csv.writer(fid)
writer.writerow(["gene","symbol","coef"])
for indx in geneInds:
    writer.writerow([[geneName[indx]][0],geneSymbol[indx],coefs[indx]])

print(np.intersect1d(geneName[geneInds],compareGenes))
    
sys.exit()
print("classification report:")
print(metrics.classification_report(y_test, y_pred, target_names=covs['dx']))



