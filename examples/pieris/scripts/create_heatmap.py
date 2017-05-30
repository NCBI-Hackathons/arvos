#!/usr/bin/env python
"""
create the manuscript figure summarizing features

volcano plot - showing sig genes/isoforms
clustering heatmap - showing how genes/isoforms cluster


"""

import sys,csv,os,getopt,re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from htsint.tools import read_matrix,read_de_results,Heatmap
from htsint.blast import BlastMapper
from htsint.database import db_connect, Taxon, Gene
from sqlalchemy.sql import select

assembly = 'dn'
transcript = 'isoforms'

homeDir = os.path.join("data")
#os.path.join(os.path.expanduser("~"),"sequencing","pieris")

if transcript not in ['genes', 'isoforms']:
    raise Exception('Invalid transcript argument')
if assembly not in ['gg','dn']:
    raise Exception('Invalid assembly argument')

## variables
fontSize = 10
fontName = 'sans-serif'
useColor = True
if useColor == True:
    myCmap = mpl.cm.gist_heat
else:
    myCmap = mpl.cm.gray

## get some dictionaries form htsint
session,engine = db_connect()
conn = engine.connect()
s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(['7227']))
_taxaQueries = conn.execute(s)
taxaQueries = _taxaQueries.fetchall()
gene2taxa,gene2desc,gene2sym = {},{},{}
for tquery in taxaQueries:
    s = select([Gene.taxa_id,Gene.ncbi_id,Gene.description,Gene.symbol],Gene.taxa_id==tquery['id'])
    _geneQueries = conn.execute(s)
    geneQueries = _geneQueries.fetchall()
    gene2taxa.update(dict([(str(r['ncbi_id']),str(r['taxa_id'])) for r in geneQueries]))
    gene2desc.update(dict([(str(r['ncbi_id']),str(r['description'])) for r in geneQueries]))
    gene2sym.update(dict([(str(r['ncbi_id']),str(r['symbol'])) for r in geneQueries]))

## load a blast map                             
bm = BlastMapper()
summaryFile1 = os.path.join(homeDir,"dn-trinity",'blast-dn-parsed_summary.csv')
summaryFile2 = os.path.join(homeDir,"dn-trinity",'blast-dm-parsed_summary.csv')
summaryFile3 = os.path.join(homeDir,"dn-trinity","blast-mc-parsed_summary.csv")
summaryFile4 = os.path.join(homeDir,"dn-trinity",'blast-dp-parsed_summary.csv')

bmapSP = bm.load_summary(summaryFile1,trinityGene=False,best=True)
bmapDM = bm.load_summary(summaryFile2,trinityGene=False,best=True)
bmapMC = bm.load_summary(summaryFile3,trinityGene=False,best=True)
bmapDP = bm.load_summary(summaryFile4,trinityGene=False,best=True)

## load ref2gene 
reader = csv.reader(open("../gene2ref.tab","r"),delimiter="\t")
ref2gene = {}
for linja in reader:
    geneId = linja[1]
    proteinAccession = linja[5]
    ref2gene[proteinAccession] = geneId

## load feature data
featuresDir = os.path.join(homeDir,"features")
deseqResultsPath = os.path.join(featuresDir,"deseq.csv")
deseqIds, deseqColumns, deseqMat = read_de_results(deseqResultsPath,tool='DESeq')
dfeMatrixPath = os.path.join(featuresDir,"deseq-samples.csv")
dfeIds,dfeColumns,dfeMat = read_matrix(dfeMatrixPath,mtype='float')

## filter out nans
padjInd = np.where(deseqColumns == 'padj')[0]
size1 = deseqIds.shape[0]
filter1 = np.where(~np.isnan(deseqMat[:,padjInd]))[0]
deseqIds = deseqIds[filter1]
deseqMat = deseqMat[filter1,:]
mask = [np.where(dfeIds == i)[0][0] for i in deseqIds]
dfeIds = dfeIds[mask]
dfeMat = dfeMat[mask,:]
print("... %s/%s transcripts pass nan filter"%(filter1.size,size1))

## filter for only the most significant transcripts (max 75) 
threshold = 0.05
size2 = deseqIds.shape[0]

## print filter info
filter3 = np.where(deseqMat[:,padjInd] <= 0.05)[0]
filter4 = np.where(deseqMat[:,padjInd] <= 0.5)[0]
print("... %s/%s transcripts pass significance (%s) filter"%(filter3.size,size2,'0.05'))
print("... %s/%s transcripts pass significance (%s) filter"%(filter4.size,size2,'0.5'))

filter2 = np.where(deseqMat[:,padjInd] <= threshold)[0]#[:65]
deseqIds = deseqIds[filter2]
deseqMat = deseqMat[filter2,:]
mask = [np.where(dfeIds == i)[0][0] for i in deseqIds]
dfeIds = dfeIds[mask]
dfeMat = dfeMat[mask,:]

print("... %s/%s transcripts pass significance (%s) filter"%(filter2.size,size2,threshold))

## map any ids to gene names 
mappedIds = []
usedInds = []
debug =0 
for t,transcriptId in enumerate(dfeIds):
    bmap,geneId = None,None
    if bmapDM.has_key(transcriptId):
        bmap = bmapDM
    elif bmapSP.has_key(transcriptId):
        bmap = bmapSP
    elif bmapMC.has_key(transcriptId):
        bmap = bmapMC
    elif bmapDP.has_key(transcriptId):
        bmap = bmapDP
        
    if bmap:
        mId = bmap[transcriptId]
        debug += 1
        if mId[1] != '-':
            geneId = mId[1]
        elif ref2gene.has_key(mId[0]):
            geneId = ref2gene[mId[0]]
                        
    if geneId:
        usedInds.append(t)
        if gene2sym.has_key(geneId):
            mappedIds.append(gene2sym[geneId])
        else:
            mappedIds.append(geneId)
    elif bmap:
        usedInds.append(t)
        mappedIds.append(bmap[transcriptId][0])
    else:
        continue
        #mappedIds.append(transcriptId)

## filter again the data based on usedInds
filter3 = np.array(usedInds)
dfeIds = dfeIds[filter3]
dfeMat = dfeMat[filter3,:]
mask3 = [np.where(deseqIds == i)[0][0] for i in dfeIds]
deseqIds = deseqIds[mask3]
deseqMat = deseqMat[mask3,:]
print("... %s/%s transcripts pass BLAST filter"%(filter3.size,len(mappedIds)))

## create the heatmap
rowLabels = mappedIds
colLabels = dfeColumns
hm = Heatmap(dfeMat[:65],rowLabels[:65],colLabels,fontName=fontName,fontSize=fontSize,vpad=0.07)
hm.draw_heatmap(cmap='uy',clabels=True,rlabels=True,rowFont=6)
hm.save("heatmap-top-significant.png",dpi=500)

hm = Heatmap(dfeMat,rowLabels,colLabels,fontName=fontName,fontSize=fontSize,vpad=0.07)
hm.draw_heatmap(cmap='uy',clabels=True,rlabels=False,rowFont=6)
hm.save("heatmap-significant.png",dpi=500)
