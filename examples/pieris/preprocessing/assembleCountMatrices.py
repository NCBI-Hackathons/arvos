#!/usr/bin/python
"""
Use normalized count matrix provided by DESeq2 to get transformed values
"""

__author__ = "Adam Richards"

import sys,os,getopt,csv,time,re,gc
import numpy as np

## declare variables
#homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")

def load_file(sample):

    ## error check
    resultsFile = os.path.join(homeDir,'features',sample,'quant.sf')
    if not os.path.exists(resultsFile):
        raise Exception("Cannot find results file %s"%resultsFile)

    ## infile
    fidin = open(resultsFile,'r')
    reader = csv.reader(fidin,delimiter="\t")

    debug = 0
    header = ['Transcript','Length','TPM','RPKM','KPKM','EstimatedNumKmers','EstimatedNumReads']
    results = {}

    for key in header:
        results[key] = []

    gc.disable()
    for linja in reader:
        if linja[0][0] == '#':
            continue
        
        results['Transcript'].append(linja[0])
        results['TPM'].append(linja[2])
        results['RPKM'].append(linja[3])
        results['EstimatedNumReads'].append(linja[6])

    gc.enable()
    return results

if __name__ == "__main__":

    sampleList = ["17", "18", "33", "46", "56", "61", "DL47", "DL61", "D163", "D178", "D185", "D239"]

    allResults = {}
    for sample in sampleList:
        allResults[sample] = load_file(sample)
    
    def write_matrix(column):
        numTranscripts = len(allResults['17']['Transcript']) 
        
        ## create a count matrix
        if column == 'EstimatedNumReads':
            outFile = os.path.join(homeDir,'features','est_counts.csv')
        else: 
            outFile = os.path.join(homeDir,'features','%s_counts.csv'%(column.lower()))
        fidout = open(outFile,'w')
        writer = csv.writer(fidout)
        writer.writerow(['transcript'] + sampleList)
        
        for row in range(numTranscripts):
            trans = allResults['17']['Transcript'][row]
            toWrite = [allResults[s][column][row] for s in sampleList]
            writer.writerow([trans] + toWrite)
    
    write_matrix('EstimatedNumReads')
    write_matrix('TPM')
    write_matrix('RPKM')

print('complete.')
