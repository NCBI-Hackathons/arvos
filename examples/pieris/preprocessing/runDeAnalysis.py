#!/usr/bin/python
"""
Run DESeq and EdgeR
"""

import os,sys,re,csv
import numpy as np
from htsint import run_subprocess


class Associative(object):
    """
    A generic class
    """

    def __init__(self):
        """
        Constructor

        """

    def create_filtered(countsPath):
        if not os.path.exists(countsPath):
            raise Exception("Cannot find counts path")
    
        filteredCountsPath = re.sub("\.csv","-filtered.csv",countsPath)
        fid1 = open(countsPath,'r')
        fid2 = open(filteredCountsPath,'w')
        reader = csv.reader(fid1)
        writer = csv.writer(fid2)
        header = reader.next()
        writer.writerow(header)

        for linja in reader:
            if np.array([float(i) for i in linja[1:]]).sum() > 1:
                writer.writerow(linja)
        fid1.close()
        fid2.close()
        return filteredCountsPath

    def run_deseq(countsPath,outFile):
        cmd = "Rscript runDESeq.R %s %s"%(countsPath,outFile)
        print("running...\n%s"%cmd)
        run_subprocess(cmd)

if __name__ == "__main__":
    
    ## specify the locations
    homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
    readsDir = os.path.join(homeDir,'reads')

    featuresDir = os.path.join(homeDir,"features")
    countsPath = os.path.join(featuresDir,"est_counts.csv")
    filteredCountsPath = create_filtered(countsPath)
    outFile = os.path.join(featuresDir,"deseq.csv")
    run_deseq(filteredCountsPath,outFile)


        


