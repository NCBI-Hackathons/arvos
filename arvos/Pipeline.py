#!/usr/bin/python
"""
Run DESeq and EdgeR


Class contains:

   General Methods
      1. 
      2. 

   Associative Methods
      1. 

"""

import os,sys,re,csv
import numpy as np
import pandas as pd
import sklearn
import subprocess

from htsint.tools import read_matrix,read_de_results


class Pipeline(object):
    """
    A pipeline
    """

    def __init__(self,results_dir):
        """
        Constructor
        """

        if not os.path.isdir(results_dir):
            print("...creating results dir")
            os.mkdir(results_dir)
        
        self.results_dir = results_dir
        

    def run_subprocess(self,cmd):
        proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stdin=subprocess.PIPE)
        try:
            outs, errs = proc.communicate(timeout=22000)
        except TimeoutExpired:
            proc.kill()
            outs, errs = proc.communicate()

    def generate_features_and_targets(self,deseq_file,deseq_matrix_file,targets_file):

        """
        input: deseq_file - ouptut from run_deseq
               deseq_matrix_file - ouptut from run_deseq
               targets_file - csv one row with a comma seperated classes (no header)
        
        samples are rows, genes are columns
        """

        ## error checking
        for fname in [deseq_file,deseq_matrix_file,targets_file]:
            if not os.path.exists(fname):
                raise Exception("file does not exist: %s"%fname)
        
        ## read in files
        deseqIds, deseqColumns, deseqMat = read_de_results(deseq_file,tool='DESeq')
        dfeIds,dfeColumns,X = read_matrix(deseq_matrix_file,mtype='float')
        y = np.loadtxt(targets_file,delimiter=",")
        print(y)
        
        ## save columns and rows
        self.deseqIds, self.deseqColumns = deseqIds, deseqColumns
        self.dfe_genes, self.dfe_samples =  dfeIds, dfeColumns
        self.deseqMat = deseqMat
        
        return X.transpose(),y


    #######################################
    # Methods for the Associative portion #
    #######################################
    def create_filtered(self, countsPath):
        if not os.path.exists(countsPath):
            raise Exception("Cannot find counts path")

        filteredCountsPath = re.sub("\.csv","-filtered.csv",countsPath)
        fid1 = open(countsPath,'r')
        fid2 = open(filteredCountsPath,'w')
        reader = csv.reader(fid1)
        writer = csv.writer(fid2)
        header = next(reader)
        writer.writerow(header)

        for linja in reader:
            if np.array([float(i) for i in linja[1:]]).sum() > 1:
                writer.writerow(linja)
        fid1.close()
        fid2.close()
        return filteredCountsPath

    def run_deseq(self,countsPath,outFile):
        cmd = "Rscript runDESeq.R %s %s"%(countsPath,outFile)
        print("running...\n%s"%cmd)
        self.run_subprocess(cmd)

    # input:
    # (self, header array from count file, path (without file name))
    # output:
    # saved file for user to do something with.
    def two_group(self, head, c_Path):
        #is outputting a csv and having them manipulate good?
        with open(c_Path+ 'vartypes.csv', 'w',newline='\n') as f:
            writer = csv.writer(f)
            writer.writerow(head)
            writer.writerow(['variable']+ ['0 or 1']*(len(head)-1))


if __name__ == "__main__":

    ## specify the locations
    homeDir = os.path.join("..")
    readsDir = os.path.join(homeDir, "reads")

    parentDir = os.path.join(homeDir,"examples", "pieris")
    countsPath = os.path.join(parentDir,"data", "est_counts.csv")
    outDirectory = os.path.join(parentDir, "results")
    outFile = os.path.join(outDirectory, "deseq.csv")

    # check if results directory exists, if !exist create one
    if not os.path.exists(outDirectory):
        os.mkdir(outDirectory)


    associative = Associative()
    filteredCountsPath = associative.create_filtered(countsPath)
    associative.run_deseq(filteredCountsPath,outFile)
