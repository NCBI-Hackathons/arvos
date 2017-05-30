#!/usr/bin/python
"""
Run STAR to align reads to transcriptome
"""

import os,sys,re,shutil,getopt,itertools
import HTSeq
from htsint import run_subprocess

## locations
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
sailfishDir =  os.path.join(homeDir,'sailfish')
readsDir = os.path.join(homeDir,'reads')
genomePath = os.path.realpath(os.path.join(homeDir,'dn-trinity','Trinity.fasta'))
#gff3Path  = os.path.realpath(os.path.join(".","Xentr7_2_Stable.gff3"))
#gtfPath  = os.path.realpath(os.path.join(sailfishDir,"Xentr7_2_Stable.gtf"))
sailfishPath = "/usr/src/Sailfish-0.6.3-Linux_x86-64/bin/sailfish"

def get_reads(sampleList):

    def get_trimmed_files(sample):
        leftFile,rightFile = None,None
        for fileName in os.listdir(readsDir):
            if re.search("unpaired",fileName):
                continue

            if re.search("^%s.*\.fq$"%sample,fileName) and re.search("left_paired",fileName):
                leftFile = os.path.realpath(os.path.join(readsDir,fileName))
            elif re.search("^%s.*.fq$"%sample,fileName) and re.search("right_paired",fileName):
                rightFile = os.path.realpath(os.path.join(readsDir,fileName))
        return leftFile,rightFile

    allLeft = []
    allRight = []
    for sample in sampleList:
        left,right = get_trimmed_files(sample)
        allLeft.append(left)
        allRight.append(right)

    ## check
    if len(allLeft) != len(sampleList):
        raise Exception("Invalid number of left sequences")
    if len(allRight) != len(sampleList):
        raise Exception("Invalid number of right sequences")

    return allLeft,allRight

if __name__ == "__main__":

    sampleList = ["17", "18", "33", "46", "56", "61", "DL47", "DL61", "D163", "D178", "D185", "D239"]

    ## export sutff
    print("export LD_LIBRARY_PATH=/usr/src/Sailfish-0.6.3-Linux_x86-64/lib:$LD_LIBRARY_PATH")
    print("export PATH=/usr/src/Sailfish-0.6.3-Linux_x86-64/bin:$PATH")

    ## copy genome to star dir
    cmd1 = "cp %s %s"%(genomePath,os.path.join(sailfishDir,'genome.fa'))
    print("\n...copy genome to sailfish dir")
    print(cmd1)
        
    ## create the index
    cmd2 = "%s index -t %s -o %s -k 20"%(sailfishPath,os.path.join(sailfishDir,'genome.fa'),sailfishDir)
    print("\n...index the genome")
    print(cmd2)

    ## align reads
    print("\n...align reads")
    allLeft,allRight = get_reads(sampleList)
    cmd3 = ""
    for s,sample in enumerate(sampleList):
        outDir = os.path.join(homeDir,'features',sample)
        if not os.path.exists(outDir):
            os.mkdir(outDir)

        if s != len(sampleList) -1:
            _cmd = '%s quant -i %s -o %s -l "T=PE:O=><:S=SA" '%(sailfishPath,sailfishDir,outDir)+\
                   '-1 %s -2 %s && '%(allLeft[s],allRight[s])
        else:
            _cmd = '%s quant -i %s -o %s -l "T=PE:O=><:S=SA" '%(sailfishPath,sailfishDir,outDir)+\
                   '-1 %s -2 %s'%(allLeft[s],allRight[s])
       
        cmd3 += _cmd
    print cmd3

    overwrite = False
    sys.exit()



    ## convert the sam files to coordinate-sorted bam and files
    print("\n...sam/bam conversions and sort")
    for s,sample in enumerate(sampleList):
        samFile = os.path.join(sailfishDir,"%s_Aligned.out.sam"%(sample))
        bamFile = os.path.join(readsDir,"%s_aligned.bam"%(sample))
        sbamFile = os.path.join(readsDir,"%s_aligned_sorted.bam"%(sample))
        ssamFile = os.path.join(readsDir,"%s_aligned_sorted.sam"%(sample))
        if not os.path.exists(samFile):
            raise Exception("cannot find sam file %s"%(samFile))

        if os.path.exists(ssamFile) and not overwrite:
            print("skipping sam to bam, align, bam to sam")
        else:
            cmd = "/usr/bin/samtools view -b -S %s > %s && "%(samFile,bamFile)
            cmd += "/usr/bin/samtools sort -n %s %s && "%(bamFile,sbamFile[:-4])
            cmd += "/usr/bin/samtools view -h %s > %s"%(sbamFile,ssamFile)
            print cmd
            run_subprocess(cmd)

    ## concat sam files and sort by coordinates
    print("\n...make single sorted bam file..")
    outBam = os.path.join(homeDir,"star_all_reads.bam")
    outSbam = os.path.join(homeDir,"star_all_reads_sorted.bam")
    cmdMerge = "/usr/bin/samtools merge %s"%(outBam)
    cmdSort = "/usr/bin/samtools sort %s %s"%(outBam,outSbam[:-4])
    for s,sample in enumerate(sampleList):
        sbamFile = os.path.join(readsDir,"%s_aligned_sorted.bam"%(sample))
        cmdMerge += " %s"%(sbamFile)
        
    cmdMergeSort = cmdMerge + " && %s"%(cmdSort)

    if os.path.exists(outSbam) and not overwrite:
        print("skipping concat bam files and sort")
    else:
        if os.path.exists(outSbam):
            os.remove(outSbam)
        if os.path.exists(outBam):
            os.remove(outBam)
        print cmdMergeSort
        run_subprocess(cmdMergeSort)
    print("\n")
