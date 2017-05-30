#!/usr/bin/Rscript
## install notes
# > source("http://bioconductor.org/biocLite.R")
# > biocLite("DESeq2")

## samples notes
# order:  17,18,33,46,56,61,DL47,DL61,D163,D178,D185,D239
#        c('R','D','D','R','R','D','R','D','D','R','D','R')
#old
# D: 17,46,56,DL61,D163,D185
# R: 18,33,61,DL47,D178,D239

# D: 18,33,46,56,DL61,D163,D185
# R: 17,33,61,DL47,D178,D239


# same as experiment 1 except we are filtering the genes
suppressMessages(library("DESeq2"))
args <-commandArgs(TRUE)

if (length(args)==0){
    print("ERROR: Did not specify counts file e.g 'Rscript runDESeq.R dn-genes-counts.csv'")
    q()
}

countsFile = args[1]
outFile = args[2]

if (!file.exists(countsFile)){
    print("ERROR: invalid counts file")
    q()
}

#### from counts file (count matrix)
countData <- read.csv(countsFile,header=TRUE,row.names=1,com='')
countData <- round(countData)
colData <- data.frame(
    row.names=colnames(countData),
    condition=c("D", "R", "R", "D", "D", "R","R","D","D","R","D","R"),
    libType=c("paired", "paired", "paired","paired","paired","paired","paired","paired","paired","paired","paired","paired")
)

## Get the normalized counts to perform filtering
dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
dds$condition <- factor(dds$condition,levels=c("R","D"))
dds <- DESeq(dds)
rld <- rlog(dds, blind=TRUE)
res <- results(dds)
use <- res$baseMean>metadata(res)$filterThreshold
print(table(use))

## run DESeq2 using filtered data
countData <- countData[use& ! is.na(res$pvalue),]
colData <- data.frame(
    row.names=colnames(countData),
    condition=c("D", "R", "R", "D", "D", "R","R","D","D","R","D","R"),
    libType=c("paired", "paired", "paired","paired","paired","paired","paired","paired","paired","paired","paired","paired")
)

dds <- DESeqDataSetFromMatrix(countData=countData,colData=colData,design=~condition)
dds$condition <- factor(dds$condition,levels=c("R","D"))
dds <- DESeq(dds)
rld <- rlog(dds, blind=TRUE)
res <- results(dds)
resOrdered <-res[order(res$padj),]
print(head(resOrdered))

## export results to csv
write.csv(as.data.frame(resOrdered),file=outFile,quote=FALSE)
write.csv(as.data.frame(assay(rld)),file=gsub("\\.csv","-samples.csv",outFile),quote=FALSE)
