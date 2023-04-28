# Load libraries (Download from BiocManager if needed)

library(BiocManager)
library(MotifDb)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(ATACseqQC)

## Put sorted BAM files into variables and make them into genomic alignments object ##
#### Indexing is only needed if you need an index file for the function (ex. estimateLibComplexity)

TS_01_sorted <- file.choose()
idx_TS_01 <- indexBam(TS_01_sorted)
GA_TS_01 <- readGAlignments(TS_01_sorted, param=ScanBamParam(what="flag"))

TS_02_sorted <- file.choose()
idx_TS_02 <- indexBam(TS_02_sorted)
GA_TS_02 <- readGAlignments(TS_02_sorted, param=ScanBamParam(what="flag"))

TS_03_sorted <- file.choose()
idx_TS_03 <- indexBam(TS_03_sorted)
GA_TS_03 <- readGAlignments(TS_03_sorted, param=ScanBamParam(what="flag"))

TS_04_sorted <- file.choose()
idx_TS_04 <- indexBam(TS_04_sorted)
GA_TS_04 <- readGAlignments(TS_04_sorted, param=ScanBamParam(what="flag"))

TS_05_sorted <- file.choose()
idx_TS_05 <- indexBam(TS_05_sorted)
GA_TS_05 <- readGAlignments(TS_05_sorted, param=ScanBamParam(what="flag"))

TS_09_sorted <- file.choose()
idx_TS_09 <- indexBam(TS_09_sorted)
GA_TS_09 <- readGAlignments(TS_09_sorted, param=ScanBamParam(what="flag"))

TS_10_sorted <- file.choose()
idx_TS_10 <- indexBam(TS_10_sorted)
GA_TS_10 <- readGAlignments(TS_10_sorted, param=ScanBamParam(what="flag"))

TS_11_sorted <- file.choose()
idx_TS_11 <- indexBam(TS_11_sorted)
GA_TS_11 <- readGAlignments(TS_11_sorted, param=ScanBamParam(what="flag"))

TS_12_sorted <- file.choose()
idx_TS_12 <- indexBam(TS_12_sorted)
GA_TS_12 <- readGAlignments(TS_12_sorted, param=ScanBamParam(what="flag"))

TS_13_sorted <- file.choose()
idx_TS_13 <- indexBam(TS_13_sorted)
GA_TS_13 <- readGAlignments(TS_13_sorted, param=ScanBamParam(what="flag"))

## Make a list of the genomic alignment objects ##
GAlist <- GAlignmentsList(c(TS_01_sorted, TS_02_sorted, TS_03_sorted, TS_04_sorted, TS_05_sorted, 
            TS_09_sorted, TS_10_sorted, TS_11_sorted, TS_12_sorted, TS_13_sorted))

#estimateLibComplexity(readsDupFreq(TS_01_sorted)) # Used to test if one of the objects would plot as the class of object we want for plotting

## Get reference genome transcripts 
txdb <- makeTxDbFromGFF(file.choose(), dataSource = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_019175285.1", organism = "Sceloporus undulatus")
txs <- transcripts(txdb)

## Make correclation plot(s)
plotCorrelation(GAlignmentsList(GAlist, txs))
