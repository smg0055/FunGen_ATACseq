## Load required librarues - Download using BiocManager if needed
library(BiocManager)
library(MotifDb)
library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFeatures)
library(ATACseqQC)

## Load all sorted and mapped BAM files
TS_01_sorted <- file.choose()
## Transform BAM files into genomic alignments object(s)
GA_TS_01 <- readGAlignments(TS_01_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_02_sorted <- file.choose()
GA_TS_02 <- readGAlignments(TS_02_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_03_sorted <- file.choose()
GA_TS_03 <- readGAlignments(TS_03_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_04_sorted <- file.choose()
GA_TS_04 <- readGAlignments(TS_04_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_05_sorted <- file.choose()
GA_TS_05 <- readGAlignments(TS_05_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_09_sorted <- file.choose()
GA_TS_09 <- readGAlignments(TS_09_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_10_sorted <- file.choose()
GA_TS_10 <- readGAlignments(TS_10_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_11_sorted <- file.choose()
GA_TS_11 <- readGAlignments(TS_11_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_12_sorted <- file.choose()
GA_TS_12 <- readGAlignments(TS_12_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

TS_13_sorted <- file.choose()
GA_TS_13 <- readGAlignments(TS_13_sorted, param=ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                            isProperPair = TRUE)))

## Create genomic alignment lists
#### Split between control and heat treated samples due to memory constraints 
GAlistC <- GAlignmentsList(c(GA_TS_01, GA_TS_03, GA_TS_10, GA_TS_11, GA_TS_12))
GAlistHT <- GAlignmentsList(c(GA_TS_05, GA_TS_02, GA_TS_09, GA_TS_04, GA_TS_13))

## Import GFF reference file as a txdb object and extract transcripts
txdb <- makeTxDbFromGFF(file.choose(), dataSource = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_019175285.1", organism = "Sceloporus undulatus")
txs <- transcripts(txdb)

## Make correlation plot using both lists - Heatmap can be replaced with PCA if preferred
plotCorrelation(GAlignmentsList(c(GAlistC, GAlistHT), txs, type = "heatmap"))
