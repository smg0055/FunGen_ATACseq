# Set working directory to the saved files
setwd("C:/Users/fmx398/Functional_genomics/NewData")

# Load required packages
library(Rsubread)
library(grid)
library(DESeq2)
library(data.table)
library(grid)
library(ggplot2)
library(tidyr)
library(edgeR)
library(pheatmap)
library(dplyr)

# Read in the peak files and standardize to smallest file size

# Control samples
TS_01 <- "TS_01.sort.bam_peaks.narrowPeak"
TS_03 <- "TS_03.sort.bam_peaks.narrowPeak"
TS_10 <- "TS_10.sort.bam_peaks.narrowPeak"
TS_11 <- "TS_11.sort.bam_peaks.narrowPeak"
TS_12 <- "TS_12.sort.bam_peaks.narrowPeak"
TS_01 <- fread(TS_01, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_01 <- head(TS_01, 14000)
TS_03 <- fread(TS_03, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_03 <- head(TS_03, 14000)
TS_10 <- fread(TS_10, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_10 <- head(TS_10, 14000)
TS_11 <- fread(TS_11, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_11 <- head(TS_11, 14000)
TS_12 <- fread(TS_12, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_12 <- head(TS_12, 14000)

# Heat Stress samples
TS_02 <- "TS_02.sort.bam_peaks.narrowPeak"
TS_04 <- "TS_04.sort.bam_peaks.narrowPeak"
TS_05 <- "TS_05.sort.bam_peaks.narrowPeak"
TS_09 <- "TS_09.sort.bam_peaks.narrowPeak"
TS_13 <- "TS_13.sort.bam_peaks.narrowPeak"
TS_02 <- fread(TS_02, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_02 <- head(TS_02, 14000)
TS_04 <- fread(TS_04, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_04 <- head(TS_04, 14000)
TS_05 <- fread(TS_05, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_05 <- head(TS_05, 14000)
TS_09 <- fread(TS_09, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_09 <- head(TS_09, 14000)
TS_13 <- fread(TS_13, header = FALSE, col.names = c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"))
TS_13 <- head(TS_13, 14000)

# Seperate by treatment groups
ctrl_peaks <- rbind(TS_01, TS_03, TS_10, TS_11, TS_12)
heat_peaks <- rbind(TS_02, TS_04, TS_05, TS_09, TS_13)

# Add treatment column to control peaks
ctrl_peaks$treatment <- "ctrl"

# Add treatment column to heat peaks
heat_peaks$treatment <- "heat"

# Combine control and heat peaks into one data set
all_peaks <- rbind(ctrl_peaks, heat_peaks)

# combine the sample names into a vector
sample_names <- c("TS_01", "TS_03", "TS_10", "TS_11", "TS_12", "TS_02", "TS_04", "TS_05", "TS_09", "TS_13")

# Extract the signal values from the dataframes for each sample
signal_values <- sapply(list(TS_01, TS_03, TS_10, TS_11, TS_12, TS_02, TS_04, TS_05, TS_09, TS_13), function(df) df$signalValue)

# Convert signal_values to a data frame
signal_values <- as.data.frame(signal_values)

# Set column names to sample_names
colnames(signal_values) <- sample_names

# create the heatmap using pheatmap based on signal value
pheatmap(signal_values, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Differential peak signal across conditions",
         show_rownames = FALSE,
         legend = TRUE)

# Make a box around control sample
grid.rect(x = unit(0.24, "npc"), y = unit(0.48, "npc"),
          width = unit(0.455, "npc"), height = unit(0.92, "npc"),
          gp=gpar(fill = "transparent",col="black", alpha = 0.5,lwd=3))

# Make a box around heat stress samples
grid.rect(x = unit(0.695, "npc"), y = unit(0.48, "npc"),
          width = unit(0.45, "npc"), height = unit(0.92, "npc"),
          gp=gpar(fill = "transparent",col="red", alpha = 0.5,lwd=3))

# Create a data frame with only the signal values as x and y
df <- data.frame(x = ctrl_peaks$signalValue, y = heat_peaks$signalValue)

# Create a scatter plot
ggplot(df, aes(x = x, y = y, color = factor(case_when(x > y ~ "x", x < y ~ "y", TRUE ~ "equal")))) +
  geom_point() +
  xlab("Control peaks") +
  ylab("Heat peaks") +
  xlim(0, 100) +
  ggtitle("Comparison of control and heat peaks") +
  scale_color_manual(values = c("blue", "red")) +
  guides(color = FALSE)

# Loading in bam files
bam_files <- list.files(pattern = "sort.bam$", full.names = TRUE)

# Creating a count matrix

## Uncomment and only run this line once (takes about 15 mins)
count_matrix <- featureCounts(bam_files, annot.ext = "C:/Users/fmx398/Functional_genomics/NewData/genomic.gtf",
                              isPairedEnd = TRUE, isGTFAnnotationFile = TRUE,
                              requireBothEndsMapped = TRUE)

# Extract the counts from the matrix and remove low-abundance peaks
counts <- count_matrix$counts
keep <- rowSums(counts) >= 10
counts <- counts[keep, ]

# Group by treatment and add average peak count by treatment

# Control counts
ctrl_counts <- counts[, c(1, 3, 7, 8, 9)]
ctrl_counts <- cbind(ctrl_counts, rowMeans(ctrl_counts))
colnames(ctrl_counts)[6] <- "rowMeans"

# Heat stress counts
heat_counts <- counts[,c(2,4,5,6,10)]
heat_counts <- cbind(heat_counts, rowMeans(heat_counts))
colnames(heat_counts)[6] <- "rowMeans"

# Combine matrices and select only rowMeans column
total_rowMeans <- as.matrix(heat_counts[, "rowMeans"])
colnames(total_rowMeans)[1] <- "Heat stress"
total_rowMeans <- cbind(total_rowMeans, ctrl_counts[, "rowMeans"])
colnames(total_rowMeans)[2] <- "Control stress"

# Normalizing peaks to heat stress to look for up-regulation
heat_diff <- total_rowMeans[,1] / total_rowMeans[,2]
heat_total_rowMeans <- total_rowMeans
heat_total_rowMeans <- cbind(heat_total_rowMeans, heat_diff)
heat_total_rowMeans <- heat_total_rowMeans[order(heat_diff, decreasing = TRUE), ]

# Subsetting top 20 peaks
heat_total_rowMeans <- heat_total_rowMeans[-c(1:4), ]
heat_sub_top_20 <- heat_total_rowMeans[c(1:20),]

# Bar plot of all peaks showing which are up regulated for heat stress
barplot(heat_total_rowMeans[,3], main="Total Peaks for Heat Stress", xlab="Genes", ylab="Number of Peaks", las=2, names.arg = rep("", nrow(heat_total_rowMeans)))

# Set the margins of the plot to move it up
par(mar = c(5, 4, 4, 4) + 1)

# Bar plot for top 20 peaks for heat stress samples
barplot(heat_sub_top_20[,3], main="Top 20 Peaks for Heat Stress", ylab="Number of Peaks", las=2, cex.names=0.75)

# Normalizing peaks to control group to look for up-regulation
ctrl_diff <- total_rowMeans[,2] / total_rowMeans[,1]
ctrl_total_rowMeans <- total_rowMeans
ctrl_total_rowMeans <- cbind(ctrl_total_rowMeans, ctrl_diff)

# Sub-setting top 20 peaks
ctrl_total_rowMeans <- ctrl_total_rowMeans[order(ctrl_diff, decreasing = TRUE), ]
ctrl_sub_top_20 <- ctrl_total_rowMeans[c(1:20),]

# Bar plot of all peaks for control samples
barplot(ctrl_total_rowMeans[,3], main="Total Peaks for Control Groups", xlab="Genes", ylab="Number of Peaks", las=2, 
        names.arg = rep("", nrow(ctrl_total_rowMeans)), cex.axis = 0.8)

# Set the margins of the plot to move it up
par(mar = c(5, 4, 4, 4) + 1)

# Bar plot for top 20 peaks for control samples
barplot(ctrl_sub_top_20[,3], main="Top 20 Peaks for Control Group", ylab="Number of Peaks", las=2, cex.names=0.75)

# Create a design matrix
treatment <- factor(c(rep("Control", 5), rep("Heat", 5)))
design <- model.matrix(~treatment)

# Create a DGEList object and perform filtering
dge <- DGEList(counts = counts, group = treatment)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
res <- topTags(qlf, n = Inf)

# Extract the significant peaks
sig_peaks <- topTags(qlf, n = Inf)$table
sig_peaks <- sig_peaks[sig_peaks$PValue < 0.05, ]

# Extract gene names from sig_peaks data frame
sig_genes <- as.data.frame(rownames(sig_peaks))

# Get the top 20 peaks for gene anaylsis based on pValue
top_peaks <- head(sig_peaks[order(-sig_peaks$PValue),], 20)

# Create MA plot to show gene regulation
plotSmear(qlf, de.tags = which(qlf$table$FDR < 0.2), main = "Significant Differentially Expressed Peaks (FDR < 0.2)", col = ifelse(qlf$table$PValue < 0.05, "blue", "black"))
abline(h = 0, col = "red")

# Create a volcano plot to show sig_peaks and FDR
ggplot(sig_peaks, aes(x=logFC, y=-log10(PValue), color=FDR)) + 
  geom_point() + 
  scale_color_gradient(low="red", high="blue") +
  labs(x="Log2 Fold Change", y="-Log10 P-value", color="FDR") +
  theme_bw() +
  ggtitle("Significant Differentially Expressed Peaks")

# Create new volcano plot looking at top 20 peaks ad genes associated witht hem
ggplot(top_peaks, aes(x=logFC, y=-log10(PValue), color=FDR)) + 
  geom_point() + 
  geom_text(aes(label=row.names(top_peaks)), hjust=.6, vjust=1.5, size=2) +
  scale_color_gradient(low="blue", high="red") +
  labs(x="Log2 Fold Change", y="-Log10 P-value", color="FDR") +
  theme_bw() +
  ggtitle("Top 20 Significant Differentially Expressed Peaks")
