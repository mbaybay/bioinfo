## CS640 Bioinformatics - Fall 2017
## Assignment   : DESeq2 InClass
## Submitted By : Melanie Baybay
## Last Modified: Oct. 31, 2017

# Read the data in "count_matrix_sub.txt" 
# and run a DESeq2 experiment comparing the transcriptome 
#   of untreated MCF-7 tumor cells to those treated with estrogen. 
# There are 7 replicates for each of the conditions (untreated, treated) 

### RUN next two lines once 
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library('DESeq2')

setwd("/Users/MelanieBaybay/Documents/USF/2017-2018Fall/Bioinformatics/R-Assignments")
# row = gene, column = subsample at various depths 
#   (i.e. first 3 cols: sample from X2012.562 at depth 2.5M, 10M, 30M)
countTable <- read.table("count_matrix_sub.txt", header=TRUE, sep="\t", row.names=1)

# These substrings of column names indicate the 7 replicates for control sample: 
ControlCistrackID <- c("2012.562","2012.563","2012.564","2012.565","2012.566","2012.568","2012.569")

# These substrings of column names indicate the 7 replicates for sample treated with estrogen:
E2CistrackID <- c("2012.570","2012.571","2012.572","2012.574","2012.575","2012.576","2012.577" )

# These substrings of column names indicate the 3 different depth of reads 
#   for each of the replicates of both untreated and treated:
readDepth <- c("2.5M", "10M", "30M") # bases sequenced per sample (actually down-sampled from 30M)
numRep <- 7   # 7 repetitions each of untreated, treated

# Produce  dataset for read depth ("2.5M")
# The dataset should include 7 replicates of untreated and 7 replicates of treated.
makeDataSet<- function ( data, controlIDs, E2IDs, readDepth){
  controlReplicates <- paste( "X", controlIDs, ".subsamp.", sep="")
  controlReplicates_readDepth <- paste(controlReplicates, readDepth, sep="")
  E2Replicates <- paste( "X", E2IDs, ".subsamp.", sep="")
  E2Replicates_readDepth <- paste(E2Replicates, readDepth, sep="")  
  dataSet <- data[, c(controlReplicates_readDepth, E2Replicates_readDepth)]
  return(dataSet)
}

# GENERATE DATA SETS
dataset_2.5 <- makeDataSet(countTable, ControlCistrackID, E2CistrackID, "2.5M")
dataset_10 <- makeDataSet(countTable, ControlCistrackID, E2CistrackID, "10M")
dataset_30 <- makeDataSet(countTable, ControlCistrackID, E2CistrackID, "30M")

sampleCondition<- c(rep("untreated",numRep), rep("treated",numRep) )

# GENERATE SAMPLE TABLES
sampleTable_2.5<-data.frame(sampleName=colnames(dataset_2.5), condition=sampleCondition)
sampleTable_10 <- data.frame(sampleName=colnames(dataset_10), condition=sampleCondition)
sampleTable_30 <- data.frame(sampleName=colnames(dataset_30), condition=sampleCondition)

# GENERATE DESEQ DATASET from MATRIX
deseq2dataset_2.5 <- DESeqDataSetFromMatrix(countData = dataset_2.5,
                                        colData = sampleTable_2.5,
                                        design = ~ condition)
deseq2dataset_10 <- DESeqDataSetFromMatrix(countData = dataset_10,
                                            colData = sampleTable_10,
                                            design = ~ condition)
deseq2dataset_30 <- DESeqDataSetFromMatrix(countData = dataset_30,
                                           colData = sampleTable_30,
                                           design = ~ condition)


# EXTRACT UP-REGULATED AND DOWN-REGULATED GENES
# --- 2.5 ---
colData(deseq2dataset_2.5)$condition<-factor(colData(deseq2dataset_2.5)$condition, levels=c("untreated","treated"))
dds_2.5<-DESeq(deseq2dataset_2.5)     # *** stats for differential expression done here
res_2.5<-results(dds_2.5)  
#head(res)
res_2.5<-res_2.5[order(res_2.5$padj),]
#head(res)
resSig_2.5 <- subset(res_2.5, res_2.5$padj < 0.05 ) # significant fold changes
resSig_2.5
# -- up-regulated -- 
genesUp_2.5 = subset(resSig_2.5, resSig_2.5$log2FoldChange >=2 )    # up-regulated genes at log2FoldChange >= 2
genesUp_2.5 <-genesUp_2.5[order(genesUp_2.5$log2FoldChange),]
genesUp_2.5_names = rownames(genesUp_2.5)
# -- down-regulated -- 
genesDown_2.5 = subset(resSig_2.5, resSig_2.5$log2FoldChange <= -1.5 )  # down-regulated genes at log2FoldChange <= -1.5
genesDown_2.5 <-genesDown_2.5[order(genesDown_2.5$log2FoldChange),]
genesDown_2.5_names = rownames(genesDown_2.5)

# --- 10 ---
colData(deseq2dataset_10)$condition<-factor(colData(deseq2dataset_10)$condition, levels=c("untreated","treated"))
dds_10<-DESeq(deseq2dataset_10)     # *** stats for differential expression done here
res_10<-results(dds_10)  
#head(res)
res_10<-res_10[order(res_10$padj),]
#head(res)
resSig_10 <- subset(res_10, res_10$padj < 0.05 ) # significant fold changes
resSig_10
# -- up-regulated -- 
genesUp_10 = subset(resSig_10, resSig_10$log2FoldChange >= 2.5 )    # up-regulated genes at log2FoldChange >= 2
genesUp_10 <-genesUp_10[order(genesUp_10$log2FoldChange),]
genesUp_10_names = rownames(genesUp_10)
# -- down-regulated -- 
genesDown_10 = subset(resSig_10, resSig_10$log2FoldChange <= -1.6 )  # down-regulated genes at log2FoldChange <= -1.5
genesDown_10 <-genesDown_10[order(genesDown_10$log2FoldChange),]
genesDown_10_names = rownames(genesDown_10)

# --- 30 ---
colData(deseq2dataset_30)$condition<-factor(colData(deseq2dataset_10)$condition, levels=c("untreated","treated"))
dds_30<-DESeq(deseq2dataset_30)     # *** stats for differential expression done here
res_30<-results(dds_30)  
#head(res)
res_30<-res_30[order(res_30$padj),]
#head(res)
resSig_30 <- subset(res_30, res_30$padj < 0.05 ) # significant fold changes
resSig_30
# -- up-regulated -- 
genesUp_30 = subset(resSig_30, resSig_30$log2FoldChange >= 3 )    # up-regulated genes at log2FoldChange >= 2
genesUp_30 <-genesUp_30[order(genesUp_30$log2FoldChange),]
genesUp_30_names = rownames(genesUp_30)
# -- down-regulated -- 
genesDown_30 = subset(resSig_30, resSig_30$log2FoldChange <= -1.8 )  # down-regulated genes at log2FoldChange <= -1.5
genesDown_30 <-genesDown_30[order(genesDown_30$log2FoldChange),]
genesDown_30_names = rownames(genesDown_30)

# --- comparing UP-REGULATED GENES based on read-depth ---
# ? Which are upregulated in all three read depths
length(genesUp_2.5_names)
length(genesUp_10_names)
length(genesUp_30_names)

upGenes = genesUp_10_names[genesUp_10_names %in% genesUp_30_names]
upGenes = genesUp_2.5_names[genesUp_2.5_names %in% upGenes]
upGenes # genes that are up-regulated in all read-depths

# --- comparing DOWN-REGULATED GENES based on read-depth ---
length(genesDown_2.5_names)
length(genesDown_10_names)
length(genesDown_30_names)

downGenes = genesDown_10_names[genesDown_10_names %in% genesDown_30_names]
downGenes = genesDown_2.5_names[genesDown_2.5_names %in% downGenes]
downGenes # genes that are up-regulated in all read-depths


## ----- VISUALIZATIONS -----
dev.new()
par(mfrow= c(3,2))
#pvalue is the Wald test p-value: condition treated vs untreated
#padj is BH (Benjamini-Hochberg) adjusted p-values <- adjusted to control the false discovery rate
hist(res_2.5$pvalue, breaks=50, col="grey", ylim=c(0,6000))
hist(res_2.5$padj, breaks=50, col="grey", ylim=c(0,6000))
hist(res_10$pvalue, breaks=50, col="red", ylim=c(0,6000))
hist(res_10$padj, breaks=50, col="red", ylim=c(0,6000))
hist(res_30$pvalue, breaks=50, col="blue", ylim=c(0,6000))
hist(res_30$padj, breaks=50, col="blue", ylim=c(0,6000))



# MA plot is application of a Bland-Altman plot for visual representation of two channel gene expression data
#  which has been transformed onto the M (log ratios) and A (mean average) scale.
# M vs. A plots of each pair (untreated, treated) is produced.
dev.new()
par(mfrow= c(1,3))
plotMA(dds_2.5,ylim=c(-3,3),main= paste("DESeq2 ", "2.5M"), alpha=.05)
plotMA(dds_10,ylim=c(-3,3),main= paste("DESeq2 ", "10M"), alpha=.05)
plotMA(dds_30,ylim=c(-3,3),main= paste("DESeq2 ", "30M"), alpha=.05)
# dev.off()


#transform the raw counts for clustering:
#choose blind so that conditions does not influence the outcome, 
# to see if conditions cluster based purely on the individual datasets, in an unbiased way.
rld_2.5 <- rlogTransformation(dds_2.5, blind=TRUE) # regularized log
vsd_2.5 <- varianceStabilizingTransformation(dds_2.5, blind=TRUE) # DESeq's variance stabilisation
rld_10 <- rlogTransformation(dds_10, blind=TRUE) 
vsd_10 <- varianceStabilizingTransformation(dds_10, blind=TRUE) 
rld_30 <- rlogTransformation(dds_30, blind=TRUE) 
vsd_30 <- varianceStabilizingTransformation(dds_30, blind=TRUE) 

################## Biotech may stop here ###############################

#see which approaches have more consistent SD across the read counts
# source("http://bioconductor.org/biocLite.R")
# biocLite("vsn")
library("vsn")

notAllZero_2.5 <- (rowSums(counts(dds_2.5))>0)
par(mfrow=c(3,1))
# log2 doesn't do well for low read counts (SD is high, varies)
meanSdPlot(log2(counts(dds_2.5,normalized=TRUE)[notAllZero_2.5,] + 1))
# regularized log (center) and DESeq's variance stabilisation (right)
#  transformations do better across the  range of counts
meanSdPlot(assay(rld_2.5[notAllZero_2.5,]))
meanSdPlot(assay(vsd_2.5[notAllZero_2.5,]))

notAllZero_10 <- (rowSums(counts(dds_10))>0)
par(mfrow=c(3,1))
# log2 doesn't do well for low read counts (SD is high, varies)
meanSdPlot(log2(counts(dds_10,normalized=TRUE)[notAllZero_10,] + 1))
# regularized log (center) and DESeq's variance stabilisation (right)
#  transformations do better across the  range of counts
meanSdPlot(assay(rld_10[notAllZero_10,]))
meanSdPlot(assay(vsd_10[notAllZero_10,]))

notAllZero_30 <- (rowSums(counts(dds_30))>0)
par(mfrow=c(3,1))
# log2 doesn't do well for low read counts (SD is high, varies)
meanSdPlot(log2(counts(dds_30,normalized=TRUE)[notAllZero_30,] + 1))
# regularized log (center) and DESeq's variance stabilisation (right)
#  transformations do better across the  range of counts
meanSdPlot(assay(rld_30[notAllZero_30,]))
meanSdPlot(assay(vsd_30[notAllZero_30,]))


par(mfrow= c(1,1))

# calculate sample to sample Euclidean distances betw lg fold change of genes < WRT what? ********
#   if dist is near zero, then very similar (dark blue is same sample zero dist)
# use rlog-transformed data to avoid domination by a few highly variable genes 
# dist calculates distances between rlog-transforms in data rows
#  our samples constitute columns so use t to transpose the matrix
distsRL_2.5 <- dist(t(assay(rld_2.5)))   
mat_2.5 <- as.matrix(distsRL_2.5)
rownames(mat_2.5) <- colnames(mat_2.5) <- with(colData(dds_2.5),
                                       paste(condition,sampleName, sep=" : "))

distsRL_10 <- dist(t(assay(rld_10)))   
mat_10 <- as.matrix(distsRL_10)
rownames(mat_10) <- colnames(mat_10) <- with(colData(dds_10),
                                               paste(condition,sampleName, sep=" : "))

distsRL_30 <- dist(t(assay(rld_30)))   
mat_30 <- as.matrix(distsRL_30)
rownames(mat_30) <- colnames(mat_30) <- with(colData(dds_30),
                                             paste(condition,sampleName, sep=" : "))

# visualizes distances and clusters:
#From the Apr 2015 vignette
# install.packages("gplots")
library("RColorBrewer")
library("gplots")
select_2.5 <- order(rowMeans(counts(dds_2.5,normalized=TRUE)),decreasing=TRUE)[1:30]
select_10 <- order(rowMeans(counts(dds_10,normalized=TRUE)),decreasing=TRUE)[1:30]
select_30 <- order(rowMeans(counts(dds_30,normalized=TRUE)),decreasing=TRUE)[1:30]

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

dev.new()
hc_2.5 <- hclust(distsRL_2.5)
heatmap.2(mat_2.5, Rowv=as.dendrogram(hc_2.5),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(18, 18))
dev.copy(png,"deseq2_2.5_heatmaps_samplebysample.png")
dev.off()

dev.new()
hc_10 <- hclust(distsRL_10)
heatmap.2(mat_10, Rowv=as.dendrogram(hc_10),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(18, 18))
dev.copy(png,"deseq2_10_heatmaps_samplebysample.png")
dev.off()

dev.new()
hc_30 <- hclust(distsRL_30)
heatmap.2(mat_30, Rowv=as.dendrogram(hc_30),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(18, 18))
dev.copy(png,"deseq2_30_heatmaps_samplebysample.png")
dev.off()


#PCA
dev.new()
colours <- c(rgb(1:3/4,0,0),rgb(0,1:3/4,0),rgb(0,0,1:3/4),rgb(1:3/4,0,1:3/4))
plotPCA( rld_2.5, intgroup = c("sampleName","condition") )
dev.copy(png,"deseq2_2.5_PCA.png")
dev.off()

dev.new()
colours <- c(rgb(1:3/4,0,0),rgb(0,1:3/4,0),rgb(0,0,1:3/4),rgb(1:3/4,0,1:3/4))
plotPCA( rld_10, intgroup = c("sampleName","condition") )
dev.copy(png,"deseq2_10_PCA.png")
dev.off()


dev.new()
colours <- c(rgb(1:3/4,0,0),rgb(0,1:3/4,0),rgb(0,0,1:3/4),rgb(1:3/4,0,1:3/4))
plotPCA( rld_30, intgroup = c("sampleName","condition") )
dev.copy(png,"deseq2_30_PCA.png")
dev.off()

