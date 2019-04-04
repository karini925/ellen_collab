#Author: Karin Isaev, karin.isaev@gmail.com 

date = Sys.Date()
print(date)

#methods based on 
#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#plot-counts

##load in packages-----------------------------------------------------------------------------

library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)
library(data.table)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(gplots)
library(DESeq2)
library(apeglm)

##Data------------------------------------------------------------------------------------------

counts = fread("Ellen_F_rawCounts.txt", data.table=F)
head(counts)

#[1] - make DGEList object for limma voom analysis! 

#make rows gene names and columns sample types 
#counts file will be converted to DGEList object 

rownames(counts) = counts$Geneid
counts = counts[,7:ncol(counts)]

#also need sample grouping 
colnames(counts) = c("Kdm6a_1", "Mll3_1", "Asxl2_1", "P53_1", "Setd2_1", "Kdm6a_2", "Setd2_2", "Scr_1", "Asxl2_2", "Scr_2", "APC_1", "rtStar_1", "APC_2", "P53_2")
z = order(colnames(counts))
counts = counts[,z]
z = which(colnames(counts) %in% c("Mll3_1", "rtStar_1", "APC_1", "APC_2"))
counts = counts[,-z]

group = sapply(colnames(counts), function(x){unlist(strsplit(x, "_"))[1]})
group = as.data.frame(group)

#make deseq object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = group,
                              design = ~ group)

#basic pre-filtering, keep rows with at least 10 counts 
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]

#select reference, in our case we are comparing everything against the scramble
dds$group <- relevel(dds$group, ref = "Scr")

#differential expression analysis 
dds <- DESeq(dds)
res <- results(dds, alpha=0.05)
res

#log fold change shrinkage for visualization and ranking 
Asxl2_lfc = lfcShrink(dds, coef="group_Asxl2_vs_Scr", type="apeglm")
Kdm6a_lfc = lfcShrink(dds, coef="group_Kdm6a_vs_Scr", type="apeglm")
Setd2_lfc = lfcShrink(dds, coef="group_Setd2_vs_Scr", type="apeglm")
P53_lfc = lfcShrink(dds, coef="group_P53_vs_Scr", type="apeglm")

#get each results 
Asxl2 = results(dds, contrast = c("group", "Asxl2", "Scr"))
Kdm6a = results(dds, contrast = c("group", "Kdm6a", "Scr"))
Setd2 = results(dds, contrast = c("group", "Setd2", "Scr"))
P53 = results(dds, contrast = c("group", "P53", "Scr"))
  
summary(Asxl2)
summary(Kdm6a)
summary(Setd2)
summary(P53)

#exploring and exporting results
plotMA(Asxl2_lfc, ylim=c(-2,2))
plotMA(Kdm6a_lfc, ylim=c(-2,2))
plotMA(Setd2_lfc, ylim=c(-2,2))
plotMA(P53_lfc, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="group")
plotCounts(dds, gene=which.max(abs(res$log2FoldChange)), intgroup="group")

#prepare to output and save (adj < 0.05 and abs(log2FC) greater than 1)
dim(Asxl2)
Asxl2 = subset(Asxl2, padj < 0.05)
dim(Asxl2)
z = which((Asxl2$log2FoldChange >=1) | (Asxl2$log2FoldChange < -1))
Asxl2 = Asxl2[z,]
Asxl2_genes = rownames(Asxl2)

#prepare to output and save (adj < 0.05 and abs(log2FC) greater than 1)
dim(Kdm6a)
Kdm6a = subset(Kdm6a, padj < 0.05)
dim(Kdm6a)
z = which((Kdm6a$log2FoldChange >=1) | (Kdm6a$log2FoldChange < -1))
Kdm6a = Kdm6a[z,]
Kdm6a_genes = rownames(Kdm6a)

#prepare to output and save (adj < 0.05 and abs(log2FC) greater than 1)
dim(Setd2)
Setd2 = subset(Setd2, padj < 0.05)
dim(Setd2)
z = which((Setd2$log2FoldChange >=1) | (Setd2$log2FoldChange < -1))
Setd2 = Setd2[z,]
Setd2_genes = rownames(Setd2)

#prepare to output and save (adj < 0.05 and abs(log2FC) greater than 1)
dim(P53)
P53 = subset(P53, padj < 0.05)
dim(P53)
z = which((P53$log2FoldChange >=1) | (P53$log2FoldChange < -1))
P53 = P53[z,]
P53_genes = rownames(P53)

library("pheatmap")

vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("group")])
rownames(df) <- colnames(assay(ntd))
colnames(df) = "group"

pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

plotPCA(vsd, intgroup=c("group"))



