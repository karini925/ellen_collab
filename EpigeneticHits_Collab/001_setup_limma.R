#Author: Karin Isaev, karin.isaev@gmail.com 

date = Sys.Date()
print(date)

#methods based on 
#https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html#unsupervised-clustering-of-samples

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

dge = DGEList(counts, group = group) #DGEList ready to go!
#take a look at dge
dge 

#[2] - transform raw counts! 

#why? want to account for library size differences. CPM (counts per million) is one such transformation.
#CPM however does not account for gene length differences
#differential expression analyses look at gene expression changes 
#between conditions rather than comparing expression across multiple genes or 
#drawing conclusions on absolute levels of expression

cpm <- cpm(dge$counts)
lcpm <- cpm(dge$counts, log=TRUE)

#get mean and median library sizes 
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6
c(L, M)

#remove lowly expressed genes 
table(rowSums(dge$counts==0)==10) 
keep.exprs <- filterByExpr(dge, group=group)
dge <- dge[keep.exprs, keep.lib.sizes=FALSE]
dim(dge$counts)

#density of log-CPM values for raw pre-filtered data (A) and post-filtered data (B) for each sample 
lcpm.cutoff <- log2(10/M + 2/L)
nsamples <- ncol(dge$counts)
col <- brewer.pal(nsamples, "Paired")
samplenames = colnames(dge$counts)

par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)

for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

lcpm <- cpm(dge, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)

for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

#[3] - normalising gene expression distributions 
x = dge
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#[4] - unsupervised clustering of samples
par(mfrow=c(1,1))

lcpm <- cpm(x, log=TRUE)
group = factor(group)
col.group <- group
levels(col.group) <-  brewer.pal(6, "Set2")
col.group <- as.character(col.group)

plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

#[5] - differential expression analysis 
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
  Asxl2vsScr = Asxl2-Scr, 
  Kdm6avsScr = Kdm6a - Scr, 
  P53vsScr = P53 - Scr, 
  Setd2vsScr = Setd2 - Scr,
  levels = colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <- voom(x, design, plot=TRUE, normalize.method = "quantile")
v

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

#tfit <- treat(vfit, lfc=1)

dt <- decideTests(efit)
summary(dt)

glMDPlot(efit, coef=1, status=dt, main=colnames(efit)[1], 
         side.main="ENTREZID", counts=lcpm, groups=group, launch=TRUE)


glMDPlot(efit, coef=2, status=dt, main=colnames(efit)[2],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=TRUE)


glMDPlot(efit, coef=3, status=dt, main=colnames(efit)[3],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=TRUE)


glMDPlot(efit, coef=4, status=dt, main=colnames(efit)[4],
         side.main="ENTREZID", counts=lcpm, groups=group, launch=TRUE)

####heatmaps

pdf(paste(date, "crispr_exps_differential_expression_KI.pdf", sep="_"), width=12, height=10)

##positive control 

tp53_scr =topTreat(efit, coef=3, n=Inf)
tp53_scr_top_genes <- rownames(tp53_scr)[1:100]
i <- which(rownames(v$E) %in% tp53_scr_top_genes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row",main = "tp53_scr",  
          labRow=rownames(v$E)[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
saveRDS(tp53_scr, file=paste(date, "DE_limma_tp53_scr.rds", sep="_"))


#actual targets 

#Asxl2vsScr
Asxl2vsScr =topTreat(efit, coef=1, n=Inf)
test = Asxl2vsScr
test = subset(test, adj.P.Val < 0.1)
test_genes <- rownames(test)[1:100]
i <- which(rownames(v$E) %in% test_genes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row", main = "Asxl2vsScr", 
          labRow=rownames(v$E)[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
saveRDS(Asxl2vsScr, file=paste(date, "DE_limma_Asxl2vsScr.rds", sep="_"))


#Kdm6avsScr
Kdm6avsScr =topTreat(efit, coef=2, n=Inf)
test = Kdm6avsScr
test = subset(test, adj.P.Val < 0.1)
test_genes <- rownames(test)[1:100]
i <- which(rownames(v$E) %in% test_genes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row", main = "Kdm6avsScr", 
          labRow=rownames(v$E)[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")
saveRDS(Kdm6avsScr, file=paste(date, "DE_limma_Kdm6avsScr.rds", sep="_"))


#Setd2vsScr
Setd2vsScr =topTreat(efit, coef=4, n=Inf)
test = Setd2vsScr
test = subset(test, adj.P.Val < 0.05)
test_genes <- rownames(test)[1:100]
i <- which(rownames(v$E) %in% test_genes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(lcpm[i,], scale="row", main = "Setd2vsScr", 
          labRow=rownames(v$E)[i], labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

dev.off()
saveRDS(Setd2vsScr, file=paste(date, "DE_limma_Setd2vsScr.rds", sep="_"))

#get list of genes for each situation 

tp53 = names(which(dt[,3] != 0))
asxl2 = names(which(dt[,1] != 0))
kdm6 = names(which(dt[,2] != 0))
setd2 = names(which(dt[,4] != 0))

#venn diagram 

pdf(paste(date, "crispr_exps_differential_expression_venn_diagrams_KI.pdf", sep="_"), width=12, height=10)

#all 
vennDiagram(dt[,1:4], circle.col=c("turquoise", "salmon", "orange", "green"))

#tp53 versus asxl2
vennDiagram(dt[,c(3,1)], circle.col=c("turquoise", "salmon", "orange", "green"))

#tp53 versus kdm6
vennDiagram(dt[,3:2], circle.col=c("turquoise", "salmon", "orange", "green"))

#tp53 versus setd2
vennDiagram(dt[,3:4], circle.col=c("turquoise", "salmon", "orange", "green"))

#comparing all three candidates
vennDiagram(dt[,c(1,2,4)], circle.col=c("turquoise", "salmon", "orange", "green"))

dev.off()



