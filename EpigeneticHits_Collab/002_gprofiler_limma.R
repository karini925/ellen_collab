#Author: Karin Isaev, karin.isaev@gmail.com 

date = Sys.Date()
print(date)

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
library(gProfileR)

##Data------------------------------------------------------------------------------------------

#tp53 versus scramble 
tp53_scr = readRDS(list.files(pattern="tp53_scr"))
tp53_scr$gene = rownames(tp53_scr)

all_genes = tp53_scr$gene #use for background cause this is all we analzyed 

tp53_scr = as.data.table(tp53_scr)
tp53_scr = as.data.table(filter(tp53_scr, adj.P.Val < 0.05))

#Asxl2 versus scramble
Asxl2vsScr = readRDS(list.files(pattern="Asxl2vsScr"))
Asxl2vsScr$gene = rownames(Asxl2vsScr)
Asxl2vsScr = as.data.table(Asxl2vsScr)
Asxl2vsScr = as.data.table(filter(Asxl2vsScr, adj.P.Val < 0.05))

#Kdm6a versus scramble
Kdm6avsScr = readRDS(list.files(pattern="Kdm6avsScr"))
Kdm6avsScr$gene = rownames(Kdm6avsScr)
Kdm6avsScr = as.data.table(Kdm6avsScr)
Kdm6avsScr = as.data.table(filter(Kdm6avsScr, adj.P.Val < 0.05))

#Setd2 versus scramble
Setd2vsScr = readRDS(list.files(pattern="Setd2vsScr"))
Setd2vsScr$gene = rownames(Setd2vsScr)
Setd2vsScr = as.data.table(Setd2vsScr)
Setd2vsScr = as.data.table(filter(Setd2vsScr, adj.P.Val < 0.05))

##Analysis--------------------------------------------------------------------------------------

Kdm6avsScr_pathways = gprofiler(Kdm6avsScr$gene, organism = "mmusculus", sort_by_structure = T,
          ordered_query = T, significant = T, exclude_iea = T, underrep = F,
          evcodes = F, region_query = F, max_p_value = 0.05, min_set_size = 5,
          max_set_size = 500, min_isect_size = 2, correction_method = "fdr",
          hier_filtering = "none", domain_size = "annotated", custom_bg = all_genes,
          include_graph = T,png_fn=T)

Asxl2vsScr_pathways = gprofiler(Asxl2vsScr$gene, organism = "mmusculus", sort_by_structure = T,
                                ordered_query = T, significant = T, exclude_iea = T, underrep = F,
                                evcodes = F, region_query = F, max_p_value = 0.05, min_set_size = 5,
                                max_set_size = 500, min_isect_size = 2, correction_method = "fdr",
                                hier_filtering = "none", domain_size = "annotated", custom_bg = all_genes,
                                include_graph = T,png_fn=T)

Setd2vsScr_pathways = gprofiler(Setd2vsScr$gene, organism = "mmusculus", sort_by_structure = T,
                                ordered_query = T, significant = T, exclude_iea = T, underrep = F,
                                evcodes = F, region_query = F, max_p_value = 0.05, min_set_size = 5,
                                max_set_size = 500, min_isect_size = 2, correction_method = "fdr",
                                hier_filtering = "none", domain_size = "annotated", custom_bg = all_genes,
                                include_graph = T,png_fn=T)

tp53_scr_pathways = gprofiler(tp53_scr$gene, organism = "mmusculus", sort_by_structure = T,
                                ordered_query = T, significant = T, exclude_iea = T, underrep = F,
                                evcodes = F, region_query = F, max_p_value = 0.05, min_set_size = 5,
                                max_set_size = 500, min_isect_size = 2, correction_method = "fdr",
                                hier_filtering = "none", domain_size = "annotated", custom_bg = all_genes,
                                include_graph = T,png_fn=T)
