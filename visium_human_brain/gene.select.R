library(SPOTlight)
library(Seurat)
library(foreach)
library(doRNG)  ## For reproducibility
library(doParallel)
library(igraph)
library(RColorBrewer)
library(ggrepel)
library(anndata)
library(parallel)
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(EnvStats))
suppressMessages(library(data.table))

# coord
setwd("data/")
all.coord = readRDS("pixel.coord.rds") 
	
# read these in
load("ct.gold.stand.RData")

# gene selection for a specific cell type
# INPUT: the cell type expression gold standard
gene_select <- function(norm.ct){
  ct.coord <- all.coord[rownames(norm.ct),] # subset spots

  # filter genes: gene expression > 0 in > 50% spots
  # use mclapply: for each gene, find % spots with > 0 expression
  ct.positive <-  mclapply(colnames(norm.ct), function(x){
    sum(sign(norm.ct[,x]) > 0)/ nrow(norm.ct)
  }, mc.cores=5)

  filter.gene <- colnames(norm.ct)[which(ct.positive > 0.5)]
  print(length(filter.gene))
  filter.E <- norm.ct[,filter.gene]
  #dim(filter.E)
  filter.E <- as.data.frame(filter.E)

  # record cv
  ct.cv <- apply(filter.E, 2, cv)
  print(ct.cv[1:2])
  return(ct.cv)
}


setwd("../")
L23.select <- gene_select(L23.matrix)
saveRDS(L23.select, "genes/L23.genes.rds")

L5.select <- gene_select(L5.matrix)
saveRDS(L5.select, "genes/L5.genes.rds")



