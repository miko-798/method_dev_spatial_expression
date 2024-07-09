
#list.of.packages <- c("ggplot2", "Rcpp")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

# library(ggsignif, lib="/hpc/home/ml509/R_libs")
# install.packages("ggpubr",lib="/hpc/home/ml509/R_libs")

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
library(doParallel)
library(optimParallel)
library(rlist)

suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(plotly))
suppressMessages(library(EnvStats))
suppressMessages(library(Rfast))
suppressMessages(library(data.table))

source("../functions.R")
setwd("/hpc/group/jilab/miko/decon_method/biccn.slideseq/eval")

# order:  "Astro"         "GABAergic"     "Glutamatergic" "Macrophage"    "Oligo" 
ct.order <- 3

# GENE SELECTION
astro.select <- readRDS("gene.select.0.05/glu.select.rds")
# select cv cutoff
smooth.astro <- astro.select[[1]]
astro.cv <- astro.select[[2]]
#print(summary(astro.cv))

i <- as.numeric(commandArgs(trailingOnly = T))
batch_size <- 200
# NOTE: add () for upper and lower limits
astro.cv <- astro.cv[ (1 + (i-1)*batch_size) : (min(i*batch_size, length(astro.cv))) ]
print(i)
print(batch_size)
print(length(astro.cv))

# no cutoff
cutoff <- 0
high.cv.astro <- names(astro.cv)
#high.cv.astro <- colnames(smooth.astro)[which(astro.cv > cutoff)]
print(paste0("number of genes: " , length(high.cv.astro)))


# GLOBAL VARIABLES

# read these in
norm.astro <- readRDS("norm.astro.rds")
norm.oligo <- readRDS("norm.oligo.rds")
norm.macro <- readRDS("norm.macro.rds")
norm.glu <- readRDS("norm.glu.rds")
norm.gaba <- readRDS("norm.gaba.rds")

# read coord
all.coord <- readRDS("all.coord.rds")
#all.coord <- as.matrix(all.coord)
#rownames(all.coord) <- c(1:n_spot)

# calculate dist matrix
dist_mtx <- as.matrix(dist(all.coord))
#dist_mtx <- dist_mtx/10^floor(log10(min(dist_mtx[dist_mtx > 0])))

# load spatial data
spatial <- readRDS("sim.spatial.rds")

# load deconvolution matrix
decon <- readRDS("rctd.prop.rds")

# OPERATIONS ON RCTD CELL TYPE PROPORTIONS:
# get (spot, cell type) positions where the proportions is close to 0, and replace them with 0
# b/c RCTD forces each cell type to have a non-zero prop, while in reality some cell types
# don't exist at all
close_to_zero <- between(decon, 0, 0.01)  # 1% as cutoff for 0 proportions
decon[close_to_zero] = 0
non_zero_prop <- which(decon > 0)

# after changing this, make the total prop as 1:
decon = sweep(decon, 1, rowSums(decon), '/')


lambda = 1; n = 5; m = 20

# make sure "df.nearest.m" is a global variable
# trim the dist matrix to keep the nearest m spots around each spot, except for itself
# get a new matrix: three cols, spot1 index, spot2 index, dist between them
dist <- c(apply(dist_mtx, 1, function(x) x[order(x)[2:(m+1)]])) # the min is the spot itself
spot1 <- rep(1:nrow(dist_mtx), each = m)
spot2 <- c(apply(dist_mtx, 1, function(x) order(x)[2:(m+1)]))
# trimmed dist matrix
df.nearest.m <- data.frame(spot1, spot2, dist, stringsAsFactors=FALSE)

### NORMALIZATION for distance matrix### 
df.nearest.m[,3] <- df.nearest.m[,3]/max(df.nearest.m[,3]) + 1


###############

n_spot <- 9852
n_ct <- 5
n_cpu <- 16

# record results
astro.result <- data.frame(gene = high.cv.astro, cv = astro.cv,
                           result = rep(0, length(high.cv.astro)))
#astro.result <- data.frame(gene = high.cv.astro, cv = astro.cv[which(astro.cv > cutoff)],
#                           result = rep(0, length(high.cv.astro)))
astro.result.output <- list()


# RUN ON ASTRO
start_time <- Sys.time()
for(gene_of_choice in high.cv.astro){ 
  print(gene_of_choice)
  
  # naming convention
  decon_input <- as.matrix(decon)
  zp <- spatial[,gene_of_choice]
  
  norm.E <- norm_E_simulate(gene_of_choice, n_spot, n_ct, norm.astro, norm.oligo, norm.macro, norm.glu, norm.gaba)
  #print(paste0("norm E: ", dim(norm.E)))
  
  # initialize
  #E_list <- list()
  #for(j in 1:ncol(decon_input)){
  #  E_list[[j]] <- create_init_E_ct(gene_of_choice, decon_input, j)
  #}
  #init.E <- do.call(cbind, E_list)
  # # subset non_zero_prop
  #init.E <- init.E[non_zero_prop]
  #print(paste0("init E: ", length(init.E)))

  #initialize
  init.E = rnorm(length(non_zero_prop))

  # convert back to matrix
  decon_input <- as.matrix(decon)
  print("Before running optim...")
  
  # Run optim, deconvolution
  run <- run_optim(exp_spot_celltype = c(init.E), lambda, n, m, decon_input, zp, ep = norm.E, n_cpu)
  print("After running optim...")
  output_exp = run[[2]]

  #output_exp_n5 = run[[2]]
  #output_exp_n1 <- run[[2]] # entire output exp
  rownames(output_exp) <- rownames(norm.E)
  
  # order:  "Astro"         "GABAergic"     "Glutamatergic" "Macrophage"    "Oligo"  
  eval <- eval_correlation_no_na(output_exp, norm.E, ct.order)
  print(eval)

  astro.result[gene_of_choice, "result"] <- eval[ct.order]
  astro.result.output <- list.append(astro.result.output, output_exp) # prediction
}

end_time <- Sys.time()
print(end_time - start_time)

setwd("/hpc/group/jilab/miko/decon_method/biccn.slideseq/eval/run_rctd_prop/batches/glu/results")
saveRDS(astro.result, paste0(i, ".glu.cor.rds"))
saveRDS(astro.result.output, paste0(i, ".glu.output.rds"))




