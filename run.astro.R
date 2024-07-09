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

setwd("../../")

# GENE SELECTION
astro.select <- readRDS("astro.select.rds")
# select cv cutoff
smooth.astro <- astro.select[[1]]
astro.cv <- astro.select[[2]]
print(summary(astro.cv))
cutoff <- 0.3
high.cv.astro <- colnames(smooth.astro)[which(astro.cv > cutoff)]
print(length(high.cv.astro))

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

# load spatial data
spatial <- readRDS("sim.spatial.rds")

# load decon
decon <- readRDS("decon.rds")
non_zero_prop <- which(decon > 0)

print(paste0("non_zero_prop: ", length(non_zero_prop)))

lambda = 8; n = 8; m = 20

# make sure "df.nearest.m" is a global variable
# trim the dist matrix to keep the nearest m spots around each spot, except for itself
# get a new matrix: three cols, spot1 index, spot2 index, dist between them
dist <- c(apply(dist_mtx, 1, function(x) x[order(x)[2:(m+1)]])) # the min is the spot itself
spot1 <- rep(1:nrow(dist_mtx), each = m)
spot2 <- c(apply(dist_mtx, 1, function(x) order(x)[2:(m+1)]))
# trimmed dist matrix
df.nearest.m <- data.frame(spot1, spot2, dist, stringsAsFactors=FALSE)



###############

n_spot <- 9852
n_ct <- 5
n_cpu <- 5

# record results
astro.result <- data.frame(gene = high.cv.astro, cv = astro.cv[which(astro.cv > cutoff)],
                           result = rep(0, length(high.cv.astro)))
astro.result.output <- list()
astro.result.loess.ep <- list()


# RUN ON ASTRO
start_time <- Sys.time()
for(gene_of_choice in high.cv.astro){ 
  print(gene_of_choice)
  
  # naming convention
  decon_input <- as.matrix(decon)
  zp <- spatial[,gene_of_choice]
  
  norm.E <- norm_E_simulate(gene_of_choice, n_spot, n_ct, norm.astro, norm.oligo, norm.macro, norm.glu, norm.gaba)
  print(paste0("norm E: ", dim(norm.E)))
  
  # initialize
  E_list <- list()
  for(i in 1:ncol(decon_input)){
    E_list[[i]] <- create_init_E_ct(gene_of_choice, decon_input, i)
  }
  init.E <- do.call(cbind, E_list)
  # # subset non_zero_prop
  init.E <- init.E[non_zero_prop]
  print(paste0("init E: ", length(init.E)))

  # convert back to matrix
  decon_input <- as.matrix(decon)
  print("Before running optim...")
  
  # Run optim, deconvolution
  run <- run_optim(exp_spot_celltype = c(init.E), lambda, n, m, decon_input, zp, ep = norm.E, n_cpu)
  print("After running optim...")
  output_exp <- run[[2]] # entire output exp
  rownames(output_exp) <- rownames(norm.E)
  
  # order:  "Astro"         "GABAergic"     "Glutamatergic" "Macrophage"    "Oligo"  
  eval <- eval_correlation(output_exp, norm.E, norm.astro, 1, pure_ct_index = c())
  
  print(eval[[1]])
  astro.result[gene_of_choice, "result"] <- eval[[1]]
  astro.result.output <- list.append(astro.result.output, output_exp) # prediction
  astro.result.loess.ep <- list.append(astro.result.loess.ep, eval[[2]]) # loess version of gold standard
}

end_time <- Sys.time()
print(end_time - start_time)


saveRDS(astro.result, "astro.cor.rds")
saveRDS(astro.result.output, "astro.output.rds")
saveRDS(astro.result.loess.ep, "astro.loess.ep.rds")

