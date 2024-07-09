# script to find correlation between simulated and real spatial data

suppressMessages(library(dplyr))
suppressMessages(library(anndata))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(parallel))


# functions
normalize <- function(mat){ return( log2(mat/rowSums(mat)* 10^5+1))}

# code from arrayMagic package
# correlation between cols
colCors = function(x, y) { 
  sqr = function(x) x*x
  if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
    stop("Please supply two matrices of equal size.")
  x   = sweep(x, 2, colMeans(x))
  y   = sweep(y, 2, colMeans(y))
  cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
  return(cor)
}


################
# LOAD real spatial data: slide-seq
spatial <- read_h5ad("../slideseq_MOp_1217.h5ad") # 9852 × 24518
print(spatial)
real.spatial <- as.data.frame(spatial$X)
# coord
all.coord <- data.frame(x_coor = spatial$obs$x, y_coor = spatial$obs$y)

# normalize
norm.real.spatial <- normalize(real.spatial)
colnames(norm.real.spatial) <- spatial$var_names # gene names

# LOAD simulated spatial data
norm.sim.spatial <- readRDS("sim.spatial.rds")

# find common genes
common.genes <- intersect(colnames(norm.sim.spatial), colnames(norm.real.spatial) )
print(length(common.genes))
# subset the matrices
norm.sim.spatial <- norm.sim.spatial[,common.genes]
norm.real.spatial <- norm.real.spatial[,common.genes]
norm.sim.spatial <- as.data.frame(norm.sim.spatial)
norm.real.spatial <- as.data.frame(norm.real.spatial)

print(dim(norm.sim.spatial)) # dim should be the same
print(dim(norm.real.spatial))


# smoothing for each gene
norm.real.spatial <- mclapply(colnames(norm.real.spatial), function(x){
  fitted(loess(norm.real.spatial[, x] ~ all.coord$x_coor + all.coord$y_coor, norm.real.spatial))}, mc.cores=5)
# smoothed expression
norm.real.spatial <- matrix(unlist(norm.real.spatial), ncol = length(norm.real.spatial))
saveRDS(norm.real.spatial, "norm.real.spatial.rds")

# smoothing
norm.sim.spatial <- mclapply(colnames(norm.sim.spatial), function(x){
  fitted(loess(norm.sim.spatial[, x] ~ all.coord$x_coor + all.coord$y_coor, norm.sim.spatial))}, mc.cores=5)
# smoothed expression
norm.sim.spatial <- matrix(unlist(norm.sim.spatial), ncol = length(norm.sim.spatial))
saveRDS(norm.sim.spatial, "norm.sim.spatial.rds")


# find cor
exp.cor <- colCors(norm.real.spatial, norm.sim.spatial)
saveRDS(exp.cor, "exp.cor.rds")







