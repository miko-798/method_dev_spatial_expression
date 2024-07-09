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
#suppressMessages(library(Rfast))
suppressMessages(library(data.table))
suppressMessages(library(gridExtra))

source("functions.R")

#args = commandArgs(trailingOnly = T)
#gene_of_choice = args[1]
#n = as.numeric(args[2])
#print(gene_of_choice)
#print(n)

## SELECT GENES TO RUN
#glu.select <- readRDS("gene.select.0.05/glu.select.rds")
# select cv cutoff
#glu.cv <- glu.select[[2]]
#glu.cv <- glu.cv[1:10]
#gene_list = names(glu.cv)
#gene_list = c("Tubb5")
#gene_of_choice = "Tubb5"

# GLOBAL VARIABLES
setwd("data/")

# read coord
all.coord <- readRDS("pixel.coord.rds")
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

n_spot <- nrow(spatial)
n_ct <- ncol(decon)
n_cpu <- 16

#grid_n = c(0.1, 0.5, 1, 5, 10, 50,100,200, 500,600,700,800, 1000, 1500,2000)

gene_list = c("FAM127A", "SLC9A6")
grid_n = c(1, 5, 10, 50,100,200)

load("ct.gold.stand.RData")
ct_list = c("L2.3.IT", "L5.IT",  "Lamp5",  "Pvalb",  "Sst",  "Vip")
gold_stand_list = list(L23.matrix, L5.matrix, Lamp.matrix, Pvalb.matrix, Sst.matrix, Vip.matrix)


setwd("../results")
# RUN METHOD
start_time <- Sys.time()

for(gene in gene_list){
  gene_of_choice = gene
  print(gene_of_choice)
  # record results
  cor.result <- data.frame(gene = gene_of_choice, n = grid_n,
                           result = matrix(rep(0, length(grid_n)*n_ct), ncol = n_ct))
  predict.result <- list()

  # naming convention
  decon_input <- as.matrix(decon)
  zp <- spatial[,gene_of_choice]

  norm.E <- norm_E_simulate(gene_of_choice, spatial, n_spot, n_ct, gold_stand_list, ct_list)
  #print(paste0("norm E: ", dim(norm.E)))
  
  #initialize
  #init.E = rnorm(length(non_zero_prop))
  overall.exp = spatial[, gene_of_choice]
  init.E = matrix(rep(overall.exp, n_ct), ncol = n_ct)
  init.E = init.E[non_zero_prop]
  # convert back to matrix
  decon_input <- as.matrix(decon)
 
  for(n in grid_n){ 
    print("Before running optim...")
    # Run optim, deconvolution
    run <- run_optim(exp_spot_celltype = c(init.E), lambda, n, m, decon_input, zp, ep = norm.E, n_cpu)
    print("After running optim...")
    output_exp = run[[2]]
    rownames(output_exp) <- rownames(norm.E)
  
    for(i in c(1:n_ct)){
      eval <- eval_correlation_no_na(output_exp, norm.E, ct.order = i)
      print(eval)
      # update results
      cor.result[cor.result$n == n, 2+i]= eval
    }	

    cor.result[cor.result$n == n,]$result = eval
    predict.result <- list.append(predict.result, output_exp) # prediction
  }
  saveRDS(cor.result, paste0(gene_of_choice, ".cor.result.rds"))
  saveRDS(predict.result, paste0(gene_of_choice, ".predict.result.rds"))

}

end_time <- Sys.time()
print(end_time - start_time)



# plotting
L23.coord <- all.coord[rownames(L23.matrix),] 
L5.coord <- all.coord[rownames(L5.matrix),]

FAM.pred = readRDS("FAM127A.predict.result.rds")
SLC.pred = readRDS("SLC9A6.predict.result.rds")

prediction.L23 = as.data.frame(FAM.pred[[6]][rownames(L23.matrix),])
colnames(prediction.L23) = colnames(decon)

prediction.L5 = as.data.frame(FAM.pred[[6]][rownames(L5.matrix),])
colnames(prediction.L5) = colnames(decon)

gold.L23 = as.data.frame(L23.matrix)
gold.L5 = as.data.frame(L5.matrix)

# gold standard
# glu expression of a gene
L23.exp = ggplot(gold.L23, aes(L23.coord$array_row, L23.coord$array_col)) + geom_point(data = gold.L23, aes(color = FAM127A), size = 1 , alpha = 10) + theme_classic() +
	scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0, 8)) + ggtitle("gold standard FAM127A expression in L23")

L5.exp = ggplot(gold.L5, aes(L5.coord$array_row, L5.coord$array_col)) + geom_point(data = gold.L5, aes(color = FAM127A), size = 1 , alpha = 10) + theme_classic() +
        scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0,8)) + ggtitle("gold standard FAM127A expression in L5")

grid.arrange(L23.exp, L5.exp, nrow = 1)


# GENE: FAM127A
###############
# prediction for L23
pred = ggplot(prediction.L23, aes(L23.coord$array_row, L23.coord$array_col)) + geom_point(aes(color = L2.3.IT), size = 1 , alpha = 10) + theme_classic() +
  scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0,8)) + ggtitle("predicted FAM127A expression in L23")
grid.arrange(L23.exp, pred, nrow = 1)


# prediction for L5
pred = ggplot(prediction.L5, aes(L5.coord$array_row, L5.coord$array_col)) + geom_point(aes(color = L5.IT), size = 1 , alpha = 10) + theme_classic() +
  scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0,8)) + ggtitle("predicted FAM127A expression in L5")
grid.arrange(L5.exp, pred, nrow = 1)

# GENE: SLC9A6
###########

# gold standard

L23.exp = ggplot(gold.L23, aes(L23.coord$array_row, L23.coord$array_col)) + geom_point(data = gold.L23, aes(color = SLC9A6), size = 1 , alpha = 10) + theme_classic() +
        scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0, 8)) + ggtitle("gold standard SLC9A6 expression in L23")

L5.exp = ggplot(gold.L5, aes(L5.coord$array_row, L5.coord$array_col)) + geom_point(data = gold.L5, aes(color = SLC9A6), size = 1 , alpha = 10) + theme_classic() +
        scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0,8)) + ggtitle("gold standard SLC9A6 expression in L5")

grid.arrange(L23.exp, L5.exp, nrow = 1)



prediction.L23 = as.data.frame(SLC.pred[[6]][rownames(L23.matrix),])
colnames(prediction.L23) = colnames(decon)

prediction.L5 = as.data.frame(SLC.pred[[6]][rownames(L5.matrix),])
colnames(prediction.L5) = colnames(decon)


# prediction for L23
pred = ggplot(prediction.L23, aes(L23.coord$array_row, L23.coord$array_col)) + geom_point(aes(color = L2.3.IT), size = 1 , alpha = 10) + theme_classic() +
  scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0,8)) + ggtitle("predicted SLC9A6 expression in L23")
grid.arrange(L23.exp, pred, nrow = 1)


# prediction for L5
pred = ggplot(prediction.L5, aes(L5.coord$array_row, L5.coord$array_col)) + geom_point(aes(color = L5.IT), size = 1 , alpha = 10) + theme_classic() +
  scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0,8)) + ggtitle("predicted SLC9A6 expression in L5")
grid.arrange(L5.exp, pred, nrow = 1)



# EVALUATE COR RESULTS
# after run, add avg correlation of two cell types
cor.result = cor.result %>% mutate(avg.cor = rowMeans(cor.result[ ,c("result.2", "result.3")]) )
ggplot(cor.result, aes(x = n, y = avg.cor)) + geom_point()

# find n that achieves the best cor result
best_n = cor.result[cor.result$avg.cor == max(cor.result$avg.cor),]$n

ggplot(cor.result, aes(x = n, y = avg.cor)) + geom_line() +
  geom_vline(xintercept = best_n,
             linetype="dotted", color = "blue", size=0.7) +
  labs(title = paste0("Gene: Tbcc, best correlation at n=", best_n ))


# BENCHMARK 1: overall exp as each of the cell type specific expression
bench1 = matrix(rep(spatial[,gene_of_choice], n_ct), ncol = n_ct)
colnames(bench1) = colnames(norm.E)
rownames(bench1) = rownames(norm.E)

 for(i in c(1:n_ct)){
      eval <- eval_correlation_no_na(bench1, norm.E, ct.order = i)
      print(eval)
 }



# plot 100 genes gold standard: find spatially variable genes

setwd("genes")
overlap = readRDS("L23.L5.genes.rds")
overlap.100 = overlap[1:100]

L23.coord <- all.coord[rownames(L23.matrix),]
gold.L23 = as.data.frame(L23.matrix)

gene_list = overlap.100
setwd("../spatial.pattern")

for(gene in gene_list){
  print(gene)
  png(filename = paste0(gene, ".L23.gold.stand.png"))
  p1 = ggplot(gold.L23, aes(L23.coord$array_row, L23.coord$array_col)) + geom_point(data = gold.L23, aes_string(color = gene), size = 1 , alpha = 10) + theme_classic() +
        scale_color_gradientn(colours = rainbow(n = 8, start=0.15, alpha = 0.5), limits = c(0, 8)) + ggtitle(paste0("gold standard expression in L23 for gene: ", gene) )
  print(p1)
  dev.off()
}








