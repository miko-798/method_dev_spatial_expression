require(tidyverse)
library(ggplot2)

source("functions.R")

setwd("results/")
df.cor <- list.files(pattern = "*.cor.result.rds") %>%
  map_dfr(readRDS)

dim(df.cor)
summary(df.cor$result)

# order:  "Astro"     "GABAergic"     "Glutamatergic" "Macrophage"    "Oligo"  
colnames(df.cor)[3:7] = c("astro", "gaba", "glu", "macro", "oligo")

# for 33 genes
df.cor = df.cor[df.cor$gene %in% prom.gene.list, ]


# create df :
# first col is all the correlation result
# second col is ct labels

cor.result = c(df.cor[,3], df.cor[,4], df.cor[,5], df.cor[,6], df.cor[,7])


ct.label <- c(rep("astro", nrow(df.cor)),
              rep("gaba", nrow(df.cor)),
              rep("glu", nrow(df.cor)),
              rep("macro", nrow(df.cor)),
              rep("oligo", nrow(df.cor)))


## LENGTH
length(cor.result) == length(ct.label)

ct.cor.result <- data.frame("cor" = cor.result,
                            "ct" = ct.label)


# predicted vs. gold standard
ggplot(ct.cor.result, aes(x=ct, y=cor, fill = ct)) + theme_minimal() +
  geom_boxplot() + xlab("Cell type") + ylab("Correlation") + ylim(0,1) + 
  ggtitle("Correlation results of in each cell type")


###
# calculate benchmark cor: predicted vs. visium 

n_ct = 5
# input: predicted matrix, spatial data (col repeated 5 times)

setwd("../data/")
spatial = readRDS("sim.spatial.rds")
n_spot = nrow(spatial)

setwd("../genes")
#genes_of_interest = readRDS("genes.50.rds")
genes_of_interest = prom.gene.list

setwd("../results/")

# Initialize the correlation result table
bench_cor = data.frame(gene = genes_of_interest, n = 10, 
		       result = matrix(rep(0, n_ct*length(genes_of_interest)), ncol = 5))
colnames(bench_cor)[3:7] = c("astro", "gaba", "glu", "macro", "oligo")


# find cor
for(gene in genes_of_interest){

print(gene)
predicted = paste0(gene, ".predict.result.rds")

benchmark1 = spatial[,gene]
benchmark1_m = matrix(rep(benchmark1, 5), ncol = 5)

output_exp = readRDS(predicted)[[1]]
rownames(benchmark1_m) = rownames(spatial)
rownames(output_exp) = rownames(spatial)


# benchmark cor result
# order:  "Astro"     "GABAergic"     "Glutamatergic" "Macrophage"    "Oligo"  
    for(i in c(1:n_ct)){
      eval <- eval_correlation_no_na(output_exp, benchmark1_m, ct.order = i)
      print(eval)
      # update results
      bench_cor[bench_cor$gene == gene, 2+i]= eval
    }
}


# create df :
# first col is all the correlation result
# second col is ct labels

bench.cor.result = c(bench_cor[,3], bench_cor[,4], bench_cor[,5], bench_cor[,6], bench_cor[,7])


bench.ct.label <- c(rep("astro", nrow(bench_cor)),
              rep("gaba", nrow(bench_cor)),
              rep("glu", nrow(bench_cor)),
              rep("macro", nrow(bench_cor)),
              rep("oligo", nrow(bench_cor)))


## LENGTH
length(bench.cor.result) == length(bench.ct.label)

plot.bench.cor.result <- data.frame("cor" = bench.cor.result,
                           	    "ct" = bench.ct.label)


# predicted vs. gold standard
ggplot(plot.bench.cor.result, aes(x=ct, y=cor, fill = ct)) + theme_minimal() +
  geom_boxplot() + xlab("Cell type") + ylab("Correlation") + ylim(0,1) +
  ggtitle("Benchmark 1: correlation results of in each cell type with spatial data")

# Boxplot side by side comparison
plot.bench.cor.result[plot.bench.cor.result$ct == "astro",]$ct = "astro.bench"
plot.bench.cor.result[plot.bench.cor.result$ct == "gaba",]$ct = "gaba.bench"
plot.bench.cor.result[plot.bench.cor.result$ct == "glu",]$ct = "glu.bench"
plot.bench.cor.result[plot.bench.cor.result$ct == "macro",]$ct = "macro.bench"
plot.bench.cor.result[plot.bench.cor.result$ct == "oligo",]$ct = "oligo.bench"

# combine benchmark result and method result
comb.result = rbind(plot.bench.cor.result, ct.cor.result)

# Plot one boxplot comparing the two
##################################
ggplot(comb.result, aes(x=ct, y=cor, fill = ct)) + theme_minimal() +
  geom_boxplot() + xlab("Cell type") + ylab("Correlation") + ylim(0,1) +
  ggtitle("Method vs. benchmark correlation results of in each cell type with spatial data")


#######
# Select genes highly expressed in both glu and gaba
setwd("genes")
glu = readRDS("glu.cv.rds")
gaba = readRDS("gaba.cv.rds")
genes.glu.gaba = intersect(names(glu), names(gaba))



# TODO: why cor is so high?? for benchmark 1?

# calculate cor of genes zp with cell type expression

setwd("../data")

norm.astro = readRDS("norm.astro.rds")
norm.oligo = readRDS("norm.oligo.rds")
norm.macro = readRDS("norm.macro.rds")
norm.glu = readRDS("norm.glu.rds")
norm.gaba = readRDS("norm.gaba.rds")


# TODO: save five of them as ONE
genes_of_interest = genes.glu.gaba

cor.zp.gold = data.frame(gene = genes_of_interest, cor= matrix(rep(0, n_ct*length(genes_of_interest)), ncol = 5))
colnames(cor.zp.gold)[2:6] = c("astro", "gaba", "glu", "macro", "oligo")


for (gene in genes_of_interest){
for (i in 1:n_ct){
  zp = spatial[, gene]
  norm.E <- norm_E_simulate(gene, spatial, n_spot, n_ct, norm.astro, norm.oligo, norm.macro, norm.glu, norm.gaba)
  ct.exp = norm.E[,i]
  index.to.remove = which(is.na(ct.exp))  # NA in norm.E, remove NA
  ct.exp = ct.exp[-index.to.remove]
  zp = zp[-index.to.remove]
  spatial.n.gold = cor(zp, ct.exp)
  cor.zp.gold[cor.zp.gold$gene == gene, 1+i] = spatial.n.gold
  #print(spatial.n.gold)
}

}

# Select genes where glu has a low correlation with overall expression
prom.genes = cor.zp.gold[cor.zp.gold$glu < 0.5,] # 33 genes
cor.zp.gold = prom.genes
prom.gene.list = prom.genes$gene

# plot the above result

# create df :
# first col is all the correlation result
# second col is ct labels

zp.gold.cor.result = c(cor.zp.gold[,2], cor.zp.gold[,3], cor.zp.gold[,4], cor.zp.gold[,5], cor.zp.gold[,6])


ct.label <- c(rep("astro", nrow(cor.zp.gold)),
              rep("gaba", nrow(cor.zp.gold)),
              rep("glu", nrow(cor.zp.gold)),
              rep("macro", nrow(cor.zp.gold)),
              rep("oligo", nrow(cor.zp.gold)))


## LENGTH
length(zp.gold.cor.result) == length(ct.label)
plot.cor.zp.gold <- data.frame("cor" = zp.gold.cor.result,
                               "ct" = ct.label)

ggplot(plot.cor.zp.gold, aes(x=ct, y=cor, fill = ct)) + theme_minimal() +
  geom_boxplot() + xlab("Cell type") + ylab("Correlation") + ylim(0,1) +
  ggtitle("correlation of gene overall expression (simulated visium) vs. gold standard")













