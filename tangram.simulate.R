library(parallel)
library(grDevices)
library(ggplot2)
suppressMessages(library(dplyr))

norm.sim.spatial <- readRDS("norm.sim.spatial.rds")
norm.real.spatial <- readRDS("norm.real.spatial.rds")
all.coord <- readRDS("all.coord.rds")

#test.spatial <- norm.sim.spatial[,1:2]
#test.coord <- all.coord.rds[,1:2]
#print(dim(test.spatial))
#print(dim(test.coord))

#test.spatial <- mclapply(colnames(test.spatial), function(x){
#   fitted(loess(test.spatial[, x] ~ test.coord$x_coor + test.coord$y_coor, test.spatial))}, mc.cores=8)
# smoothed expression
#test.spatial <- matrix(unlist(test.spatial), ncol = length(test.spatial))
#saveRDS(test.spatial, file = "test.after.smooth.rds")

# smooth curve on simulation data, using multi core
#norm.sim.spatial <- mclapply(colnames(norm.sim.spatial), function(x){
#   fitted(loess(norm.sim.spatial[, x] ~ all.coord$x_coor + all.coord$y_coor, norm.sim.spatial))}, mc.cores=16)
# smoothed expression
#norm.sim.spatial <- matrix(unlist(norm.sim.spatial), ncol = length(norm.sim.spatial))
#saveRDS(norm.sim.spatial, file = "after.smooth.rds")

# smooth curve on real data, using multi core
#norm.real.spatial <- mclapply(colnames(norm.real.spatial), function(x){
#   fitted(loess(norm.real.spatial[, x] ~ all.coord$x_coor + all.coord$y_coor, norm.real.spatial))}, mc.cores=16)
# smoothed expression
#norm.real.spatial <- matrix(unlist(norm.real.spatial), ncol = length(norm.real.spatial))
#saveRDS(norm.real.spatial, file = "after.smooth.real.rds")


smooth.sim.spatial <- readRDS("after.smooth.rds")
smooth.real.spatial <- readRDS("after.smooth.real.rds")
colnames(smooth.sim.spatial) <- colnames(norm.sim.spatial)
colnames(smooth.real.spatial) <- colnames(norm.real.spatial)

# Gene Correlation 
###################
# use spearman correlation

smooth.sim.spatial <- as.data.frame(smooth.sim.spatial)
smooth.real.spatial <- as.data.frame(smooth.real.spatial)

smooth.cor <- mapply(function(x,y) cor(x,y,method = "spearman"), as.list(smooth.sim.spatial), as.list(smooth.real.spatial))

smooth.cor <- as.data.frame(smooth.cor)

saveRDS(smooth.cor, file = "smooth.cor.rds")

options(bitmapType='cairo')
png("gene.cor.png")
ggplot(smooth.cor, aes(y=smooth.cor)) + geom_histogram(binwidth = 0.05, col="grey", fill="green", alpha = .2) + 
      labs(title="Gene correlation after smoothing",
           x="Number of genes", y="Correlation") + theme(axis.text=element_text(size=18),
                                                         axis.title=element_text(size=18),
                                                         title = element_text(size = 18)) +
      ylim(-0.5, 1)
dev.off()

high.cor.gene <- smooth.cor %>% filter(smooth.cor > 0.5)
#high.cor.gene

canonical_gene <- c("alc1a3", "aldh1l1", "apoe", "aqp4", "gfap", "gja1", "glul", 
                    "ndrg2", "sox9", "vim", "aldoc", "npy", "rbfox3", "snap25", 
                    "sncb", "stmn2", "tubb3", "satb2", "slc17a7", "slc17a6",
                    "mog", "mag", "omg", "plp1", "sox10", "sox8", "car8", "mbp",
                    "mobp", "olig2", "olig1", "cspg4", "pcdh15", "pdgfra", "sox6",
                    "adgrg1","cd163", "f13a1", "ifitm2", "ifitm3", "tagln2", "tgfbi",
                    "cldn5", "vtn", "naaa", "nkain4", "aldoc", "gja1", "slc1a2", "slc1a3",
                    "tmem119", "itgam","tgfbi")

# add to smooth.cor a col: indicate whether it's a canonical gene
new.col <- rownames(smooth.cor) %in% canonical_gene
smooth.cor <- smooth.cor %>% mutate(canonical = new.col)

# Boxplot: correlation of canonical and non-canonical genes
options(bitmapType='cairo')
png("boxplot.canonical.png", res = 80)
ggplot(smooth.cor, aes(x=as.factor(canonical), y=smooth.cor)) + 
    geom_boxplot(fill="slateblue", alpha=0.2) + 
    labs(title="Gene correlation",
        x="Canonical", y="Correlation")
dev.off()








