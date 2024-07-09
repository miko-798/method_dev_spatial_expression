## FUNCTIONS

# Function to simulate spot by cell type expression for a particular gene: 
# *** cell type in the SAME order as decon matrix from SPOTlight ***
# norm.astro: these are pre-simulated, spot by gene matrix for each cell type
norm_E_simulate <- function(gene_of_choice, n_spot, n_ct, norm.astro, norm.oligo, norm.macro, norm.glu, norm.gaba){
  
  # initialization
  norm.E <- matrix(NA, nrow = n_spot, ncol = n_ct)
  colnames(norm.E) <- c("Astro",  "GABAergic", "Glutamatergic", "Macrophage", "Oligo")
  rownames(norm.E) <- c(1:n_spot)
  
  rows <- rownames(norm.E)[rownames(norm.E) %in% rownames(norm.astro)]
  norm.E[rows,1] <- norm.astro[rows, gene_of_choice]
  
  rows <- rownames(norm.E)[rownames(norm.E) %in% rownames(norm.gaba)]
  norm.E[rows,2] <- norm.gaba[rows, gene_of_choice]
  
  rows <- rownames(norm.E)[rownames(norm.E) %in% rownames(norm.glu)]
  norm.E[rows,3] <- norm.glu[rows, gene_of_choice]
  
  rows <- rownames(norm.E)[rownames(norm.E) %in% rownames(norm.macro)]
  norm.E[rows,4] <- norm.macro[rows, gene_of_choice]
  
  rows <- rownames(norm.E)[rownames(norm.E) %in% rownames(norm.oligo)]
  norm.E[rows,5] <- norm.oligo[rows, gene_of_choice]
  
  return(norm.E)
}


### Create initialization E values for each cell type
# find spots with high abundance of cell type, to predict the rest
create_init_E_ct <- function(gene_of_choice, decon_input, ct_index){
  
  decon_input <- as.data.frame(decon_input)
  # create fitting data
  spatial.coord <- as.data.frame(cbind(spatial[,gene_of_choice], all.coord))
  colnames(spatial.coord)[1] <- gene_of_choice
  
  high.ct.index <- which(decon_input[,as.numeric(ct_index)] > 0.5) # spots where more than half is this cell type
  #print(length(high.ct.index))
  fit.spot <- spatial.coord[high.ct.index,]
  loess.fit <- loess(fit.spot[,1] ~ x_coor + y_coor, data = fit.spot, control = loess.control(surface = "direct"))
  ## predict a larger set
  ct.init <- predict(loess.fit, newdata = all.coord)
  ct.init[which(ct.init < 0)] <- 0 # get rid of negative numbers
  return(ct.init)
}

# m: take the nearest m spots from each spot to calculate spatial continuity
min.obj.fn <- function(exp_spot_celltype, lambda, n, m, n_spot, n_ct){

  #print(paste0("decon_input: ",length(decon_input)))
  #print(paste0("non_zero_prop: ",length(non_zero_prop)))
  #print(paste0("exp_spot_celltype: : ",length(exp_spot_celltype)))

  # get the entire matrix
  input_exp <-  matrix(0, nrow = nrow(decon_input), ncol = ncol(decon_input))
  #print(paste0("input_exp: ", dim(input_exp)))

  input_exp[non_zero_prop] <- exp_spot_celltype
  #print(paste0("after assign, input_exp: ",length(input_exp)))

  # PART 1
  # predicted gene expression for each spot
  predicted <- rowSums(decon_input*input_exp)
  # residual sum of squares
  rss = mean((zp - predicted)^2)

  # PART 2
  # difference in expression between two spots, across all cell types (# rowSums: all cell type)
  exp_diff <- rowSums((input_exp[df.nearest.m[,1],]-input_exp[df.nearest.m[,2],])^2)
  spatial_cont <- lambda * mean( exp_diff / df.nearest.m[,3]^n )

  return(rss + spatial_cont)
}



# spot by gene matrix, gene in the cols
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

# gradient function
grr <- function(exp_spot_celltype, lambda, n, m, n_spot, n_ct){
  # get the entire matrix
  input_exp <-  matrix(0, nrow = nrow(decon_input), ncol = ncol(decon_input))
  input_exp[non_zero_prop] <- exp_spot_celltype
  
  gradient_term1 <- -2*decon_input*(zp-rowSums(decon_input*input_exp))
  #print(class(gradient_term1))
  #print(class(non_zero_prop))
  
  gradient_term1 <- gradient_term1[non_zero_prop]
  #print(dim(gradient_term1))
  gradient_term1 <- gradient_term1/nrow(decon_input)  

  term2 <- (input_exp[df.nearest.m[,1],] - input_exp[df.nearest.m[,2],]) / df.nearest.m[,3]^n
  
  term2c1 <- rowsum(term2, group =  df.nearest.m[,1])
  term2c2 <- rowsum(term2, group =  df.nearest.m[,2])  ##
  ## term2 is the addition of term2c1 and term2c2, but need to align their names
  
  gradient_term2 <-  matrix(0, nrow = nrow(decon_input), ncol = ncol(decon_input))
  gradient_term2[as.numeric(rownames(term2c1)),] <- term2c1
  gradient_term2[as.numeric(rownames(term2c2)),] <- gradient_term2[as.numeric(rownames(term2c2)),] - term2c2
  
  gradient_term2 <- 2 * lambda * gradient_term2
  #print(dim(gradient_term2))
  gradient_term2 <- gradient_term2[non_zero_prop]
  #print(dim(gradient_term2))
  gradient_term2 <- gradient_term2/nrow(df.nearest.m) # (n_ct*m)

  # return a spot x cell type matrix: the two terms both should be such a matrix
  return(gradient_term1 + gradient_term2)
}


### NO PARALLEL
run_optim_no_parallel <- function(exp_spot_celltype, lambda, n, m, decon_input, zp, ep, n_cpu){
  #print("0")

  cl <- makeCluster(n_cpu) # set the number of processor cores
  setDefaultCluster(cl=cl) # set 'cl' as default cluster
  
  # export objects for parallel process to view them
  clusterExport(cl=cl, c('decon_input', 'zp', 'df.nearest.m', 'non_zero_prop'))
  #print("1")
  
  set.seed(123)
  result <- optim(par = exp_spot_celltype,
                        fn = min.obj.fn,
                        gr = grr,
                        lambda = lambda,
                        n = n,
                        m = m,
			n_spot = n_spot,
			n_ct = n_ct,
                        method = "L-BFGS-B",
                        lower = 0,
                        upper = Inf,
                        control = list(maxit = 100))
  
  setDefaultCluster(cl=NULL); stopCluster(cl)
  #print("2")
  
  # get the entire matrix
  output_exp <-  matrix(0, nrow = nrow(decon_input), ncol = ncol(decon_input))
  output_exp[non_zero_prop] <- result$par
  
  return(list(result, output_exp))
}



### eval method : eval only on the spots containing the cell type
### UPDATED: 8.10, NO SMOOTHING
# RETURN: correlation for that gene

eval_correlation <- function(output_exp, norm.E, norm_cell_type, pure_ct_index){
  # cell type coord
  output_exp <- output_exp[rownames(norm_cell_type),]
  output_exp <- as.matrix(output_exp)

  norm.E <- as.data.frame(norm.E)
  norm.E <- norm.E[rownames(norm_cell_type),]
  norm.E <- as.matrix(norm.E)
  
  if( length(pure_ct_index) != 0 ){
    cor_ct <- colCors(norm.E[-pure_ct_index,], output_exp[-pure_ct_index,])
  }
  else(
    cor_ct <- colCors(norm.E, output_exp)
  )
  return(cor_ct)
}


# Evaluate correlation between method output and gold standard, using all the common rows that are non-NA
eval_correlation_no_na <- function(output_exp, norm.E, ct.order){

  # remove rows where the ct is NA
  norm.E <- norm.E[!is.na(norm.E[,ct.order]) ,]
  output_exp <- output_exp[!is.na(output_exp[,ct.order]) ,]

  # find common rows between gold standard and method output
  a <- rownames(output_exp)
  b <- rownames(norm.E)
  norm.E <- norm.E[intersect(a,b),]
  output_exp <- output_exp[intersect(a,b),]

  # find cor
  output_exp <- as.matrix(output_exp)
  norm.E <- as.matrix(norm.E)

  cor(norm.E[,ct.order], output_exp[,ct.order])

}




# main function for optimization
# returns:
# 1) output value from optimParallel;
# 2) predicted cell type expression
run_optim <- function(exp_spot_celltype, lambda, n, m, decon_input, zp, ep, n_cpu){
  #print("0")

  ## Use optimParallel
  cl <- makeCluster(n_cpu) # set the number of processor cores
  setDefaultCluster(cl=cl) # set 'cl' as default cluster

  # export objects for parallel process to view them
  clusterExport(cl=cl, c('decon_input', 'zp', 'df.nearest.m', 'non_zero_prop'))
  #print("1")

  set.seed(123)
  result <- optimParallel(par = exp_spot_celltype,
                          fn = min.obj.fn,
                          gr = grr,
                          lambda = lambda,
                          n = n,
                          m = m,
			  n_spot = n_spot,
                          n_ct = n_ct,
                          method = "L-BFGS-B",
                          lower = c(0.001,0.001),
                          upper = 10000,
                          control = list(maxit = 100),
                          sleep=0,
                          verbose=TRUE,
                          parallel=list(loginfo=TRUE))




  setDefaultCluster(cl=NULL); stopCluster(cl)
  #print("2")

  # get the entire matrix
  output_exp <-  matrix(0, nrow = nrow(decon_input), ncol = ncol(decon_input))
  output_exp[non_zero_prop] <- result$par

  return(list(result, output_exp))
}
## CHECK GRADIENT IS CORRECT
# grr.result = grr(exp_spot_celltype, lambda, n, m, n_spot, n_ct)
# fn.result = grad(func = min.obj.fn, x=init.E, lambda, n, m, n_spot, n_ct, method="Richardson", side=NULL, method.args=list())
# summary(abs(grr.result-fn.result))



