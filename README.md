# Method development for deconvoluting cell-type-specific spatial gene expression pattern. 

Project at Duke University, 2020-2021

Spatial continuity assumption: Spatially neighboring cells with the same cell type have similar gene expression. 

We designed an objective function to recover the observed overall gene expression aggregating the cell-type-specific gene expression, while following the spatial continuity assumption. 

We used cross validation approach to select the hyperparameters in the objective function. optimParallel method in R was used for optimizing the objective function. 

Some example code scripts are included in this repo. 
