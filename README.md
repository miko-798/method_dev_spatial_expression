# Method development for deconvoluting cell-type-specific spatial gene expression pattern  

Project at Duke University, 2020-2021

**Goal:**  
Deconvolute cell-type-specific expression from spatial gene expression data

**Input:**  
Spatial gene expression + Cell type proportions for each spot

**Output:**  
Cell-type-specific gene expression for each spot 

<img width="901" alt="Screenshot 2024-07-14 at 6 06 25 PM" src="https://github.com/user-attachments/assets/3a0c5ec0-83ca-4c42-81a5-25dc53a10033">

<br>
  
**Spatial continuity assumption:** Spatially neighboring cells with the same cell type have similar gene expression. 

We designed an objective function to recover the observed overall gene expression aggregating the cell-type-specific gene expression, while following the spatial continuity assumption. 

We used cross validation approach to select the hyperparameters in the objective function. optimParallel method in R was used for optimizing the objective function. 

<img width="913" alt="Screenshot 2024-07-14 at 6 07 15 PM" src="https://github.com/user-attachments/assets/64cf8407-abb7-4e9c-92fe-e3a6b162dec5">


**Hyperparameters:**  
n: balance between local and global expression continuity, controls how similar the expression should be for adjacent spots.  
Lambda: balance the two terms in the objective function.

Some example code scripts are included in this repo. 
