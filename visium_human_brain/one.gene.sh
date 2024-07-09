#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err
#SBATCH -c 16
#SBATCH --mem=500G
#SBATCH --partition=scavenger
Rscript one.gene.R
