#!/bin/bash
#SBATCH -o slurm.out
#SBATCH -e slurm.err
#SBATCH -c 16
#SBATCH --mem=500G
#SBATCH --partition=scavenger
Rscript run.glu.R $1 >& logs/$1.txt
