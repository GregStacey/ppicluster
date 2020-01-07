#!/bin/bash
#SBATCH -o test_clust_pertub.out
#SBATCH -e test_clust_pertub.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=16000M

module load r/3.4.3

Rscript test_clust_pertub.R

