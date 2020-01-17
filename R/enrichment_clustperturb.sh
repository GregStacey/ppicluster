#!/bin/bash
#SBATCH -o enrichment_clustperturb.out
#SBATCH -e enrichment_clustperturb.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=71:55:00
#SBATCH --mem=16000M

module load r/3.4.3

Rscript enrichment_clustperturb.R

