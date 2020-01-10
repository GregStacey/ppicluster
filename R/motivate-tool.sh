#!/bin/bash
#SBATCH -o motivate-tool.out
#SBATCH -e motivate-tool.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=4000M

module load r/3.4.3

Rscript motivate-tool.R

