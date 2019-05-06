#!/bin/bash
#SBATCH -o supFig01.out
#SBATCH -e supFig01.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=02:55:00
#SBATCH --mem=16000M

module load r/3.4.3

Rscript suppFigure01.R

