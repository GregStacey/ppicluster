#!/bin/bash
#SBATCH -o supFig03.out
#SBATCH -e supFig03.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=16000M

module load r/3.4.3

Rscript suppFigure02.R

