#!/bin/bash
#SBATCH -o facebook_netw.out
#SBATCH -e facebook_netw.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=72:55:00
#SBATCH --mem=16000M

module load r/3.4.3

Rscript facebook_netw_analysis.R

