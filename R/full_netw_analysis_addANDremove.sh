#!/bin/bash
#SBATCH -o full_netw_addANDremove.out
#SBATCH -e full_netw_addANDremove.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=16000M

module load r/3.4.3

Rscript full_netw_analysis_addANDremove.R

