#!/bin/bash
#SBATCH -o implied_edges.out
#SBATCH -e implied_edges.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=4000M
#SBATCH --array=1-27

module load r/3.4.3

Rscript implied_edges.R $SLURM_ARRAY_TASK_ID

