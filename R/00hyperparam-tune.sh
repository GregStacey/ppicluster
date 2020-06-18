#!/bin/bash
#SBATCH -o 00hyperparam-tune.out
#SBATCH -e 00hyperparam-tune.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=02:55:00
#SBATCH --mem=8G
#SBATCH --array=1-150
#SBATCH --account=rrg-ljfoster-ab

module load gcc/7.3.0
module load r/3.6.0

Rscript 00hyperparam-tune.R $SLURM_ARRAY_TASK_ID