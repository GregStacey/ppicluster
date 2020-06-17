#!/bin/bash
#SBATCH -o 00hyperparam-tune.out
#SBATCH -e 00hyperparam-tune.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=023:55:00
#SBATCH --mem=8G
#SBATCH --array=1-150

module load r/3.6
module load gcc/5.4.0

Rscript experiment-architecture-ms.R $SLURM_ARRAY_TASK_ID