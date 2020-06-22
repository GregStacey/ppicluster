#!/bin/bash
#SBATCH -o calculate_outcomes_MCP.out
#SBATCH -e calculate_outcomes_MCP.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=02:55:00
#SBATCH --mem=4G
#SBATCH --array=1-56
#SBATCH --account=rrg-ljfoster-ab

module load gcc/7.3.0
module load r/3.6.0
source ~/ENV3/bin/activate

Rscript calculate_outcomes_MCP.R $SLURM_ARRAY_TASK_ID