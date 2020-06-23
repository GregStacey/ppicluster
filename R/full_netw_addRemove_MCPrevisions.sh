#!/bin/bash
#SBATCH -o full_netw_addRemove_MCPrevisions.out
#SBATCH -e full_netw_addRemove_MCPrevisions.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=8G
#SBATCH --array=1-324
#SBATCH --account=rrg-ljfoster-ab

module load gcc/7.3.0
module load r/3.6.0
source ~/ENV3/bin/activate

Rscript full_netw_addRemove_MCPrevisions.R $SLURM_ARRAY_TASK_ID
