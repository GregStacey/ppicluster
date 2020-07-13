#!/bin/bash
#SBATCH -o full_netw_analysis_MCPrevisions_large.out
#SBATCH -e full_netw_analysis_MCPrevisions_large.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=71:55:00
#SBATCH --mem=64G
#SBATCH --array=1-216
#SBATCH --account=rrg-ljfoster-ab

module load gcc/7.3.0
module load r/3.6.0
source ~/ENV3/bin/activate

Rscript full_netw_analysis_MCPrevisions.R $SLURM_ARRAY_TASK_ID large
