#!/bin/bash
#PBS -N full_netw_analysis_MCPrevisions_large
#PBS -l walltime=71:55:00
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M richard.greg.stacey@gmail.com
#PBS -l mem=64gb   

#SBATCH -o full_netw_analysis_MCPrevisions_large.out
#SBATCH -e full_netw_analysis_MCPrevisions_large.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=71:55:00
#SBATCH --mem=64G
#SBATCH --array=1-216
#SBATCH --account=rrg-ljfoster-ab

module load gcc/9.1.0
module load python/3.7.3
module load r/3.6.2-py3.7.3
source ~/environments/ppicluster_venv/bin/activate

Rscript full_netw_analysis_MCPrevisions.R $SLURM_ARRAY_TASK_ID large
