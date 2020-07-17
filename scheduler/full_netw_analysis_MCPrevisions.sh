#!/bin/bash
#PBS -N full_netw_analysis_MCPrevisions_large
#PBS -l walltime=71:55:00
#PBS -e /scratch/st-ljfoster-1/ppicluster/scheduler/5full_netw_analysis_MCPrevisions_large.err
#PBS -o /scratch/st-ljfoster-1/ppicluster/scheduler/5full_netw_analysis_MCPrevisions_large.out
#PBS -m abe
#PBS -M richard.greg.stacey@gmail.com
#PBS -l mem=8gb   
#PBS -A st-ljfoster-1
#PBS -J 1-6


module load gcc/9.1.0
module load python/3.7.3
module load r/3.6.2-py3.7.3
source ~/environments/ppicluster_venv/bin/activate

echo 'the array number is'
echo $PBS_ARRAY_INDEX

Rscript ~/projects/ppicluster/R/full_netw_analysis_MCPrevisions.R \
  ~/projects/ppicluster/data/jobs_huri_louvain.txt \
  $PBS_ARRAY_INDEX
