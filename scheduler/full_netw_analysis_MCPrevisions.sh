#!/bin/bash
#PBS -N full_netw_analysis_MCPrevisions_large
#PBS -l walltime=71:55:00
#PBS -e /scratch/st-ljfoster-1/ppicluster/scheduler/full_netw_analysis_MCPrevisions_large.err
#PBS -o /scratch/st-ljfoster-1/ppicluster/scheduler/full_netw_analysis_MCPrevisions_large.out
#PBS -m abe
#PBS -M richard.greg.stacey@gmail.com
#PBS -l mem=64gb   
#PBS -J 1-36
#PBS -A st-ljfoster-1

module load gcc/9.1.0
module load python/3.7.3
module load r/3.6.2-py3.7.3
source ~/environments/ppicluster_venv/bin/activate

Rscript full_netw_analysis_MCPrevisions.R $PBS_ARRAYID large
