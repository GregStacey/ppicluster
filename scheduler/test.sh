#!/bin/bash
#PBS -N full_netw_analysis_MCPrevisions_large
#PBS -l walltime=00:01:00
#PBS -e /scratch/st-ljfoster-1/ppicluster/scheduler/test.err
#PBS -o /scratch/st-ljfoster-1/ppicluster/scheduler/test.out
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M richard.greg.stacey@gmail.com
#PBS -J 0-201:3
#PBS -A st-ljfoster-1

Rscript test.R
