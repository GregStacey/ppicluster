#!/bin/bash
#PBS -N full_netw_analysis_MCPrevisions_large
#PBS -l walltime=00:01:00
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M richard.greg.stacey@gmail.com
#PBS –t 201   

Rscript test.R
