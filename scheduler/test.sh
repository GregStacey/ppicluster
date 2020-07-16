#!/bin/bash
#PBS -N full_netw_analysis_MCPrevisions_large
#PBS -l walltime=00:01:00
#PBS -e /home/staceyri/projects/ppicluster/scheduler/full_netw_analysis_MCPrevisions_large.err
#PBS -o /home/staceyri/projects/ppicluster/scheduler/full_netw_analysis_MCPrevisions_large.out
#PBS -j oe
#PBS -k eod
#PBS -m abe
#PBS -M richard.greg.stacey@gmail.com
#PBS -J 0-12:3
#PBS -A st-ljfoster-1

Rscript test.R
