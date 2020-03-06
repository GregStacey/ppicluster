#!/bin/bash
#SBATCH -o Fig08.out
#SBATCH -e Fig08.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=16000M

module load gcc/7.3.0
module load nixpkgs/16.09
module load netcdf/4.6.1
module load r/3.6.0

Rscript Figure08.R