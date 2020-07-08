#!/bin/bash
#SBATCH -o cofrac_MCPrevisions.out
#SBATCH -e cofrac_MCPrevisions.err
#SBATCH --mail-user=richard.greg.stacey@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=23:55:00
#SBATCH --mem=4G
#SBATCH --array=1-56
#SBATCH --account=rrg-ljfoster-ab

module load gcc/7.3.0
module load r/3.6.0
source ~/ENV3/bin/activate

Rscript cofrac_MCPrevisions.R