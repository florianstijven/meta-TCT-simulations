#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=96
#SBATCH --time=1:00:00
#SBATCH --cluster=wice
#SBATCH --partition=batch_sapphirerapids
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

module use /apps/leuven/rocky8/sapphirerapids/2024a/modules/all
module load GSL
module load CMake
module load R-bundle-CRAN

make simulation



