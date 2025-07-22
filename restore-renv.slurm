#!/bin/bash
#
#SBATCH --job-name=restore-renv
#SBATCH --ntasks=1 --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --cluster=wice
#SBATCH --partition=batch_sapphirerapids
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

module use /apps/leuven/rocky8/sapphirerapids/2024a/modules/all
module load GSL
module load CMake
module load R

Rscript -e "renv::status()" -e "renv::restore(exclude = 'mmrm')" - "renv::status()"