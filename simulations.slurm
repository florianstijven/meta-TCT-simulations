#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=36 --mem-per-cpu=5000M
#SBATCH --time=00:00:15
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module purge
module load GSL
module load CMake
module load R

Rscript simulation-study.R



