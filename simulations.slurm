#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=36 --mem-per-cpu=5000M
#SBATCH --time=00:10:00
#SBATCH --cluster=wice
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1
module list


Rscript simulation-study.R 35



