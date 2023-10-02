#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=18 --mem-per-cpu=5000M
#SBATCH --time=00:20:00
#SBATCH --cluster=genius
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module unuse /apps/leuven/skylake/2018a/modules/all
module use /apps/leuven/skylake/2021a/modules/all

module purge
module load GSL
module load CMake
module load R

Rscript simulation-study.R 18



