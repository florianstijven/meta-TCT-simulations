#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=18 --mem-per-cpu=5000M
#SBATCH --time=00:08:00
#SBATCH --cluster=genius
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module load R

Rscript simulation-study.R 17



