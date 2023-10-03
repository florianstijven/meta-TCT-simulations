#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=36 --mem-per-cpu=5000M
#SBATCH --time=00:20:00
#SBATCH --cluster=genius
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

module load R/4.2.1-foss-2021a

Rscript simulation-study.R 35



