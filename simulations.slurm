#!/bin/bash
#
#SBATCH --job-name=meta-TCT-simulations
#SBATCH --ntasks=1 --cpus-per-task=96
#SBATCH --time=48:00:00
#SBATCH --cluster=wice
#SBATCH --partition=batch_sapphirerapids
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -A lp_doctoralresearch

export OMP_NUM_THREADS=1

module load cluster/wice/batch_sapphirerapids
module load R/4.4.2-gfbf-2024a
module load R-bundle-CRAN/2024.11-foss-2024a

make simulation



