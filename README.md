# Meta TCT: Simulations and Application

This repo contains the code associated with Stijven et al. (2024). In this
paper, the theoretical and empirical properties of meta time-component tests
(TCT) are studied. Below, the structure of this repo is explained. 

This project relies on the `renv` package for managing R packages and ensuring a
reproducible environment. When opening this project in RStudio, `renv` is
automatically initialized and will suggest to install/update packages. To ensure
reproducibility, it is advised to follow this suggestion.

# Simulations

## Data Generation and Analysis

The data are simulated and analyzed in `simulation-study.R`. Because the
simulations are computer intensive, this script was run on a high-performance
computing (HPC) cluster, provided by the VSC (Vlaams Supercomputer Centrum). The
`simulations.slurm` file is the job script that was used to run `simulation-study.R`
on the HPC cluster; the corresponding console output was saved to
`slurm-61418490.out`. The results of the simulation study were saved to
`results_simulation_lean.rds`. This files is not present in this repo because it
is too large, but it can be downloaded from [here](https://kuleuven-my.sharepoint.com/:u:/g/personal/florian_stijven_kuleuven_be/EWUgNdCHNnpMsjRvngtcfIYBavy6Gjvhh3FPaf2Dk4yoxA?e=7MNmG5).

## Results

The figures presented in Stijven et al. (2024) can be reproduced by running
`figures-and-tables.R`; the figures will then automatically be saved into
`figures\`. To reproduce these figures `results_simulation_lean.rds` needs to be
present in the repo; hence, it should first be downloaded (using the link above)
and saved in this repo's main directory.

# Application

The `application\real-data-analysis\` directory contains an RMarkdown file (and
knitted html files) in which the Lighthouse data are analyzed. The corresponding
results are presented in Stijven et al. (2024). In these files, every step in
the analysis is explained - from the data exploration to the interpretation of
the meta-TCT estimates. The original data cannot be shared and are, therefore,
not present in this repo.

The `application\mock-analysis\` directory contains the same files as
`application\real-data-analysis\` and, additionally, a simulated data set. The
simulation of these data was based on the Lighthouse data. These files serve as
a practical demonstration of the meta TCT estimators implemented in the `TCT` R
package.

# References

Stijven F, Mallinckrodt C, Molenberghs G, Alonso A, Dickson S, Hendrix S. (2024) 
Meta TCT: Towards Interpretable Treatment Effects in Clinical Trials for 
Progressive Diseases
