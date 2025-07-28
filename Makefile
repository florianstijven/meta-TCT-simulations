.PHONY: all application simulation synthetic-application
all: simulation \
	application-synthetic
	
simulation: results/raw-results/simulations/results_simulation_full.rds \
	R/simulations/simulations.Rout
	
application: results/reports/application/lighthouse-analysis-with-code.html \
	results/reports/application/lighthouse-analysis-wo-code.html
	
synthetic-application: results/reports/synthetic-application/synthetic-application-with-code.html \
	results/reports/synthetic-application/synthetic-application-wo-code.html


application-synthetic: R/application/data-exploration-synthetic.Rout \
	results/raw-results/application-synthetic/ipd_surr_indices_tbl.rds \
	results/raw-results/application-synthetic/ma_trt_effects_tbl.rds \
	R/application/meta_analysis-synthetic.Rout \
	results/raw-results/application-synthetic/bayesian_ma_results.rds \
	R/application/processing-results-synthetic.Rout
	
	
results/raw-results/simulations/results_simulation_full.rds: R/simulations/simulation-study.R
	Rscript R/simulations/simulation-study.R 20 > R/simulations/simulation-study.Rout 2> R/simulations/simulation-study.Rout
	
R/simulations/processing.Rout: R/simulations/processing.R R/simulations/simulation-study.R
	Rscript R/simulations/processing.R > $@ 2> $@
	
	
results/reports/application/lighthouse-analysis-with-code.html: R/application/real-data-analysis/lighthouse-analysis.Rmd
	Rscript -e 'library(rmarkdown); rmarkdown::render(input = "R/application/real-data-analysis/lighthouse-analysis.Rmd", output_file = "lighthouse-analysis-with-code.html", output_dir = "results/reports/application", params = list(printcode = TRUE))'

results/reports/application/lighthouse-analysis-wo-code.html: R/application/real-data-analysis/lighthouse-analysis.Rmd
	Rscript -e 'library(rmarkdown); rmarkdown::render(input = "R/application/real-data-analysis/lighthouse-analysis.Rmd", output_file = "lighthouse-analysis-wo-code.html", output_dir = "results/reports/application", params = list(printcode = FALSE))'



results/reports/synthetic-application/synthetic-application-with-code.html: R/application/real-data-analysis/lighthouse-analysis.Rmd
	Rscript -e 'library(rmarkdown); rmarkdown::render(input = "R/application/mock-analysis/analysis.Rmd", output_file = "synthetic-application-with-code.html", output_dir = "results/reports/synthetic-application", params = list(printcode = TRUE))'

results/reports/synthetic-application/synthetic-application-wo-code.html: R/application/real-data-analysis/lighthouse-analysis.Rmd
	Rscript -e 'library(rmarkdown); rmarkdown::render(input = "R/application/mock-analysis/analysis.Rmd", output_file = "synthetic-application-wo-code.html", output_dir = "results/reports/synthetic-application", params = list(printcode = FALSE))'

