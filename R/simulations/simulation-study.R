#!/usr/bin/env Rscript
.libPaths()
args = commandArgs(trailingOnly=TRUE)
options(RENV_CONFIG_SANDBOX_ENABLED = FALSE)
# Print R-version for better reproducibility.
print(version)
# Ensure that the state of the R packages is correct.
renv::restore()

# Directory to save results in.
dir_results = here::here("results", "raw-results", "simulations")

# Number of independent replications for each setting. This is specified as a 
# command line argument.
N_trials = as.integer(args[1])

# Load the required packages.
library(TCT)
library(mmrm)
library(dplyr)
library(future)
library(furrr)

# Start tracking time.
a = Sys.time()

# Set up parallel computing. We use the `future` package to enable parallel
# computing. The number of workers is set to the number of cores minus one. By
# default, we try to use multicore processing (which relies on forking, this is
# not possible on Winfows systems). Else, we use multisession which creates N
# separate R processes.
if (parallelly::supportsMulticore()) {
  plan("multicore", workers = parallel::detectCores() - 1)
} else {
  plan(multisession, workers = parallel::detectCores() - 1)
}

# Variance-Covariance matrix as reported by Raket (2022, doi: 10.1002/sim.9581) in
# section 5.1.
vcov_ref = matrix(c(45.1, 40.0, 45.1, 54.9, 53.6, 53.6, 60.8,
                    40.0, 57.8, 54.4, 66.3, 64.1, 64.1, 74.7,
                    45.1, 54.4, 72.0, 80.0, 77.6, 77.6, 93.1, 
                    54.9, 66.3, 80.0, 109.8, 99.3, 99.3, 121.7, 
                    53.6, 64.1, 77.6, 99.3, 111.4, 99.1, 127.8,
                    53.6, 64.1, 77.6, 99.3, 99.1, 111.4, 127.8,
                    60.8, 74.7, 93.1, 121.7, 127.8, 127.8, 191.4), nrow = 7, byrow = TRUE)
# Mean values as reported by Raket (2022): `normal`. The `fast` vector
# represents an additional scenario where progression goes faster and/or the
# study duration is longer.
ref_means_list = list(
  normal = c(19.6, 20.5, 20.9, 22.7, 23.8, 25.8, 27.4),
  fast = c(18, 19.7, 20.9, 22.7, 24.7, 27.1, 29.2)
)

# In `settings_tbl`, we define all settings for which we will generate data. We
# consider all possible combinations of trial settings.
settings_tbl = tidyr::expand_grid(
  progression = c("normal", "fast"),
  gamma_slowing = c(1, 0.9, 0.75, 0.5),
  n = c(50, 200, 500, 1000),
  time_points_chr = c("24 Months", 
                      "36 Months", 
                      "36(-30) Months")
) 


settings_tbl = settings_tbl %>%
  left_join(tibble(
    time_points_chr = c("24 Months",
                        "36 Months",
                        "36(-30) Months"),
    time_points = list(c(0, 6, 12, 18, 24),
                       c(0, 6, 12, 18, 24, 30, 36),
                       c(0, 6, 12, 18, 24, 36))
  ))

# Set the seed for reproducibility. Reproducibility of the random number
# generator in the parallel computations below is handled by the `future` and 
# `furrr` packages.
set.seed(1)


#-----
# We next define a set of helper functions. 

# The analyze_mmrm_new() function fits the MMRM to a given data set using the
# mmrm::mmrm() function. The model fit object returned by the latter function is
# returned by this function as well.
analyze_mmrm_new = function(data_trial, type, reml = TRUE) {
  formula = formula(outcome~arm_time + 0)
  if (type == "null") {
    formula = update.formula(old = formula,
                             .~. - arm_time + as.factor(time_int))
  }
  data_trial$SubjId = as.factor(data_trial$SubjId)
  data_trial$time_int = as.factor(data_trial$time_int)
  formula = update.formula(old = formula, .~. + us(time_int | SubjId))
  mmrm_fit = mmrm::mmrm(
    formula = formula,
    data = data_trial,
    reml = reml,
    control = mmrm::mmrm_control(method = ifelse(reml, "Kenward-Roger", "Satterthwaite"))
  )
  # By default, the actual data set is save into the object returned by
  # mmrm::mmrm(). This causes a significant overhead. We therefore manually
  # remove the data saved into `fit`.
  mmrm_fit$data = NULL
  mmrm_fit$tmb_data = NULL
  mmrm_fit$tmb_object = NULL
  gc()
  return(mmrm_fit)
}

# Function to generate the data.
simulate_data = function(n,
                         time_points,
                         ref_means,
                         trt_means,
                         vcov) {
  # Character vector for naming the various measurement occasions.
  time_names = paste0("month_", time_points)
  
  tibble(as_tibble(
    # Generate data for the control group from a multivariate normal
    # distribution.
    mvtnorm::rmvnorm(
      n = n / 2,
      mean = ref_means,
      sigma = vcov
    ) %>% `colnames<-`(time_names)
  ),
  arm = 0L) %>%
    bind_rows(tibble(# Generate data for the active treatment group.
      as_tibble(
        mvtnorm::rmvnorm(
          n = n / 2,
          mean = trt_means,
          sigma = vcov
        ) %>% `colnames<-`(time_names)
      ),
      arm = 1L)) %>%
    mutate(# Add patient identifier
      SubjId = 1:n) %>%
    # Convert the data, generated in a wide format, to the long format.
    tidyr::pivot_longer(
      cols = all_of(time_names),
      names_to = "time_chr",
      values_to = "outcome"
    ) %>%
    # Add an integer variable indicting the measurement occasion. Also add a
    # categorical variable for each time-arm combination. Except for the
    # baseline measurement occasion.
    mutate(
      time_int = as.integer(stringr::str_remove_all(time_chr, "month_")),
      arm_time = ifelse(
        time_int == 0L,
        "baseline",
        paste0(
          arm,
          ":",
          stringr::str_pad(
            time_int,
            side = "left",
            pad = "0",
            width = 2
          )
        )
      )
    )
}

# Define function where the data are generated from a multivariate normal
# distribution and analyzed with the MMRM.
simulate_and_analyze = function(n,
                                time_points,
                                ref_means,
                                trt_means,
                                vcov, 
                                reml = TRUE) {
  simulated_data = simulate_data(
    n = n,
    time_points = time_points,
    ref_means = ref_means,
    trt_means = trt_means,
    vcov = vcov
  )
  analyze_mmrm_new(simulated_data, type = "full", reml = reml)
}

# Define function where we extract the relevant quantities from the fitted MMRM.
extract_mmrm = function(mmrm_fitted) {
  # The p-value for the F-test with approximate degrees of freedom is computed
  # first.
  # Compute number of post-randomization measurement occasions.
  K = (length(coef(mmrm_fitted)) - 1) / 2
  contrast = matrix(data = 0,
                    nrow = K,
                    ncol = length(coef(mmrm_fitted)))
  contrast[1:K, 1:K] = diag(1, nrow = K)
  contrast[1:K, (K + 1):(2 * K)] = diag(-1, nrow = K)
  p = mmrm::df_md(mmrm_fitted, contrast)$p_val
  
  # The estimated parameters and corresponding estimated variance-covariance
  # matrix are extracted.
  # Compute number of post-randomization measurement occasions.
  estimates = coef(mmrm_fitted)
  vcov_mmrm = vcov(mmrm_fitted)
  # Everyting is returned in a list.
  return(
    list(
      p = p,
      estimates = estimates,
      vcov_mmrm = vcov_mmrm
    )
  )
}

# Define function that handles everyhting from simulating the data to extracting
# the relevant quantities from the MMRM object.
simulate_to_extract = function(n,
                               time_points,
                               ref_means,
                               trt_means,
                               vcov,
                               reml = TRUE) {
  mmrm_fitted = simulate_and_analyze(n,
                                     time_points,
                                     ref_means,
                                     trt_means,
                                     vcov,
                                     reml = TRUE)
  extract_mmrm(mmrm_fitted)
}

#---------

settings_tbl = settings_tbl %>%
  # Group by each row.
  rowwise(everything()) %>%
  # Compute a set of variables that will be useful further on.
  summarise(
    # The number of measurement occasions.
    K = length(time_points),
    # Numbers for the measurement occasions that are included, starting from 1
    # for the first measurement occasion.
    K_vec = list((time_points %/% 6) + 1),
    # Character vector for the naming the various measurement occasions.
    time_names = list(paste0("month_", time_points)),
    # Mean vector in the reference group.
    ref_means = list(ref_means_list[[progression]][K_vec]),
    # Mean vector in the active treatment group.
    trt_means = list(
      spline(
        x = time_points,
        y = ref_means,
        xout = gamma_slowing * time_points,
        method = "natural"
      )$y
    ),
    vcov = list(vcov_ref[K_vec, K_vec]),
  ) %>%
  ungroup()

# Duplicate rows N_trial times. We add a trial identifier.
settings_tbl = settings_tbl %>%
  cross_join(tibble(trial_number = 1:N_trials))

# Fit the MMRM model for each generated data set.
results_tbl = settings_tbl %>%
  mutate(
    mmrm_extracted = future_pmap(
      .l = list(n = n,
                time_points = time_points,
                ref_means = ref_means,
                trt_mean = trt_means,
                vcov = vcov),
      .f = simulate_to_extract,
      .options = furrr_options(
        seed = TRUE,
        stdout = FALSE,
        conditions = character()
      )
    )
  )

# `results_tbl` is not yet in an easy to process format because the mmrm_extracted
# variable is a list of lists. We put the list elements into separate variables.
# We also drop variables that are not used further on.
results_tbl = results_tbl %>%
  select(-time_names, -ref_means, -trt_means, -vcov) %>%
  rowwise(everything()) %>%
  reframe(
    tibble(
      p_value_mmrm = mmrm_extracted$p,
      vcov_mmrm = list(mmrm_extracted$vcov_mmrm),
      coef_mmrm = list(mmrm_extracted$estimates)
    )
  ) %>% # drop list-valued column.
  select(-mmrm_extracted)

print(Sys.time() - a)
print("MMRMs fitted")


# Apply meta-TCT methodology.
results_tbl = results_tbl %>%
  # We first add additional columns that specify the inference options.
  cross_join(
    tidyr::expand_grid(
      drop_first_occasions = 0,
      constraints = c(FALSE),
      interpolation = c("spline", "linear"),
      inference = c("contrast", "least-squares"),
      B = c(0, 5e2)
    )
  ) %>%
  # We only use the bootstrap for the normal progression scenario.
  filter(!(B == 5e2 & (progression == "fast" | interpolation == "linear"))) %>%
  # If for a given setting, we consider bootstrap replications, then we can
  # delete the scenarios corresponding to B == 0 because we would be fitting the
  # same model twice.
  filter(!(B == 0 & !(progression == "fast" | interpolation == "linear")))


# Estimate the time-specific and common acceleration factors based on the
# summary-level information from the MMRMs. We immediately extract the
# estimates, SEs, and CIs we need and leave everything else out to save memory.
results_tbl = future_pmap(
  .l = list(
    coef_mmrm = results_tbl$coef_mmrm,
    vcov_mmrm = results_tbl$vcov_mmrm,
    constraints = results_tbl$constraints,
    interpolation = results_tbl$interpolation,
    K = results_tbl$K,
    time_points = results_tbl$time_points,
    drop_first_occasions = results_tbl$drop_first_occasions,
    inference = results_tbl$inference,
    B = results_tbl$B,
    gamma_slowing = results_tbl$gamma_slowing
  ),
  .f = function(coef_mmrm,
                vcov_mmrm,
                constraints,
                interpolation,
                K,
                time_points,
                drop_first_occasions,
                inference,
                B,
                gamma_slowing) {
    # The TCT_meta() function is called within a TryCatch() expression to
    # prevent the code from failing when one simulation fails. If TCT_meta()
    # fails, an NA is returned.
    TCT_meta_fit = tryCatch(
      TCT_meta(
        time_points = time_points,
        exp_estimates = coef_mmrm[K:(2 * K - 2)],
        ctrl_estimates = coef_mmrm[c(2 * K - 1, 1:(K - 1))],
        vcov = vcov_mmrm[c(2 * K - 1, 1:(K - 1), K:(2 * K - 2)), c(2 * K - 1, 1:(K - 1), K:(2 * K - 2))],
        interpolation = interpolation,
        inference = "wald",
        B = 0,
        constraints = constraints
      ),
      error = function(cond) {
        return(NA)
      }
    )
    if (!is.list(TCT_meta_fit)) {
      return(NA)
    }
    # Estimate common acceleration factor.
    TCT_meta_common_fit = tryCatch(
      TCT_meta_common(
        TCT_Fit = TCT_meta_fit,
        inference = inference,
        B = B,
        select_coef = (drop_first_occasions + 1):length(coef(TCT_meta_fit)),
        constraints = constraints,
        start_gamma = gamma_slowing
      ),
      error = function(cond) {
        return(NA)
      }
    )
    if (!is.list(TCT_meta_common_fit)) {
      return(NA)
    }
    
    # Compute confidence intervals and p-values. First, we need to call the summary
    # method on the object returned by TCT_meta_common(). This method computes the
    # p-value and confidence intervals.
    summary_TCT_common = tryCatch(
      summary(TCT_meta_common_fit),
      error = function(cond) {
        return(NA)
      }
    )
    
    if (!is.list(summary_TCT_common)) {
      return(NA)
    }
    
    return(
      data.frame(
        p_value_TCT_common = summary_TCT_common$p_value,
        estimate = coef(TCT_meta_common_fit),
        conf_int_TCT_common_lower = summary_TCT_common$gamma_common_ci[1],
        conf_int_TCT_common_upper = summary_TCT_common$gamma_common_ci[2],
        se_TCT_common = summary_TCT_common$gamma_common_se,
        se_TCT_common_bs = ifelse(
          is.null(summary_TCT_common$se_bootstrap),
          NA,
          summary_TCT_common$se_bootstrap
        ),
        conf_int_TCT_common_lower_bs = ifelse(
          is.null(summary_TCT_common$se_bootstrap),
          NA,
          summary_TCT_common$ci_bootstrap[1]
        ),
        conf_int_TCT_common_upper_bs = ifelse(
          is.null(summary_TCT_common$se_bootstrap),
          NA,
          summary_TCT_common$ci_bootstrap[2]
        )
      )
    )
  },
  .options = furrr_options(seed = TRUE)
) %>%
  # Convert the list of lists to a tibble.
  purrr::list_rbind() %>%
  # Add the results to the results_tbl.
  bind_cols(results_tbl, .)

# Drop variables that are not needed further on.
results_tbl = results_tbl %>%
  select(
    -time_points, 
    -K_vec
  )

saveRDS(
  object = results_tbl,
  file = here::here(dir_results, "results_simulation.rds")
)

print(Sys.time() - a)
