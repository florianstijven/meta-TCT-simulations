#!/usr/bin/env Rscript
.libPaths()
args = commandArgs(trailingOnly=TRUE)
options(RENV_CONFIG_SANDBOX_ENABLED = FALSE)
Sys.setenv(TZ='Europe/Brussels')
ncores = as.integer(args[1])
# Print R-version for better reproducibility.
print(version)
# Ensure to the state of packages is up-to-date.
renv::restore()
# Load the required packages. 
library(TCT)
library(mmrm)
library(dplyr)
a = Sys.time()
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

# In the settings data frame, we define all settings for which we will generate
# data. We consider all possible combinations of trial settings. 
settings = tidyr::expand_grid(
  progression = c("normal", "fast"),
  gamma_slowing = c(1, 0.9, 0.75, 0.5),
  n = c(50, 200, 500, 1000),
  time_points_chr = c("24 Months", 
                      "36 Months", 
                      "36(-30) Months")
) 

settings = settings %>%
  left_join(tibble(
    time_points_chr = c("24 Months",
                        "36 Months",
                        "36(-30) Months"),
    time_points = list(c(0, 6, 12, 18, 24),
                       c(0, 6, 12, 18, 24, 30, 36),
                       c(0, 6, 12, 18, 24, 36))
  ))

# Number of independent replications for each setting.
N_trials = 5e3
# Set the seed for reproducibility.
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
                         vcov,
                         seed) {
  # Character vector for naming the various measurement occasions.
  time_names = paste0("month_", time_points)
  
  # Generate data.
  withr::local_seed(seed)
  
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
                                seed, 
                                reml = TRUE) {
  simulated_data = simulate_data(
    n = n,
    time_points = time_points,
    ref_means = ref_means,
    trt_means = trt_means,
    vcov = vcov,
    seed = seed
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
                               seed, 
                               reml = TRUE) {
  mmrm_fitted = simulate_and_analyze(n,
                       time_points,
                       ref_means,
                       trt_means,
                       vcov,
                       seed, 
                       reml = TRUE)
  extract_mmrm(mmrm_fitted)
}

#---------

settings = settings %>%
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
        xout = gamma_slowing * time_points
      )$y
    ),
    vcov = list(vcov_ref[K_vec, K_vec]),
  ) %>%
  ungroup()

# Duplicate rows N_trial times. We add a trial identifier which further serves
# as the seed to generate the trial data.
settings = settings %>%
  cross_join(tibble(trial_number = 1:N_trials))

# Fit the MMRM model for each generated data set.
cl = parallel::makeCluster(ncores)
parallel::clusterEvalQ(cl, library(TCT))
parallel::clusterEvalQ(cl, library(dplyr))
parallel::clusterExport(cl,
                        c(
                          "simulate_data",
                          "analyze_mmrm_new",
                          "simulate_and_analyze",
                          "extract_mmrm"
                        ))
results_tbl = settings %>%
  mutate(
    mmrm_extracted = parallel::clusterMap(
      cl = cl,
      n = n,
      time_points = time_points,
      ref_means = ref_means,
      trt_mean = trt_means,
      vcov = vcov,
      seed = trial_number,
      fun = simulate_to_extract
    )
  )

# results_tbl is not yet in an easy to process format because the mmrm_extracted
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
      constraints = c(TRUE, FALSE),
      interpolation = c("spline", "linear"),
      inference = c("wald", "score", "least-squares"),
      B = c(0, 5e2)
    ) %>%
      # We do not consider pre-selected measurement occasions if we adaptively
      # estimate the weights because this is contradictory. Either the weights are
      # selected a priori, or the weights are determined adaptively.
      filter(!(
        drop_first_occasions > 0 & inference == "score"
      ))
  ) %>%
  # We only use the bootstrap for least-squares with nore-estimation under 
  # constraints and no measurements dropped.
  filter(!(B == 5e2 &
             (constraints | drop_first_occasions == 1 | inference != "least-squares" | progression == "fast"))) %>%
  # If for a given setting, we consider bootstrap replications, then we can
  # delete the scenarios corresponding to B == 0 because we would be fitting the
  # same model twice.
  filter(!(B == 0 &
             !(constraints | drop_first_occasions == 1 | inference != "least-squares" | progression == "fast")))

results_tbl$TCT_meta_fit = parallel::clusterMap(
  cl = cl,
  coef_mmrm = results_tbl$coef_mmrm,
  vcov_mmrm = results_tbl$vcov_mmrm,
  constraints = results_tbl$constraints,
  interpolation = results_tbl$interpolation,
  K = results_tbl$K,
  time_points = results_tbl$time_points,
  fun = function(coef_mmrm,
                 vcov_mmrm,
                 constraints,
                 interpolation,
                 K,
                 time_points) {
    # The TCT_meta() function is called within a TryCatch() expression to
    # prevent the code from failing when one simulation fails. If TCT_meta()
    # fails, an NA is returned. 
    out = tryCatch(
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
    return(out)
  },
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE,
  RECYCLE = FALSE,
  .scheduling = "dynamic"
)

attr(results_tbl$TCT_meta_fit, "split_type") = NULL
attr(results_tbl$TCT_meta_fit, "split_labels") = NULL
str(attributes(results_tbl$TCT_meta_fit))

print(Sys.time() - a)
print("meta-TCT finished")


# Estimate common acceleration factor.
results_tbl$TCT_meta_common_fit = parallel::clusterMap(
  cl = cl,
  TCT_meta_fit = results_tbl$TCT_meta_fit,
  drop_first_occasions = results_tbl$drop_first_occasions,
  inference = results_tbl$inference,
  constraints = results_tbl$constraints,
  B = results_tbl$B,
  gamma_slowing = results_tbl$gamma_slowing,
  fun = function(TCT_meta_fit,
                 drop_first_occasions,
                 inference,
                 constraints,
                 B, 
                 gamma_slowing) {
    # Return NA if TCT_meta_fit is a missing value itself.
    if (!is.list(TCT_meta_fit)) return(NA)
    
    type = NULL
    if (inference == "score")
      type = "custom"
    # Some analyses use the bootstrap, and are thus stochastic. We therefore set
    # the seed next to avoid differences depending on the parallel setup.
    set.seed(1)
    out = tryCatch(
      TCT_meta_common(
        TCT_Fit = TCT_meta_fit,
        inference = inference,
        B = B,
        select_coef = (drop_first_occasions + 1):length(coef(TCT_meta_fit)),
        constraints = constraints,
        type = type,
        start_gamma = gamma_slowing
      ),
      error = function(cond) {
        return(NA)
      }
    )
    return(out)
  },
  SIMPLIFY = FALSE,
  USE.NAMES = TRUE,
  RECYCLE = FALSE,
  .scheduling = "dynamic"
)
print(Sys.time() - a)
print("meta-TCT-common finished")

# Compute confidence intervals and p-values. First, we need to call the summary
# method on the object returned by TCT_meta_common(). This method computes the
# p-value and confidence intervals. Because this can take some time (confidence
# intervals are computed numerically), we first save the summary-object to the
# tibble.
results_tbl$summary_TCT_common = parallel::parLapplyLB(
  cl = cl,
  X = results_tbl$TCT_meta_common_fit,
  fun = function(x) {
    if (is.list(x))
      return(summary(x))
    else
      return(NA)
  }
)
parallel::stopCluster(cl)

print("Common acceleration factors estimated")

# Compute the confidence intervals and p-values.
results_tbl = results_tbl %>%
  mutate(
    p_value_TCT_common = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = function(x) {
        if (is.list(x))
          return(x$p_value)
        else
          return(NA)
      }
    ),
    estimate = purrr::map_dbl(
      .x = TCT_meta_common_fit,
      .f = function(x) {
        if (is.list(x))
          return(coef(x))
        else
          return(NA)
      }
    ),
    conf_int_TCT_common_lower = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = function(x) {
        if (is.list(x))
          return(x$gamma_common_ci[1])
        else
          return(NA)
      }
    ),
    conf_int_TCT_common_upper = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = function(x) {
        if (is.list(x))
          return(x$gamma_common_ci[2])
        else
          return(NA)
      }
    ),
    se_TCT_common = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = function(x) {
        if (is.list(x))
          return(x$gamma_common_se)
        else
          return(NA)
      }
    ),
    se_TCT_common_bs = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = function(x) {
        if (!is.list(x))
          return(NA)
        else {
          if (is.null(x$se_bootstrap))
            return(NA)
          else
            return(x$se_bootstrap)
        }
        
      }
    ),
    conf_int_TCT_common_lower_bs = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = function(x) {
        if (!is.list(x))
          return(NA)
        else {
          if (is.null(x$se_bootstrap))
            return(NA)
          else
            return(x$ci_bootstrap[1])
        }
        
      }
    ),
    conf_int_TCT_common_upper_bs = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = function(x) {
        if (!is.list(x))
          return(NA)
        else{
          if (is.null(x$se_bootstrap))
            return(NA)
          else
            return(x$ci_bootstrap[2])
        }
        
      }
    )
  )

# Two versions of the simulation results are saved. First, a version containing
# all information, except the original simulated data. The latter cannot be
# included as the file size would be extremely large. This data set can be used
# to reproduce all results, starting from the fitted MMRM models. Second, a lean
# version of the simulation results that will be used for further processing. In
# the latter data set, only necessary information for processing is included.
results_lean_tbl = results_tbl %>%
  select(
    -vcov_mmrm,
    -coef_mmrm,
    -TCT_meta_fit,
    -summary_TCT_common,
    -TCT_meta_common_fit,
    -time_points, 
    -K_vec
  )

saveRDS(object = results_tbl, file = "results_simulation_full.rds")
saveRDS(object = results_lean_tbl, file = "results_simulation_lean.rds")
write.csv(x = results_lean_tbl,
          file = "results_simulation_lean.csv")

print(Sys.time() - a)

