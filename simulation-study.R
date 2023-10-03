#!/usr/bin/env Rscript
.libPaths()
args = commandArgs(trailingOnly=TRUE)
print(args)
Sys.setenv(TZ='Europe/Brussels')
ncores = as.integer(args[1])
# Ensure to the state of packages is up-to-date.
renv::restore()
# Load the required packages. 
library(TCT)
library(mmrm)
library(dplyr)
library(doParallel)
a = Sys.time()
# Variance-Covariance matrix as reported by Raket (2022, doi: 10.1002/sim.9581) in
# section 5.1.
vcov_ref = matrix(c(45.1, 40.0, 45.1, 54.9, 53.6, 60.8,
                    40.0, 57.8, 54.4, 66.3, 64.1, 74.7,
                    45.1, 54.4, 72.0, 80.0, 77.6, 93.1, 
                    54.9, 66.3, 80.0, 109.8, 99.3, 121.7,
                    53.6, 64.1, 77.6, 99.3, 111.4, 127.8, 
                    60.8, 74.7, 93.1, 121.7, 127.8, 191.4), nrow = 6, byrow = TRUE)
# Mean values as reported by Raket (2022): `normal`. The `fast` vector
# represents an additional scenario where progression goes faster and/or the
# study duration is longer.
ref_means_list = list(
  normal = c(19.6, 20.5, 20.9, 22.7, 23.8, 27.4),
  fast = c(18, 19.7, 20.9, 22.7, 24.7, 29.2)
)

# In the settings data frame, we define all settings for which we will generate
# data. We consider all possible combinations of trial settings, except settings
# with extremely large or small power. 
settings = tidyr::expand_grid(
  progression = c("normal", "fast"),
  gamma_slowing = c(1, 0.75, 0.5),
  n = c(50, 200, 500, 1000),
  time_points = list(
    c(0, 6, 12, 18, 24),
    c(0, 6, 12, 18, 24, 36)
  )
) %>%
  # Exclude scenario with power almost equal to 1.
  filter(!(n == 1000 & gamma_slowing == 0.5)) %>%
  # Exclude scenario with very small power.
  filter(!(n == 50 & gamma_slowing == 0.75)) %>%
  # Exclude additional "less interesting" scenario to limit computational load.
  filter(!(n %in% c(200, 1000) & length(unlist(time_points) == 6)))

# Number of independent replications for each setting.
N_trials = 5
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
  # mmmrm::mmrm(). This causes a significant overhead. We therefore manually
  # remove the data saved into `fit`.
  mmrm_fit$data = NULL
  mmrm_fit$tmb_data = NULL
  mmrm_fit$tmb_object = NULL
  gc()
  return(mmrm_fit)
}

#---------


simulated_data_tbl = settings %>%
  # Group by each row.
  rowwise(everything()) %>%
  # Compute a set of variables that will be useful further on.
  summarise(
    # The number of measurement occasions.
    K = length(time_points),
    # Character vector for the naming the various measurement occasions.
    time_names = list(paste0("month_", time_points)),
    # Mean vector in the reference group.
    ref_means = list(ref_means_list[[progression]][1:K]),
    # Mean vector in the active treatment group.
    trt_means = list(
      spline(
        x = time_points,
        y = ref_means,
        xout = gamma_slowing * time_points
      )$y
    )
  ) %>%
  ungroup() %>%
  # For each row in the original tibble, the data will be generated for N_trial
  # trials.
  rowwise(c(
    "progression",
    "gamma_slowing",
    "n",
    "time_points",
    "K",
    "time_points"
  )) %>%
  reframe(
    tibble(as_tibble(
      # Generate data for the control group from a multivariate normal
      # distribution.
      mvtnorm::rmvnorm(
        n = (n / 2) * N_trials,
        mean = unlist(ref_means),
        sigma = vcov_ref[1:K, 1:K]
      ) %>% `colnames<-`(c(unlist(time_names)))
    ) ,
    arm = 0L) %>%
      bind_rows(tibble(
        # Generate data for the active treatment group.
        as_tibble(
          mvtnorm::rmvnorm(
            n = (n / 2) * N_trials,
            mean = unlist(trt_means),
            sigma = vcov_ref[1:K, 1:K]
          ) %>% `colnames<-`(c(unlist(time_names)))
        ),
        arm = 1L
      )) %>%
      mutate(
        # Add identifier for each generated trial.
        trial_number = rep(1:N_trials, times = n),
        SubjId = 1:(n * N_trials)
      ) %>%
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
        arm_time = ifelse(time_int == 0L,
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
                          ))
      )
  ) 

# The data are not saved very efficiently. For instance, variables indicating
# trial characteristics are duplicated for each outcome. Therefore, we construct
# a tibble where there is one row for each generated data set corresponding to a
# single trial.
simulated_data_tbl = simulated_data_tbl %>%
  group_by(progression, gamma_slowing, n, K, trial_number, time_points) %>%
  summarise(trial_data = list(pick(
    c("arm_time", "outcome", "time_int", "SubjId")
  ))) %>%
  ungroup()

print("Data Simulations Done")

# Fit the MMRM model for each generated data set.
cl = parallel::makeCluster(ncores)
results_tbl = simulated_data_tbl %>%
  mutate(mmrm_fit = parallel::parLapply(
    cl = cl,
    X = simulated_data_tbl$trial_data,
    fun = analyze_mmrm_new,
    type = "FULL"
  ))

print(Sys.time() - a)
print("MMRMs fitted")

# We now no longer need `simulated_data_tbl`. To free up space, the object
# containing the simulated data is removed.
rm("simulated_data_tbl")

# Compute p-value based on MMRM for each fitted model.
results_tbl = results_tbl %>%
  mutate(p_value_mmrm = sapply(
    X = mmrm_fit,
    FUN = function(x) {
      # Compute number of post-randomization measurement occasions.
      K = (length(coef(x)) - 1) / 2
      contrast = matrix(data = 0,
                        nrow = K,
                        ncol = length(coef(x)))
      contrast[1:K, 1:K] = diag(1, nrow = K)
      contrast[1:K, (K + 1):(2 * K)] = diag(-1, nrow = K)
      return(mmrm::df_md(x, contrast)$p_val)
    }
  ))

# Compute non-linear regression model.

# To save space, we now remove the generated data, and model fit objects. From
# the latter, we first extract the estimated variance-covariance matrix and
# estimated coefficients.
results_tbl = results_tbl %>%
  mutate(
    vcov_mmrm = lapply(
      X = mmrm_fit,
      FUN = vcov
    ),
    coef_mmrm = lapply(
      X = mmrm_fit,
      FUN = coef
    )
  ) %>%
  select(-trial_data, -mmrm_fit)

# Apply meta-TCT methodology.
results_tbl = results_tbl %>%
  # We first add additional columns that specify the inference options.
  cross_join(
    tidyr::expand_grid(
      drop_first_occasions = 0:2,
      constraints = c(TRUE, FALSE),
      interpolation = c("spline"),
      inference = c("wald", "score", "least-squares")
    ) %>%
      # We do not consider pre-selected measurement occasions if we adaptively
      # estimate the weights because this is contradictory. Either the weights are
      # selected a priori, or the weights are determined adaptively.
      filter(!(
        drop_first_occasions > 0 & inference == "score"
      ))
  ) %>%
  # Do not consider re-estimation under constraints for n = 1000 as violations
  # of the constraints are extremely unlikely for large sample sizes.
  filter(!(constraints & n == 1000)) %>%
  mutate(TCT_meta_fit = parallel::clusterMap(
    cl = cl,
    coef_mmrm = coef_mmrm,
    vcov_mmrm = vcov_mmrm,
    constraints = constraints,
    interpolation = interpolation,
    K = K, 
    time_points = time_points,
    fun = function(coef_mmrm,
                   vcov_mmrm,
                   constraints,
                   interpolation,
                   K,
                   time_points) {
      TCT::TCT_meta(
        time_points = time_points,
        exp_estimates = coef_mmrm[K:(2 * K - 2)],
        ctrl_estimates = coef_mmrm[c(2 * K - 1, 1:(K - 1))],
        vcov = vcov_mmrm,
        interpolation = interpolation,
        inference = "score",
        B = 0,
        constraints = constraints
      )
    }
  ))

attr(results_tbl$TCT_meta_fit, "split_type") = NULL
attr(results_tbl$TCT_meta_fit, "split_labels") = NULL
str(attributes(results_tbl$TCT_meta_fit))

print("meta-TCT finished")


# Estimate common acceleration factor.
results_tbl = results_tbl %>%
  mutate(TCT_meta_fit = parallel::clusterMap(
    cl = cl,
    TCT_meta_fit = TCT_meta_fit,
    drop_first_occasions = drop_first_occasions,
    inference = inference,
    constraints = constraints,
    fun = function(TCT_meta_fit,
                   drop_first_occasions,
                   inference,
                   constraints) {
      type = NULL
      if (inference == "score")
        type = "custom"
      TCT::TCT_meta_common(
        TCT_Fit = TCT_meta_fit,
        inference = inference,
        B = 0,
        select_coef = (drop_first_occasions + 1):length(coef(TCT_meta_fit)),
        constraints = constraints,
        type = type
      )
    }
  ))

# Compute confidence intervals and p-values. First, we need to call the summary
# method on the object returned by TCT_meta_common(). This method computes the
# p-value and confidence intervals. Because this can take some time (confidence
# intervals are computed numerically), we first save the summary-object to the
# tibble.
results_tbl = results_tbl %>%
  mutate(
    summary_TCT_common = parallel::parLapply(
      cl = cl,
      X = TCT_meta_fit,
      fun = summary
    )
  )

print("Common acceleration factors estimated")

#Compute the confidence intervals and p-values.
results_tbl = results_tbl %>%
  mutate(
    p_value_TCT_common = purrr::map_dbl(
      .x = summary_TCT_common, 
      .f = "p_value"
    ),
    conf_int_TCT_common = purrr::map(
      .x = summary_TCT_common,
      .f = "gamma_common_ci"
    ),
    se_TCT_common = purrr::map_dbl(
      .x = summary_TCT_common,
      .f = "gamma_common_se"
    )
  )

# Two versions of the simulation results are saved. First, a version containing
# all information, except the original simulated data. The latter cannot be
# included as the file size would be extremely large. This data set can be used
# to reproduce all results, starting from the fitted MMRM models. Second, a lean
# version of the simulation results that will be used for further processing. In
# the latter data set, only necessary information for processing is included.
results_lean_tbl = results_tbl %>%
  select(-vcov_mmrm, -coef_mmrm, -TCT_meta_fit, -summary_TCT_common)

saveRDS(object = results_tbl, file = "results_simulation_full.rds")
saveRDS(object = results_lean_tbl, file = "results_simulation_lean.rds")

print(Sys.time() - a)

