#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
.libPaths(new = "/data/leuven/351/vsc35197/R")
print(args)
Sys.setenv(TZ='Europe/Brussels')
ncores = as.integer(args[1])
save_to = args[2]
# Load the required packages. 
library(TCT)
library(mmrm)
library(tidyverse)
library(nlme)
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
# data. We consider all possible combinations of trial settings. In total, there
# are 48 different settings. 
settings = tidyr::expand_grid(
  progression = c("normal", "fast"),
  gamma_slowing = c(1, 0.75, 0.5),
  n = c(50, 200, 500, 1000),
  time_points = list(
    c(0, 6, 12, 18, 24),
    c(0, 6, 12, 18, 24, 36)
  )
)
results = list()

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
  formula = formula(ADAScog_integer~arm_time + 0)
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

get_tct_est = function(data_trial_wide, index = 1:nrow(data_trial_wide)) {
  #convert data from wide to long format
  data_trial_long = data_trial_wide[index, ] %>%
    pivot_longer(cols = 4:8, names_to = "time_int", values_to = "ADAScog_integer") %>%
    mutate(SubjId = rep(1:nrow(data_trial_wide[index, ]), each = 5),
           time_int = as.integer(time_int),
           arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
  mmrm_fit = analyze_mmrm_new(data_trial_long, type = "full")
  TCT_fit = TCT(
    time_points = 0:4,
    ctrl_estimates = coef(mmrm_fit)[c(9, 1:4)],
    exp_estimates = coef(mmrm_fit)[5:8],
    vcov = vcov(mmrm_fit)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
    interpolation = "spline",
    B = 0
  )
  TCT_common_fit = TCT_common(
    TCT_Fit = TCT_fit,
    B = 0,
    bs_fix_vcov = FALSE,
    select_coef = 1:4
  )
  return(TCT_common_fit$coefficients)
}

bs_mmrm_tct = function(data_trial, R = 50) {
  data_trial_wide = data_trial %>%
    pivot_wider(id_cols = c("arm", "trial_number", "SubjId"), 
                names_from = "time_int", values_from = "ADAScog_integer")
  boot(data = data_trial_wide, 
       statistic = get_tct_est,
       R = R)
}

# Function to fit proportional slowing model directly. Code adapted from Raket
# (2022).

# Mean function based on cubic splines.
TPMRM <- function(t, v0, v1, v2, v3, v4,
                  log_gamma) {
  gamma = exp(log_gamma)
  months <- 0:4
  t_out <- gamma * t 
  
  spline(x = months,
         y = c(v0[1], v1[1], v2[1], v3[1], v4[1]),
         method = 'natural',
         xout = t_out)$y
}

proportional_slowing_nlm = function(data, start_vec) {
  # The first time point should be zero.
  t0 = min(data$time_int)
  data = data %>%
    mutate(time_int = time_int - t0)
    
  # Fit model
  pst_pmrm <-
    nlme::gnls(
      model = ADAScog_integer ~ TPMRM(time_int, v0, v1, v2, v3, v4,
                                      log_gamma),
      data = data,
      params = list(v0 + v1 + v2 + v3 + v4 ~ 1,
                    log_gamma ~ arm + 0),
      correlation = corSymm(form = ~ I(time_int + 1) | SubjId),
      weights = varIdent(form = ~ 1 | I(time_int + 1)),
      start = start_vec,
      control = gnlsControl(nlsTol = 1)
    )
  return(pst_pmrm)
}

#---------


settings_tbl = settings[1:8, ] %>%
  rowwise(everything()) %>%
  summarise(
    K = length(time_points),
    time_names = list(paste0("month_", time_points)),
    ref_means = list(ref_means_list[[progression]][1:K]),
    trt_means = list(
      spline(
        x = time_points,
        y = ref_means,
        xout = gamma_slowing * time_points
      )$y
    )
  ) %>%
  ungroup() %>%
  rowwise(c("progression", "gamma_slowing", "n", "time_points", "K")) %>%
  reframe(
    tibble(as_tibble(
      mvtnorm::rmvnorm(
        n = (n / 2) * N_trials,
        mean = unlist(ref_means),
        sigma = vcov_ref[1:K, 1:K]
      )
    ) %>% `colnames<-`(c(unlist(
      time_names
    ))),
    arm = 0L) %>%
      bind_rows(tibble(
        as_tibble(
          mvtnorm::rmvnorm(
            n = (n / 2) * N_trials,
            mean = unlist(trt_means),
            sigma = vcov_ref[1:K, 1:K]
          )
        ) %>% `colnames<-`(c(unlist(time_names))),
        arm = 1L
      )) %>%
      mutate(
        trial_number = rep(1:N_trials, times = n),
        SubjId = 1:(n * N_trials)
      ) %>%
      pivot_longer(cols = time_names, 
                   names_to = "time_arm",
                   values_to = "outcome") %>%
      mutate(time_int = rep(1:K, N_trials * n))
  ) 

settings_tbl = settings_tbl %>% 
  group_by(progression, gamma_slowing, n, K, trial_number) %>%
  summarise(trial_data = list(pick(c("time_arm", "outcome", "time_int", "SubjId"))))
  

# We iterate over all possible settings. 
for (i in 1:nrow(settings)) {
  # simulate data
  n = settings_tbl$n[i]
  data_trial = rbind(
    mvtnorm::rmvnorm(n = (n / 2) * N_trials, 
                     mean = time_means, 
                     sigma = vcov_ref),
    mvtnorm::rmvnorm(n = (n / 2) * N_trials, 
                     mean = time_means_trt, 
                     sigma = vcov_ref)
  )
  colnames(data_trial) = c("w0", "w25", "w50", "w75", "w100")
  data_trial = data.frame(data_trial)
  data_trial$arm = rep(0:1, each = (n / 2) * N_trials)
  data_trial$trial_number = rep(1:N_trials, times = n)
  data_trial$SubjId = 1:(n * N_trials)
  data_trial = data_trial %>%
    pivot_longer(cols = c("w0", "w25", "w50", "w75", "w100"), 
                 values_to = "ADAScog_integer")
  data_trial$time_int = rep(1:5, N_trials * n)
  data_trial = data_trial %>%
    mutate(arm_time = ifelse(time_int == 1L,
                             "baseline",
                             paste0(arm, ":", time_int)))
  
  list_of_trials = base::split(x = as.data.frame(data_trial), f = data_trial$trial_number)
  rm("data_trial")
  gc()
  print("simulated data")
  #------
  # MMRM
  cl = makeCluster(ncores)
  clusterEvalQ(cl, .libPaths(new = "/data/leuven/351/vsc35197/R"))
  print("sss")
  mmrm_fits = parallel::parLapply(cl = cl,
                                  X = list_of_trials,
                                  fun = analyze_mmrm_new,
                                  type = "full")
  mmrm_fits_null = parallel::parLapply(cl = cl,
                                  X = list_of_trials,
                                  fun = analyze_mmrm_new,
                                  type = "null", 
                                  reml = FALSE)
  stopCluster(cl)
  
  
  print("mmrm fitted")
  
  cl = makeCluster(ncores)
  clusterEvalQ(cl, .libPaths(new = "/data/leuven/351/vsc35197/R"))
  p1 = parSapply(
    X = mmrm_fits,
    FUN = function(x) {
      contrast = matrix(data = 0, nrow = 4, ncol = length(mmrm::component(x, "beta_est")))
      contrast[1:4, 1:4] = diag(1, nrow = 4)
      contrast[1:4, 5:8] = diag(-1, nrow = 4)
      return(mmrm::df_md(x, contrast)$p_val)
    },
    cl = cl
  )
  stopCluster(cl)
  print("mmrm p-values")

  
  # TCT
  cl = makeCluster(ncores)
  clusterEvalQ(cl, .libPaths(new = "/data/leuven/351/vsc35197/R"))
  clusterEvalQ(cl, library(TCT))
  clusterEvalQ(cl, library(mmrm))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(stats))
  clusterEvalQ(cl, library(mvtnorm))
  spline_common_1_4_results = parallel::parLapply(
    cl = cl,
    X = mmrm_fits,
    fun = function(x) {
      failed_return_value = list(
        p_bootstrap = NA,
        ci_bootstrap = NA,
        coefficients = NA,
        bootstrap_estimates = NA
      )
      error = TRUE
      try({
        TCT_fit = TCT(
          time_points = 0:4,
          ctrl_estimates = coef(x)[c(9, 1:4)],
          exp_estimates = coef(x)[5:8],
          vcov = vcov(x)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
          interpolation = "spline",
          B = 0
        )
        TCT_common_fit = TCT_common(
          TCT_Fit = TCT_fit,
          B = 2e3,
          bs_fix_vcov = FALSE,
          select_coef = 1:4
        )
        error = FALSE
      })
      if (error)
        return(failed_return_value)
      else
        return(summary(TCT_common_fit))
    }
  )
  spline_common_2_4_results = parallel::parLapply(
    cl = cl,
    X = mmrm_fits,
    fun = function(x) {
      failed_return_value = list(
        p_bootstrap = NA,
        ci_bootstrap = NA,
        coefficients = NA,
        bootstrap_estimates = NA
      )
      error = TRUE
      try({
        TCT_fit = TCT(
          time_points = 0:4,
          ctrl_estimates = coef(x)[c(9, 1:4)],
          exp_estimates = coef(x)[5:8],
          vcov = vcov(x)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
          interpolation = "spline",
          B = 0
        )
        TCT_common_fit = TCT_common(
          TCT_Fit = TCT_fit,
          B = 2e3,
          bs_fix_vcov = FALSE,
          select_coef = 2:4
        )
        error = FALSE
      })
      if (error)
        return(failed_return_value)
      else
        return(summary(TCT_common_fit))
    }
  )
  spline_common_3_4_results = parallel::parLapply(
    cl = cl,
    X = mmrm_fits,
    fun = function(x) {
      failed_return_value = list(
        p_bootstrap = NA,
        ci_bootstrap = NA,
        coefficients = NA,
        bootstrap_estimates = NA
      )
      error = TRUE
      try({
        TCT_fit = TCT(
          time_points = 0:4,
          ctrl_estimates = coef(x)[c(9, 1:4)],
          exp_estimates = coef(x)[5:8],
          vcov = vcov(x)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
          interpolation = "spline",
          B = 0
        )
        TCT_common_fit = TCT_common(
          TCT_Fit = TCT_fit,
          B = 2e3,
          bs_fix_vcov = FALSE,
          select_coef = 3:4
        )
        error = FALSE
      })
      if (error)
        return(failed_return_value)
      else
        return(summary(TCT_common_fit))
    }
  )
  stopCluster(cl)
  print("TCT MMRM")
  
  # TCT with constraints
  cl = makeCluster(ncores)
  clusterEvalQ(cl, .libPaths(new = "/data/leuven/351/vsc35197/R"))
  clusterEvalQ(cl, library(TCT))
  clusterEvalQ(cl, library(mmrm))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(stats))
  clusterEvalQ(cl, library(mvtnorm))
  spline_common_1_4_results_constraint = parallel::parLapply(
    cl = cl,
    X = mmrm_fits,
    fun = function(x) {
      failed_return_value = list(
        p_bootstrap = NA,
        ci_bootstrap = NA,
        coefficients = NA,
        bootstrap_estimates = NA
      )
      error = TRUE
      try({
        TCT_fit = TCT(
          time_points = 0:4,
          ctrl_estimates = coef(x)[c(9, 1:4)],
          exp_estimates = coef(x)[5:8],
          vcov = vcov(x)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
          interpolation = "spline",
          B = 0,
          constraints = TRUE
        )
        TCT_common_fit = TCT_common(
          TCT_Fit = TCT_fit,
          B = 2e3,
          bs_fix_vcov = FALSE,
          select_coef = 1:4
        )
        error = FALSE
      })
      if (error)
        return(failed_return_value)
      else
        return(summary(TCT_common_fit))
    }
  )
  spline_common_2_4_results_constraint = parallel::parLapply(
    cl = cl,
    X = mmrm_fits,
    fun = function(x) {
      failed_return_value = list(
        p_bootstrap = NA,
        ci_bootstrap = NA,
        coefficients = NA,
        bootstrap_estimates = NA
      )
      error = TRUE
      try({
        TCT_fit = TCT(
          time_points = 0:4,
          ctrl_estimates = coef(x)[c(9, 1:4)],
          exp_estimates = coef(x)[5:8],
          vcov = vcov(x)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
          interpolation = "spline",
          B = 0,
          constraints = TRUE
        )
        TCT_common_fit = TCT_common(
          TCT_Fit = TCT_fit,
          B = 2e3,
          bs_fix_vcov = FALSE,
          select_coef = 2:4
        )
        error = FALSE
      })
      if (error)
        return(failed_return_value)
      else
        return(summary(TCT_common_fit))
    }
  )
  spline_common_3_4_results_constraint = parallel::parLapply(
    cl = cl,
    X = mmrm_fits,
    fun = function(x) {
      failed_return_value = list(
        p_bootstrap = NA,
        ci_bootstrap = NA,
        coefficients = NA,
        bootstrap_estimates = NA
      )
      error = TRUE
      try({
        TCT_fit = TCT(
          time_points = 0:4,
          ctrl_estimates = coef(x)[c(9, 1:4)],
          exp_estimates = coef(x)[5:8],
          vcov = vcov(x)[c(9, 1:4, 5:8), c(9, 1:4, 5:8)],
          interpolation = "spline",
          B = 0,
          constraints = TRUE
        )
        TCT_common_fit = TCT_common(
          TCT_Fit = TCT_fit,
          B = 2e3,
          bs_fix_vcov = FALSE,
          select_coef = 3:4
        )
        error = FALSE
      })
      if (error)
        return(failed_return_value)
      else
        return(summary(TCT_common_fit))
    }
  )
  stopCluster(cl)
  print("TCT MMRM (constraints)")
  
  # Non-linear model with IPD.
  cl = parallel::makeCluster(ncores)
  clusterEvalQ(cl, .libPaths(new = "/data/leuven/351/vsc35197/R"))
  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(nlme))
  clusterExport(cl, c("proportional_slowing_nlm", "TPMRM"))

  ps_nlm_fits = parallel::parLapply(
    cl = cl,
    X = list_of_trials,
    fun = function(x, start_vec) {
      failed_return_value = list(
        tTable = matrix(NA, ncol = 5, nrow = 6)
      )
      error = TRUE
      try({
        nlm_fit = proportional_slowing_nlm(x, start_vec)
        error = FALSE
      })
      if (error)
        return(failed_return_value)
      else
        return(summary(nlm_fit))
    },
    start_vec = c(time_means, log(gamma_slowing_i))
  )
  stopCluster(cl)
  print("nlm fitted")


  results[[i]] = list(
    p1 = p1,
    nlm_min2loglik = purrr::map2_dbl(
      .x = mmrm_fits_null,
      .y = ps_nlm_fits,
      .f = function(.x, .y) {
        min2loglik = NA
        try({min2loglik = -2 * (logLik(.x) - logLik(.y))})
        return(min2loglik)
      }
    ), 
    mmrm_coef = purrr::map(mmrm_fits, stats::coef),
    nlm_gamma = purrr::map(ps_nlm_fits, function(x) return(x$tTable[6, ])),
    TCT_results = list(
      list(
        p_bs = purrr::map_dbl(spline_common_1_4_results, "p_bootstrap"),
        ci_bs = purrr::map(spline_common_1_4_results, "ci_bootstrap"),
        estimate = purrr::map_dbl(spline_common_1_4_results, "coefficients"),
        se_bs = purrr::map_dbl(spline_common_1_4_results,
                               function(x) {
                                 sd(x$bootstrap_estimates[[1]])
                               }),
        bs_estimates = purrr::map(spline_common_1_4_results, "bootstrap_estimates"),
        # bs_estimates_null = purrr::map(spline_common_1_4_results, "bootstrap_estimates_null"),
        type = paste(duration_i, as.character(gamma_slowing_i)),
        estimator = "1-4 unconstrained"
      ),
      list(
        p_bs = purrr::map_dbl(spline_common_2_4_results, "p_bootstrap"),
        ci_bs = purrr::map(spline_common_2_4_results, "ci_bootstrap"),
        estimate = purrr::map_dbl(spline_common_2_4_results, "coefficients"),
        se_bs = purrr::map_dbl(spline_common_2_4_results,
                               function(x) {
                                 sd(x$bootstrap_estimates[[1]])
                               }),
        bs_estimates = purrr::map(spline_common_2_4_results, "bootstrap_estimates"),
        # bs_estimates_null = purrr::map(spline_common_2_4_results, "bootstrap_estimates_null"),
        type = paste(duration_i, as.character(gamma_slowing_i)),
        estimator = "2-4 unconstrained"
      ),
      list(
        p_bs = purrr::map_dbl(spline_common_3_4_results, "p_bootstrap"),
        ci_bs = purrr::map(spline_common_3_4_results, "ci_bootstrap"),
        estimate = purrr::map_dbl(spline_common_3_4_results, "coefficients"),
        se_bs = purrr::map_dbl(spline_common_3_4_results,
                               function(x) {
                                 sd(x$bootstrap_estimates[[1]])
                               }),
        bs_estimates = purrr::map(spline_common_3_4_results, "bootstrap_estimates"),
        # bs_estimates_null = purrr::map(spline_common_3_4_results, "bootstrap_estimates_null"),
        type = paste(duration_i, as.character(gamma_slowing_i)),
        estimator = "3-4 unconstrained"
      ),
      list(
        p_bs = purrr::map_dbl(spline_common_1_4_results_constraint, "p_bootstrap"),
        ci_bs = purrr::map(spline_common_1_4_results_constraint, "ci_bootstrap"),
        estimate = purrr::map_dbl(spline_common_1_4_results_constraint, "coefficients"),
        se_bs = purrr::map_dbl(spline_common_1_4_results_constraint,
                               function(x) {
                                 sd(x$bootstrap_estimates[[1]])
                               }),
        bs_estimates = purrr::map(spline_common_1_4_results_constraint, "bootstrap_estimates"),
        # bs_estimates_null = purrr::map(spline_common_1_4_results_constraint, "bootstrap_estimates_null"),
        type = paste(duration_i, as.character(gamma_slowing_i)),
        estimator = "1-4 constrained"
      ),
      list(
        p_bs = purrr::map_dbl(spline_common_2_4_results_constraint, "p_bootstrap"),
        ci_bs = purrr::map(spline_common_2_4_results_constraint, "ci_bootstrap"),
        estimate = purrr::map_dbl(spline_common_2_4_results_constraint, "coefficients"),
        se_bs = purrr::map_dbl(spline_common_2_4_results_constraint,
                               function(x) {
                                 sd(x$bootstrap_estimates[[1]])
                               }),
        bs_estimates = purrr::map(spline_common_2_4_results_constraint, "bootstrap_estimates"),
        # bs_estimates_null = purrr::map(spline_common_2_4_results_constraint, "bootstrap_estimates_null"),
        type = paste(duration_i, as.character(gamma_slowing_i)),
        estimator = "2-4 constrained"
      ),
      list(
        p_bs = purrr::map_dbl(spline_common_3_4_results_constraint, "p_bootstrap"),
        ci_bs = purrr::map(spline_common_3_4_results_constraint, "ci_bootstrap"),
        estimate = purrr::map_dbl(spline_common_3_4_results_constraint, "coefficients"),
        se_bs = purrr::map_dbl(spline_common_3_4_results_constraint,
                               function(x) {
                                 sd(x$bootstrap_estimates[[1]])
                               }),
        bs_estimates = purrr::map(spline_common_3_4_results_constraint, "bootstrap_estimates"),
        # bs_estimates_null = purrr::map(spline_common_3_4_results_constraint, "bootstrap_estimates_null"),
        type = paste(duration_i, as.character(gamma_slowing_i)),
        estimator = "3-4 constrained"
      )
    )
  )
  
}

results_df = tibble(
  settings,
  results
)

print(Sys.time() - a)

file_path = paste0(save_to, "/results_df-08-03.rds")
saveRDS(object = results_df, file = file_path)
print(file_path)

