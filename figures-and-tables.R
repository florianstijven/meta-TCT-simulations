# In this R-script, the figures and tables, that are used in the paper, are
# produced.

# Preliminaries and Data Preparation ----
library(tidyverse)

# Size parameter for saving plots to disk.
single_width = 9
double_width = 14
single_height = 8.2
double_height = 12.8
res = 600

# Set the theme for all plots.
theme_set(theme(legend.position = "bottom",
                legend.margin = margin(l = -1), 
                legend.direction = "horizontal", 
                legend.box = "vertical"))


# Read data set with results of the simulation study.
results_tbl = readRDS(file = "results_simulation_lean.rds")

# Compute all summary measures regarding the performance of the estimators.
results_tbl_estimation = results_tbl %>%
  # Group by all simulation scenarios.
  group_by(
    progression,
    gamma_slowing,
    n,
    time_points_chr,
    drop_first_occasions,
    constraints,
    inference,
    interpolation
  ) %>%
  summarise(
    mean_estimate = mean(estimate),
    mean_estimate_se = sd(estimate) / sqrt(5e3),
    empirical_sd = sd(estimate),
    mean_se = mean(se_TCT_common),
    median_se = median(se_TCT_common),
    mean_se_bs = mean(se_TCT_common_bs),
    median_se_bs = median(se_TCT_common_bs),
    mse = mean((estimate - gamma_slowing) ** 2)
  ) %>%
  ungroup() %>%
  rename("Follow Up" = time_points_chr) %>%
  mutate(
    interpolation = forcats::fct_recode(
      interpolation,
      "Linear" = "linear",
      "Natural Cubic Spline" = "spline"
    ),
    progression = forcats::fct_recode(
      progression,
      "Fast" = "fast",
      "Normal" = "normal"
    )
  )

results_tbl_inference = results_tbl %>%
  mutate(
    # Determine whether the CI contain the true value.
    covered = (gamma_slowing <= conf_int_TCT_common_upper) &
      (gamma_slowing >= conf_int_TCT_common_lower),
    covered_bs = (gamma_slowing <= conf_int_TCT_common_upper_bs) &
      (gamma_slowing >= conf_int_TCT_common_lower_bs),
    # Determine the conclusion based on a bootstrap-based CI test.
    rejection_bs = !((1 <= conf_int_TCT_common_upper_bs) &
                       (1 >= conf_int_TCT_common_lower_bs)
    ),
    # Determine whether the null has been rejected
    rejection = p_value_TCT_common <= 0.05,
    rejection_mmrm = p_value_mmrm <= 0.05
  ) %>%
  group_by(progression,
           gamma_slowing,
           n,
           time_points_chr,
           drop_first_occasions,
           constraints,
           inference,
           interpolation) %>%
  summarise(
    coverage = mean(covered, na.rm = TRUE),
    prop_failed_ci = mean(is.na(covered)),
    coverage_bs = mean(covered_bs),
    rejection_rate_bs = mean(rejection_bs),
    rejection_rate = mean(rejection),
    rejection_rate_mmrm = mean(rejection_mmrm)
  ) %>%
  ungroup() %>%
  rename("Follow Up" = time_points_chr) %>%
  mutate(
    interpolation = forcats::fct_recode(
      interpolation,
      "Linear" = "linear",
      "Natural Cubic Spline" = "spline"
    ),
    progression = forcats::fct_recode(
      progression,
      "Fast" = "fast",
      "Normal" = "normal"
    )
  )



# Simulation Results with Inference Based on the Least-Squares Criterion ----

## Properties of the Estimator ----

# Define dummy tibble that allows us to set the axis limits in each facet
# separately.
facet_lims_tbl = tibble(
  n = NA,
  gamma_slowing = rep(c(0.5, 0.75, 0.9, 1), each = 2),
  mean_estimate = c(0.40, 0.60, 0.65, 0.85, 0.80, 1.00, 0.90, 1.10),
  rejection_rate = c(0, 1, 0, 1, 0, 1, 0, 0.25),
  `Follow Up` = c("24 Months"),
  progression = NA,
  interpolation = NA
)
results_tbl_estimation %>%
  filter(constraints == FALSE,
         inference == "least-squares",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = mean_estimate,
    color = interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_hline(mapping = aes(yintercept = gamma_slowing),
             alpha = 0.5) +
  geom_hline(mapping = aes(yintercept = gamma_slowing + 0.05),
             alpha = 0.30) +
  geom_hline(mapping = aes(yintercept = gamma_slowing - 0.05),
             alpha = 0.30) +
  geom_blank(data = facet_lims_tbl, aes(gamma_slowing, mean_estimate)) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("E(\\hat{\\gamma})")) +
  scale_linetype(name = "Progression Rate") +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months")),
    scales = "free"
  ) +
  guides(color = guide_legend(label.position = "right")) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")  +
  scale_y_continuous(n.breaks = 3)
ggsave(filename = "figures/main-text/bias-nl-gls.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_estimation %>%
  filter(constraints == FALSE,
         inference == "least-squares",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = mse,
    color = `interpolation`,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  scale_y_continuous(trans = "log10", name = "MSE", limits = c(0.0005, 5)) +
  scale_x_continuous(trans = "log10", name = "Sample Size (n)") +
  scale_linetype(name = "Progression Rate") +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/web-appendix/mse-nl-gls.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_estimation %>%
  filter(constraints == FALSE,
         inference == "least-squares",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = empirical_sd,
    color = interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_point(
    data = results_tbl_estimation %>%
      filter(constraints == FALSE,
             inference == "least-squares",
             drop_first_occasions == 0),
    mapping = aes(x = n, y = median_se),
    shape = 2
  ) +
  scale_y_continuous(trans = "log10", 
                     name = "Empirical SD and Median Estimated SE",
                     limits = c(0.025, 1.9)) +
  scale_x_continuous(trans = "log10", name = "Sample Size (n)") +
  scale_linetype(name = "Progression Rate") +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/main-text/se-nl-gls.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

## Inferential Properties ----

results_tbl_inference %>%
  filter(constraints == FALSE,
         inference == "least-squares",
         drop_first_occasions == 0) %>%
  ggplot(
    aes(
      x = n,
      y = rejection_rate,
      color = interpolation,
      linetype = progression
    )
  ) +
  geom_point() +
  geom_line() +
  geom_point(
    data = results_tbl_inference %>%
      filter(constraints == FALSE,
             inference == "least-squares",
             drop_first_occasions == 0),
    mapping = aes(x = n,
                  y = rejection_rate_mmrm),
    alpha = 0.5,
    color = "gray"
  ) +
  geom_line(
    data = results_tbl_inference %>%
      filter(constraints == FALSE,
             inference == "least-squares",
             drop_first_occasions == 0),
    mapping = aes(x = n,
                  y = rejection_rate_mmrm,
                  linetype = progression),
    alpha = 0.75,
    color = "gray"
  ) +
  geom_hline(
    # Add horizontal lines to indicate the 0.05 nominal error rate in the null
    # settings.
    data = tidyr::expand_grid(
      rejection_rate = 0.05,
      drop_first_occasions = 0,
      gamma_slowing = 1
    ),
    mapping = aes(yintercept = rejection_rate),
    alpha = 0.5
  ) +
  geom_blank(data = facet_lims_tbl, aes(gamma_slowing, rejection_rate)) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("Empirical Power or Type 1 Error Rate")) +
  scale_linetype(name = "Progression Rate") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months")), scales = "free") +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/web-appendix/power-type1-error-nl-gls.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_inference %>%
  filter(constraints == FALSE,
         inference == "least-squares",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = coverage,
    color = interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_hline(
    # Add horizontal lines to indicate the 0.95 nominal coverage.
    data = tidyr::expand_grid(
      coverage = 0.95,
      drop_first_occasions = 0:1,
      gamma_slowing = c(0.5, 0.75, 0.9, 1)
    ),
    mapping = aes(yintercept = coverage),
    alpha = 0.5
  ) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("Empirical Coverage")) +
  ylim(c(0.75, 1)) +
  scale_linetype(name = "Progression Rate") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/main-text/coverage-nl-gls.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

## Simulation Results with Parametric Bootstrap ----

results_tbl_estimation %>%
  # Change progression from character to factor variable. This ensure that the
  # line types are consistent with the previous graphs.
  mutate(progression = factor(progression)) %>%
  filter(
    constraints == FALSE,
    inference == "least-squares",
    progression == "Normal",
    drop_first_occasions == 0
  ) %>%
  ggplot(aes(
    x = n,
    y = empirical_sd,
    color =interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_point(
    data = results_tbl_estimation %>%
      filter(
        constraints == FALSE,
        inference == "least-squares",
        progression == "Normal",
        drop_first_occasions == 0
      ),
    mapping = aes(x = n, y = median_se_bs),
    shape = 2
  ) +
  scale_y_continuous(trans = "log10", name = "Empirical SD and Median Estimated SE") +
  scale_x_continuous(trans = "log10", name = "Sample Size (n)") +
  scale_linetype(name = "Progression Rate", drop = FALSE) +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/web-appendix/se-nl-gls-bs.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_inference %>%
  # Change progression from character to factor variable. This ensure that the
  # line types are consistent with the previous graphs.
  mutate(progression = as.factor(progression)) %>%
  filter(
    constraints == FALSE,
    inference == "least-squares",
    progression == "Normal",
    drop_first_occasions == 0
  ) %>%
  ggplot(aes(
    x = n,
    y = coverage_bs,
    color = interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_hline(
    # Add horizontal lines to indicate the 0.95 nominal coverage.
    data = tidyr::expand_grid(
      coverage = 0.95,
      gamma_slowing = c(0.5, 0.75, 0.9, 1)
    ),
    mapping = aes(yintercept = coverage),
    alpha = 0.5
  ) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("Empirical Coverage")) +
  ylim(c(0.75, 1)) +
  scale_linetype(name = "Progression Rate", drop = FALSE) +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(`Follow Up`,
                                    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/main-text/coverage-nl-gls-bs.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_inference %>%
  # Change progression from character to factor variable. This ensure that the
  # line types are consistent with the previous graphs.
  mutate(progression = as.factor(progression)) %>%
  filter(
    constraints == FALSE,
    inference == "least-squares",
    drop_first_occasions == 0,
    progression == "Normal"
  ) %>%
  ggplot(
    aes(
      x = n,
      y = rejection_rate_bs,
      color = interpolation,
      linetype = progression
    )
  ) +
  geom_point() +
  geom_line() +
  geom_point(
    data = results_tbl_inference %>%
      filter(
        constraints == FALSE,
        inference == "least-squares",
        drop_first_occasions == 0,
        progression == "Normal"
      ),
    mapping = aes(x = n,
                  y = rejection_rate_mmrm),
    alpha = 0.5,
    color = "gray"
  ) +
  geom_line(
    data = results_tbl_inference %>%
      filter(
        constraints == FALSE,
        inference == "least-squares",
        drop_first_occasions == 0,
        progression == "Normal"
      ),
    mapping = aes(x = n,
                  y = rejection_rate_mmrm,
                  linetype = progression),
    alpha = 0.75,
    color = "gray"
  ) +
  geom_hline(
    # Add horizontal lines to indicate the 0.05 nominal error rate in the null
    # settings.
    data = tidyr::expand_grid(
      rejection_rate = 0.05,
      drop_first_occasions = 0,
      gamma_slowing = 1
    ),
    mapping = aes(yintercept = rejection_rate),
    alpha = 0.5
  ) +
  geom_blank(data = facet_lims_tbl, aes(gamma_slowing, rejection_rate)) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("Empirical Power or Type 1 Error Rate")) +
  scale_linetype(name = "Progression Rate", drop = FALSE) +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(`Follow Up`,
                                    levels = c("24 Months", "36(-30) Months", "36 Months")), scales = "free") +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/web-appendix/error-rates-nl-gls-bs.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)



# Simulation Results with the adaptive score-based estimator ----

## Properties of the Estimator ----

# Define dummy tibble that allows us to set the axis limits in each facet
# separately.
results_tbl_estimation %>%
  filter(constraints == FALSE,
         inference == "score",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = mean_estimate,
    color = interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_hline(mapping = aes(yintercept = gamma_slowing),
             alpha = 0.5) +
  geom_hline(mapping = aes(yintercept = gamma_slowing + 0.05),
             alpha = 0.30) +
  geom_hline(mapping = aes(yintercept = gamma_slowing - 0.05),
             alpha = 0.30) +
  geom_blank(data = facet_lims_tbl, aes(gamma_slowing, mean_estimate)) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("E(\\hat{\\gamma})")) +
  scale_linetype(name = "Progression Rate") +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months")),
    scales = "free"
  ) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation") +
  scale_y_continuous(n.breaks = 3)
ggsave(filename = "figures/web-appendix/bias-score-adap.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_estimation %>%
  filter(constraints == FALSE,
         inference == "score",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = mse,
    color = `interpolation`,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  scale_y_continuous(trans = "log10", name = "MSE", limits = c(0.0005, 5)) +
  scale_x_continuous(trans = "log10", name = "Sample Size (n)") +
  scale_linetype(name = "Progression Rate") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/web-appendix/mse-score-adap.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_estimation %>%
  filter(constraints == FALSE,
         inference == "score",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = empirical_sd,
    color = interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_point(
    data = results_tbl_estimation %>%
      filter(constraints == FALSE,
             inference == "score",
             drop_first_occasions == 0),
    mapping = aes(x = n, y = median_se),
    shape = 2
  ) +
  scale_y_continuous(trans = "log10", 
                     name = "Empirical SD and Median Estimated SE",
                     limits = c(0.025, 1.9)) +
  scale_x_continuous(trans = "log10", name = "Sample Size (n)") +
  scale_linetype(name = "Progression Rate") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/web-appendix/se-score-adap.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

## Inferential Properties ----

results_tbl_inference %>%
  filter(constraints == FALSE,
         inference == "score",
         drop_first_occasions == 0) %>%
  ggplot(
    aes(
      x = n,
      y = rejection_rate,
      color = interpolation,
      linetype = progression
    )
  ) +
  geom_point() +
  geom_line() +
  geom_point(
    data = results_tbl_inference %>%
      filter(constraints == FALSE,
             inference == "score",
             drop_first_occasions == 0),
    mapping = aes(x = n,
                  y = rejection_rate_mmrm),
    alpha = 0.5,
    color = "gray"
  ) +
  geom_line(
    data = results_tbl_inference %>%
      filter(constraints == FALSE,
             inference == "score",
             drop_first_occasions == 0),
    mapping = aes(x = n,
                  y = rejection_rate_mmrm,
                  linetype = progression),
    alpha = 0.75,
    color = "gray"
  ) +
  geom_hline(
    # Add horizontal lines to indicate the 0.05 nominal error rate in the null
    # settings.
    data = tidyr::expand_grid(
      rejection_rate = 0.05,
      drop_first_occasions = 0,
      gamma_slowing = 1
    ),
    mapping = aes(yintercept = rejection_rate),
    alpha = 0.5
  ) +
  geom_blank(data = facet_lims_tbl, aes(gamma_slowing, rejection_rate)) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("Empirical Power or Type 1 Error Rate")) +
  scale_linetype(name = "Progression Rate") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months")), scales = "free") +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
  ggsave(filename = "figures/web-appendix/power-type1-error-score-adap.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

results_tbl_inference %>%
  filter(constraints == FALSE,
         inference == "score",
         drop_first_occasions == 0) %>%
  ggplot(aes(
    x = n,
    y = coverage,
    color = interpolation,
    linetype = progression
  )) +
  geom_point() +
  geom_line() +
  geom_hline(
    # Add horizontal lines to indicate the 0.95 nominal coverage.
    data = tidyr::expand_grid(
      coverage = 0.95,
      drop_first_occasions = 0:1,
      gamma_slowing = c(0.5, 0.75, 0.9, 1)
    ),
    mapping = aes(yintercept = coverage),
    alpha = 0.5
  ) +
  xlab("Sample Size (n)") +
  ylab(latex2exp::TeX("Empirical Coverage")) +
  ylim(c(0.75, 1)) +
  scale_linetype(name = "Progression Rate") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid(gamma_slowing ~ factor(
    `Follow Up`,
    levels = c("24 Months", "36(-30) Months", "36 Months"))) +
  scale_color_brewer(type = "qual", palette = 2, name = "Interpolation")
ggsave(filename = "figures/web-appendix/coverage-score-adap.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

# For the previous plot, there were a few settings where the confidence
# intervals could not always be computed. However, this proportion is maximum
# 0.2%. This also only happened in the n = 50 settings.

# Data Generating Mechanism ----

ref_means_list = list(
  "Normal Progression" = c(19.6, 20.5, 20.9, 22.7, 23.8, 25.8, 27.4),
  "Fast Progression" = c(18, 19.7, 20.9, 22.7, 24.7, 27.1, 29.2)
)
time_points = c(0, 6, 12, 18, 24, 30, 36)
time_grid = seq(from = 0,
                to = 36,
                length.out = 3000)
dgm_settings = tidyr::expand_grid(
  progression = c("Normal Progression", "Fast Progression"),
  gamma_slowing = c(1, 0.9, 0.75, 0.5),
  time_points = list(time_points)
)
trajectory_observed_tbl = dgm_settings %>%
  rowwise(everything()) %>%
  summarise(ref_means = ref_means_list[progression],
            trajectory_points = list(
              spline(
                x = time_points,
                y = ref_means,
                xout = gamma_slowing * unlist(time_points)
              )$y
            )) %>%
  ungroup() %>%
  mutate(treatment = ifelse(gamma_slowing == 1, "Control Treatment", "Active Treatment"))
trajectory_observed_tbl = trajectory_observed_tbl %>%
  rowwise(everything()) %>%
  reframe(tibble(
    trajectory_points = unlist(trajectory_points),
    time_points = unlist(time_points)
  )) %>%
  ungroup() %>%
  mutate(gamma_slowing = as.factor(gamma_slowing))

trajectory_interpolated_tbl = dgm_settings %>%
  rowwise(everything()) %>%
  summarise(
    ref_means = ref_means_list[progression],
    trajectory_interpolated = list(
      spline(
        x = time_points,
        y = ref_means,
        xout = gamma_slowing * time_grid
      )$y
    ),
    time_grid = list(time_grid)
  ) %>%
  ungroup() %>%
  mutate(treatment = ifelse(gamma_slowing == 1, "Control Treatment", "Active Treatment"))
trajectory_interpolated_tbl = trajectory_interpolated_tbl %>%
  rowwise(everything()) %>%
  reframe(tibble(
    trajectory_interpolated = unlist(trajectory_interpolated),
    time_grid = unlist(time_grid)
  )) %>%
  ungroup() %>%
  mutate(gamma_slowing = as.factor(gamma_slowing))

trajectory_observed_tbl %>%
  ggplot(aes(x = time_points, y = trajectory_points)) +
  geom_point(alpha = 0.25) +
  geom_line(
    data = trajectory_interpolated_tbl,
    mapping = aes(
      x = time_grid,
      y = trajectory_interpolated,
      group = interaction(gamma_slowing, treatment, progression),
      color = gamma_slowing
    )
  ) +
  scale_x_continuous(breaks = time_points) +
  scale_color_brewer(type = "qual", 
                     palette = 6,
                     name = latex2exp::TeX("Acceleration Factor")) +
  xlab("Months Since Randomization") +
  ylab("Mean Trajectory") +
  theme(legend.position = "bottom") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1)) +
  facet_grid( ~ progression)
ggsave(filename = "figures/main-text/dgm-mean-trajectories.png",
       device = "png",
       width = double_width,
       height = double_height,
       units = "cm",
       dpi = res)

# Consider the "true" trajectories for the problematic setting
trajectory_interpolated_problematic_tbl = tibble(
  progression = c("Normal Progression"),
  gamma_slowing = c(0.9, 1),
  time_points = list(time_points[-6])
) %>%
  rowwise(everything()) %>%
  summarise(
    ref_means = ref_means_list[progression],
    trajectory_interpolated = list(
      spline(
        x = time_points,
        y = ref_means[-6],
        xout = gamma_slowing * time_grid
      )$y
    ),
    time_grid = list(time_grid)
  ) %>%
  ungroup() %>%
  mutate(treatment = ifelse(gamma_slowing == 1, "Control Treatment", "Active Treatment"))

trajectory_interpolated_problematic_tbl = trajectory_interpolated_problematic_tbl %>%
  rowwise(everything()) %>%
  reframe(tibble(
    trajectory_interpolated = unlist(trajectory_interpolated),
    time_grid = unlist(time_grid)
  )) %>%
  ungroup() %>%
  mutate(gamma_slowing = as.factor(gamma_slowing))

bind_rows(trajectory_interpolated_problematic_tbl %>%
            mutate("Measurement Pattern" = "36(-30) Months"),
          trajectory_interpolated_tbl %>%
            filter(progression == "Normal Progression",
                   gamma_slowing %in% c(0.9, 1)) %>%
            mutate("Measurement Pattern" = "36 Months")
          ) %>%
ggplot(aes(
  x = time_grid,
  y = trajectory_interpolated,
  color = gamma_slowing,
  linetype = `Measurement Pattern`
)) +
  geom_line() +
  scale_x_continuous(breaks = time_points) +
  scale_color_brewer(type = "qual", 
                     palette = 6,
                     name = latex2exp::TeX("Acceleration Factor")) +
  xlab("Months Since Randomization") +
  ylab("Interpolated Mean Trajectory") +
  theme(legend.position = "bottom") +
  theme(legend.position = "bottom", legend.margin = margin(l = -1), 
        legend.direction = "vertical", legend.box = "horizontal") +
  facet_grid( ~ progression)
ggsave(filename = "figures/main-text/problematic-mean-trajectories.png",
       device = "png",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)
