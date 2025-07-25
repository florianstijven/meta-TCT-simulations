# In this R-script, we make additional figures that illustrate the simulations.

# Preliminaries ----
library(tidyverse)

# Set directories where the plots and tables are saved.
dir_figures = here::here("results", "figures", "simulations")

theme_set(
  theme_get() +
    theme(
      legend.position = "bottom",
      legend.margin = margin(l = -1),
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.spacing.y = unit(0.1, "cm"),
      legend.box.spacing = unit(0.1, "cm")
    )
)

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
                xout = gamma_slowing * unlist(time_points), 
                method = "natural"
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
        xout = gamma_slowing * time_grid,
        method = "natural"
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
                     name = latex2exp::TeX("\\gamma"), 
                     direction = -1) +
  xlab("t") +
  ylab(expr(f[0](gamma %.% t))) +
  facet_grid(progression ~ .)
ggsave(filename = "dgm-mean-trajectories.pdf",
       path = dir_figures,
       device = "pdf",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Illustration of Misspecification with Correct Interpolation----

ref_means = c(19.6, 20.5, 20.9, 22.7, 23.8, 25.8, 27.4)
time_points_means_tbl = tibble(
  pattern = c("36 Months", "24 Months", "36(-30) Months"),
  time_points = list(
    c(0, 6, 12, 18, 24, 30, 36),
    c(0, 6, 12, 18, 24),
    c(0, 6, 12, 18, 24, 36)
  ),
  ref_means = list(
    ref_means,
    ref_means[-c(6, 7)],
    ref_means[-c(6)]
  )
)
time_grid = seq(from = 0,
                to = 36,
                length.out = 5e3)

trajectory_observed_tbl = time_points_means_tbl %>%
  cross_join(tibble(gamma = c(1, 0.9))) %>%
  rowwise(everything()) %>%
  summarise(trajectory_points = spline(
                x = time_points,
                y = ref_means,
                xout = gamma * unlist(time_grid), 
                method = "natural"
              )$y %>%
              list() 
            ) %>%
  ungroup() %>%
  select(-time_points, -ref_means) %>%
  rowwise(everything()) %>%
  reframe(tibble(
    trajectory_points = unlist(trajectory_points),
    time_grid = unlist(time_grid)
  ))

time_points_means_tbl_trt = time_points_means_tbl %>%
  cross_join(tibble(gamma = c(1, 0.9))) %>%
  rowwise(everything()) %>%
  summarise(
    trajectory_points = spline(
      x = time_points,
      y = ref_means,
      xout = gamma * unlist(time_points),
      method = "natural"
    )$y %>%
      list()
  ) %>%
  reframe(tibble(
    time_points = unlist(time_points),
    means = unlist(trajectory_points)
  ))

trajectory_observed_tbl %>%
  filter(pattern %in% c("36 Months", "36(-30) Months")) %>%
  ggplot(aes(x = time_grid, y = trajectory_points, linetype = pattern)) +
  geom_line(aes(color = as.factor(gamma))) +
  geom_point(
    data = time_points_means_tbl_trt %>%
      filter(pattern %in% c("36 Months")),
    aes(x = time_points, y = means)
  ) +
  scale_x_continuous(breaks = seq(0, 36, by = 6)) +
  scale_color_brewer(
    type = "qual",
    palette = 6,
    name = latex2exp::TeX("\\gamma"),
    direction = -1
  ) +
  # drop empty factor levels from legend
  scale_linetype(drop = TRUE) +
  xlab("t") +
  ylab(expr(f[0](gamma %.% t)))
ggsave(
  filename = "misspecification-mean-trajectories.pdf",
  path = dir_figures,
  device = "pdf",
  width = single_width,
  height = single_height,
  units = "cm",
  dpi = res
)


