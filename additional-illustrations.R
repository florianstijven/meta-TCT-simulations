library(tidyverse)

save_to = "figures/additional-illustrations/"
# Time Mapping ------------------------------------------------------------

ref_means = c(17.5, 20.2, 24.7, 29.2)
trt_means = c(18.5, 23, 27)
time_points = c(0, 10, 22, 36)
time_grid = seq(from = 0,
                to = 36,
                length.out = 3000)
trajectory_points = spline(
  x = time_points,
  y = ref_means,
  xout = time_grid,
  method = "natural"
)$y
trt_mapped_times = TCT:::get_new_time(
  y_ref = ref_means,
  x_ref = time_points,
  y_obs = trt_means,
  method = "spline"
)

# Construct the base plot for illustrating the time mappings. In this plot, we
# add the control and active group means, and the reference trajectory (but now
# arrows).
time_mapping_base_plot = tibble(time_grid, trajectory_points) %>%
  ggplot(aes(x = time_grid, y = trajectory_points)) +
  geom_line(color = RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]) +
  geom_point(
    data = tibble(time_points, ref_means),
    aes(x = time_points, y = ref_means),
    alpha = 0.25
  ) +
  geom_point(data = tibble(time_points = time_points[-1], trt_means),
             aes(x = time_points, y = trt_means)) +
  scale_x_continuous(
    name = "Time Since Randomization",
    breaks = c(time_points, trt_mapped_times),
    labels = c(
      "0",
      latex2exp::TeX("$t_1$"),
      latex2exp::TeX("$t_2$"),
      latex2exp::TeX("$t_3$"),
      latex2exp::TeX("$\\gamma_1 t_1$"),
      latex2exp::TeX("$\\gamma_2 t_2$"),
      latex2exp::TeX("$\\gamma_3 t_3$")
    ),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    name = "Mean Outcome",
    breaks = c(ref_means, trt_means),
    labels = c(
      latex2exp::TeX("$\\alpha_0$"),
      latex2exp::TeX("$\\alpha_1$"),
      latex2exp::TeX("$\\alpha_2$"),
      latex2exp::TeX("$\\alpha_3$"),
      latex2exp::TeX("$\\beta_1$"),
      latex2exp::TeX("$\\beta_2$"),
      latex2exp::TeX("$\\beta_3$")
    ),
    minor_breaks = NULL
  )

# Add vertical arrows to the base plot. This illustrates classic comparisons in
# mixed models.
time_mapping_base_plot + geom_segment(
  data = tibble(
    time_points = time_points[-1],
    trt_mapped_times,
    ref_means = ref_means[-1],
    trt_means
  ),
  aes(
    x = time_points,
    xend = time_points,
    y = trt_means,
    yend = ref_means
  ),
  arrow = arrow(length = unit(0.25, 'cm'))
) 

ggsave(filename = paste0(save_to, "time-mapping-vertical.pdf"),
       device = "pdf",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Add vertical arrows to the base plot. This illustrates classic comparisons in
# mixed models.
time_mapping_base_plot + geom_segment(
  data = tibble(
    time_points = time_points[-1],
    trt_mapped_times,
    ref_means = ref_means[-1],
    trt_means
  ),
  aes(
    x = time_points,
    xend = trt_mapped_times,
    y = trt_means,
    yend = trt_means
  ),
  arrow = arrow(length = unit(0.25, 'cm'))
) 

ggsave(filename = paste0(save_to, "time-mapping-horizontal.pdf"),
       device = "pdf",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

# Time-Mapping Function ---------------------------------------------------

## Proportional Slowing ---------------------------------------------------

ref_means = c(17.5, 20.2, 24.7, 29.2)
trt_means = c(18.5, 23, 27)
time_points = c(0, 10, 22, 36)
time_grid = seq(from = 0,
                to = 36,
                length.out = 3000)
trajectory_points = spline(
  x = time_points,
  y = ref_means,
  xout = time_grid,
  method = "natural"
)$y
trt_mapped_times = TCT:::get_new_time(
  y_ref = ref_means,
  x_ref = time_points,
  y_obs = trt_means,
  method = "spline"
)

tibble(time_grid, trajectory_points) %>%
  ggplot(aes(x = time_grid, y = trajectory_points)) +
  geom_line(color = RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]) +
  geom_point(
    data = tibble(time_points, ref_means),
    aes(x = time_points, y = ref_means),
    alpha = 0.25
  ) +
  geom_point(data = tibble(time_points = time_points[-1], trt_means),
             aes(x = time_points, y = trt_means)) +
  geom_segment(
    data = tibble(
      time_points = time_points[-1],
      trt_mapped_times,
      ref_means = ref_means[-1],
      trt_means
    ),
    aes(
      x = time_points,
      xend = trt_mapped_times,
      y = trt_means,
      yend = trt_means
    ),
    arrow = arrow(length = unit(0.25, 'cm'))
  ) +
  scale_x_continuous(
    name = "Time Since Randomization",
    breaks = c(time_points, trt_mapped_times),
    labels = c(
      "0",
      latex2exp::TeX("$t_1$"),
      latex2exp::TeX("$t_2$"),
      latex2exp::TeX("$t_3$"),
      latex2exp::TeX("$\\gamma_1 t_1$"),
      latex2exp::TeX("$\\gamma_2 t_2$"),
      latex2exp::TeX("$\\gamma_3 t_3$")
    ),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    name = "Mean Outcome",
    breaks = c(ref_means, trt_means),
    labels = c(
      latex2exp::TeX("$\\alpha_0$"),
      latex2exp::TeX("$\\alpha_1$"),
      latex2exp::TeX("$\\alpha_2$"),
      latex2exp::TeX("$\\alpha_3$"),
      latex2exp::TeX("$\\beta_1$"),
      latex2exp::TeX("$\\beta_2$"),
      latex2exp::TeX("$\\beta_3$")
    ),
    minor_breaks = NULL
  )

ggsave(filename = "figures/additional-illustrations/proportional-slowing.pdf",
       device = "pdf",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

## Proportional Slowing with an Offset ----------------------------------

ref_means = c(17.5, 20.2, 24.7, 29.2)
trt_means = c(18.5, 23, 27)
time_points = c(0, 10, 22, 36)
time_grid = seq(from = 0,
                to = 36,
                length.out = 3000)
trajectory_points = spline(
  x = time_points,
  y = ref_means,
  xout = time_grid,
  method = "natural"
)$y
trt_mapped_times = TCT:::get_new_time(
  y_ref = ref_means,
  x_ref = time_points,
  y_obs = trt_means,
  method = "spline"
)

tibble(time_grid, trajectory_points) %>%
  ggplot(aes(x = time_grid, y = trajectory_points)) +
  geom_line(color = RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]) +
  geom_point(
    data = tibble(time_points, ref_means),
    aes(x = time_points, y = ref_means),
    alpha = 0.25
  ) +
  geom_point(data = tibble(time_points = time_points[-1], trt_means),
             aes(x = time_points, y = trt_means)) +
  geom_segment(
    data = tibble(
      time_points = time_points[-1],
      trt_mapped_times,
      ref_means = ref_means[-1],
      trt_means
    ),
    aes(
      x = time_points,
      xend = trt_mapped_times,
      y = trt_means,
      yend = trt_means
    ),
    arrow = arrow(length = unit(0.25, 'cm'))
  ) +
  scale_x_continuous(
    name = "Time Since Randomization",
    breaks = c(time_points, trt_mapped_times),
    labels = c(
      "0",
      latex2exp::TeX("$t_1$"),
      latex2exp::TeX("$t_2$"),
      latex2exp::TeX("$t_3$"),
      latex2exp::TeX("$\\gamma_1 t_1$"),
      latex2exp::TeX("$\\gamma_2 t_2$"),
      latex2exp::TeX("$\\gamma_3 t_3$")
    ),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    name = "Mean Outcome",
    breaks = c(ref_means, trt_means),
    labels = c(
      latex2exp::TeX("$\\alpha_0$"),
      latex2exp::TeX("$\\alpha_1$"),
      latex2exp::TeX("$\\alpha_2$"),
      latex2exp::TeX("$\\alpha_3$"),
      latex2exp::TeX("$\\beta_1$"),
      latex2exp::TeX("$\\beta_2$"),
      latex2exp::TeX("$\\beta_3$")
    ),
    minor_breaks = NULL
  )

ggsave(filename = "figures/additional-illustrations/proportional-slowing.pdf",
       device = "pdf",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)
