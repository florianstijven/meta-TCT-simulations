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
    color = "red"
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
      latex2exp::TeX("$g(t_1; \\gamma)$"),
      latex2exp::TeX("$g(t_2; \\gamma)$"),
      latex2exp::TeX("$g(t_3; \\gamma)$")
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

tibble(
  gamma = c(1.5, 0.75, 0.5),
  effect = c("acceleration", "slowing", "slowing")
) %>% 
  ggplot() +
  geom_abline(aes(slope = gamma, intercept = 0, color = as.factor(gamma))) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed") +
  scale_x_continuous(name = "Time in Active Group (Years)", breaks = 0:4, limits = c(0, 4)) +
  scale_y_continuous(name = latex2exp::TeX("$g(t; \\gamma) = \\gamma \\cdot t$"), breaks = -1:4, limits = c(-1, 4)) +
  scale_color_brewer(type = "qual", 
                     palette = 6,
                     name = NULL, 
                     labels = parse(text = expression(gamma == 1.5, gamma == 0.75, gamma == 0.5))) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 4),expand = FALSE) +
  ggtitle("Proportional Slowing")


ggsave(filename = "figures/additional-illustrations/proportional-slowing.pdf",
       device = "pdf",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)

## Proportional Slowing with an Offset ----------------------------------

tibble(
  x = c(0, 1, 0, 1, 0, 1),
  y = c(0, 1, 0, 0, 0, -1),
  xend = c(1, 4, 1, 4, 1, 4),
  yend = c(1, 3, 0, 3, -1, 4),
  effect = c("a", "a", "b", "b", "c", "c")
) %>%
  ggplot() +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, color = effect)) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed") +
  scale_x_continuous(name = "Time in Active Group (Years)", breaks = 0:4, limits = c(0, 4)) +
  scale_y_continuous(name = latex2exp::TeX("$g(t; \\gamma)$"), breaks = -1:4, limits = c(-1, 4)) +
  scale_color_brewer(type = "qual", 
                     palette = 6,
                     name = NULL) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 4), expand = FALSE) +
  ggtitle("Complex Time Mappings")


ggsave(filename = "figures/additional-illustrations/proportional-slowing-offset.pdf",
       device = "pdf",
       width = single_width,
       height = single_height,
       units = "cm",
       dpi = res)
