library(dplyr)
library(purrr)
library(tibble)

#source("helper_new_stopping_rule.R")
source("helper.R")

run_multi_p_study <- function(
    p_true_values = c(0.6, 0.9),
    hpd_width_thresholds = c(0.1, 0.01, 0.001),
    n_reps = 10,
    alpha_prior = 1,
    beta_prior = 1,
    cred_mass = 0.95,
    max_samples = 1000,
    seed = 123) {
  
  set.seed(seed)
  
  scenario_grid <- tidyr::expand_grid(
    p_true = p_true_values,
    hpd_width_threshold = hpd_width_thresholds
  )
  
  results <- purrr::pmap_dfr(
    scenario_grid,
    \(p_true, hpd_width_threshold) {
      run_simulation_study(
        n_reps = n_reps,
        p_true = p_true,
        alpha_prior = alpha_prior,
        beta_prior = beta_prior,
        hpd_width_threshold = hpd_width_threshold,
        cred_mass = cred_mass,
        max_samples = max_samples,
        seed = NULL
      )$results |>
        dplyr::mutate(
          p_true = p_true,
          hpd_width_threshold = hpd_width_threshold
        )
    }
  )
  
  summary <- results |>
    dplyr::group_by(p_true, hpd_width_threshold) |>
    dplyr::summarise(
      mean_sample_size = mean(n_samples),
      median_sample_size = median(n_samples),
      coverage_rate = mean(covered),
      mean_hpd_width = mean(hpd_width),
      mean_absolute_error = mean(abs(posterior_mean - p_true)),
      median_absolute_error = median(abs(posterior_mean - p_true)),
      mean_relative_error = mean(relative_error),
      mean_successes = mean(successes),
      stop_rate = mean(stopped_by_rule),
      .groups = "drop"
    )
  
  list(
    results = results,
    summary = summary
  )
}


multi_p_out <- run_multi_p_study(
  p_true_values = c(0.5, 0.1, 0.01, 0.001, 0.0001),
  hpd_width_thresholds = c(0.1, 0.01, 0.001, 0.0001),
  n_reps = 1000,
  max_samples = 10000,
  seed = 123
)

#multi_p_out$summary


plot_assessment <- multi_p_out$summary |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = p_true,
      y = mean_sample_size,
      color = factor(hpd_width_threshold),
      group = hpd_width_threshold
    )
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    title = "Samples needed to stop across prevalence levels",
    x = "True prevalence",
    y = "Mean samples to stop",
    color = "HPD width threshold"
  ) +
  ggplot2::theme_minimal(base_size = 13)
  
ggplot2::ggsave(
  filename = "figures/assessments.png",
  plot = plot_assessment,
  width = 6,
  height = 4,
  dpi = 300
)

plot_coverage <- multi_p_out$summary |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = p_true,
      y = coverage_rate,
      color = factor(hpd_width_threshold),
      group = hpd_width_threshold
    )
  ) +
  ggplot2::geom_hline(
    yintercept = 0.95,
    linetype = "dashed",
    color = "gray40"
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_log10() +
  ggplot2::coord_cartesian(ylim = c(0.7, 1)) +
  ggplot2::labs(
    title = "Coverage across prevalence levels",
    x = "True prevalence",
    y = "Coverage rate",
    color = "HPD width threshold"
  ) +
  ggplot2::theme_minimal(base_size = 13)


ggplot2::ggsave(
  filename = "figures/coverage.png",
  plot = plot_coverage,
  width = 6,
  height = 4,
  dpi = 300
)


plot_absolute_error <- multi_p_out$summary |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = p_true,
      y = mean_absolute_error,
      color = factor(hpd_width_threshold),
      group = hpd_width_threshold
    )
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Estimation error at stopping",
    x = "True prevalence",
    y = "Mean absolute error",
    color = "HPD threshold"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  filename = "figures/absolute_error.png",
  plot = plot_absolute_error,
  width = 6,
  height = 4,
  dpi = 300
)


plot_relative_error <- multi_p_out$summary |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = p_true,
      y = mean_relative_error,
      color = factor(hpd_width_threshold),
      group = hpd_width_threshold
    )
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Estimation error at stopping",
    x = "True prevalence",
    y = "Mean relative error",
    color = "HPD threshold"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  filename = "figures/relative_error.png",
  plot = plot_relative_error,
  width = 6,
  height = 4,
  dpi = 300
)


plot_successes <- multi_p_out$summary |>
  ggplot2::ggplot(
    ggplot2::aes(
      x = p_true,
      y = mean_successes,
      color = factor(hpd_width_threshold),
      group = hpd_width_threshold
    )
  ) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    title = "Estimation successes at stopping",
    x = "True prevalence",
    y = "Mean successes",
    color = "HPD threshold"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  filename = "figures/successes.png",
  plot = plot_successes,
  width = 6,
  height = 4,
  dpi = 300
)


readr::write_csv(
  multi_p_out$summary,
  "results/summary.csv"
)


save.image(file='session.RData')
