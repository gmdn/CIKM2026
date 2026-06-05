run_single_simulation <- function(
    p_true,
    alpha_prior = 1,
    beta_prior = 1,
    hpd_width_threshold = 0.1,
    cred_mass = 0.95,
    max_samples = 1000) {
  successes <- 0L
  failures <- 0L
  n_samples <- 0L
  
  repeat {
    y_new <- rbinom(n = 1, size = 1, prob = p_true)
    
    n_samples <- n_samples + 1L
    
    if (y_new == 1) {
      successes <- successes + 1L
    } else {
      failures <- failures + 1L
    }
    
    alpha_post <- alpha_prior + successes
    beta_post <- beta_prior + failures
    
    hpd <- compute_hpd_beta(
      shape1 = alpha_post,
      shape2 = beta_post,
      cred_mass = cred_mass
    )
    
    if (hpd$width <= hpd_width_threshold || n_samples >= max_samples) {
      break
    }
  }
  
  tibble(
    p_true = p_true,
    alpha_prior = alpha_prior,
    beta_prior = beta_prior,
    n_samples = n_samples,
    successes = successes,
    failures = failures,
    alpha_post = alpha_post,
    beta_post = beta_post,
    posterior_mean = alpha_post / (alpha_post + beta_post),
    relative_error = abs(posterior_mean - p_true) / p_true,
    hpd_lower = hpd$lower,
    hpd_upper = hpd$upper,
    hpd_width = hpd$width,
    covered = p_true >= hpd$lower && p_true <= hpd$upper,
    stopped_by_rule = hpd$width <= hpd_width_threshold
  )
}

run_simulation_study <- function(
    n_reps,
    p_true,
    alpha_prior = 1,
    beta_prior = 1,
    hpd_width_threshold = 0.1,
    cred_mass = 0.95,
    max_samples = 1000,
    seed = 123) {
  set.seed(seed)
  
  results <- map_dfr(
    seq_len(n_reps),
    \(i) {
      run_single_simulation(
        p_true = p_true,
        alpha_prior = alpha_prior,
        beta_prior = beta_prior,
        hpd_width_threshold = hpd_width_threshold,
        cred_mass = cred_mass,
        max_samples = max_samples
      ) |>
        mutate(replication = i)
    }
  )
  
  summary <- results |>
    summarise(
      n_reps = n(),
      p_true = first(p_true),
      alpha_prior = first(alpha_prior),
      beta_prior = first(beta_prior),
      mean_sample_size = mean(n_samples),
      median_sample_size = median(n_samples),
      min_sample_size = min(n_samples),
      max_sample_size = max(n_samples),
      sd_sample_size = sd(n_samples),
      coverage_rate = mean(covered),
      stop_rate = mean(stopped_by_rule),
      mean_hpd_width = mean(hpd_width)
    )
  
  list(
    results = results,
    summary = summary
  )
}

compute_hpd_beta <- function(shape1, shape2, cred_mass = 0.95) {
  interval_width <- function(prob_left) {
    qbeta(prob_left + cred_mass, shape1 = shape1, shape2 = shape2) -
      qbeta(prob_left, shape1 = shape1, shape2 = shape2)
  }
  
  opt <- optimize(
    f = interval_width,
    interval = c(0, 1 - cred_mass)
  )
  
  lower <- qbeta(
    opt$minimum,
    shape1 = shape1,
    shape2 = shape2
  )
  
  upper <- qbeta(
    opt$minimum + cred_mass,
    shape1 = shape1,
    shape2 = shape2
  )
  
  list(
    lower = lower,
    upper = upper,
    width = upper - lower
  )
}
