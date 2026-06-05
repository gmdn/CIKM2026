# CIKM2026

This repository contains the source code used to reproduce the experiments reported in the paper:

> **Adaptive Bayesian Prevalence Estimation Under Limited Assessment Budgets**

The study investigates Bayesian sequential prevalence estimation using a Beta--Bernoulli model and Highest Posterior Density (HPD) intervals. The objective is to estimate the prevalence of a binary property while minimizing assessment effort through uncertainty-aware stopping rules.

## Overview

Estimating the prevalence of a target class is a common problem in knowledge graph quality assessment, information retrieval, systematic reviews, legal discovery, and other large-scale assessment tasks where labels are expensive to obtain.

This repository implements a Bayesian sequential framework in which:

* observations are modeled as Bernoulli trials;
* prevalence is assigned a Beta prior;
* posterior uncertainty is quantified through HPD intervals;
* assessment stops automatically when posterior uncertainty falls below a predefined threshold.

The experiments evaluate how prevalence affects:

1. Assessment effort.
2. Posterior coverage.
3. Estimation accuracy.
4. Number of positive observations collected before stopping.

## Repository Contents

```text
.
├── simulation.R      # Main simulation script
├── helper.R          # Supporting functions
├── figures/          # Generated figures
├── results/          # Simulation outputs
└── README.md
```

## Requirements

The code was developed in R (version 4.3+ recommended).

Required packages include:

```r
dplyr
purrr
tibble
tidyr
ggplot2
```

Install them with:

```r
install.packages(
  c(
    "dplyr",
    "purrr",
    "tibble",
    "tidyr",
    "ggplot2"
  )
)
```

If an `renv.lock` file is provided, the recommended installation procedure is:

```r
renv::restore()
```

## Reproducing the Experiments

Run:

```r
source("simulation.R")
```

The script executes the complete simulation study.

The experiments evaluate:

```r
p_true_values = c(
  0.5,
  0.1,
  0.01,
  0.001,
  0.0001
)

hpd_width_thresholds = c(
  0.1,
  0.01,
  0.001,
  0.0001
)
```

using:

```r
n_reps = 1000
max_samples = 10000
```

for each parameter combination.

## Output

The simulation produces summary statistics including:

* Mean number of assessments required to stop.
* Coverage of 95% HPD intervals.
* Mean HPD width.
* Absolute estimation error.
* Relative estimation error.
* Mean number of positive observations.

The paper reports results using:

* Assessment effort.
* Coverage.
* Relative error.
* Positive observations at stopping.

## Main Findings

The experiments show that:

* Assessment effort is strongly affected by prevalence.
* Highly imbalanced prevalence levels often require fewer assessments.
* Coverage remains close to the nominal 95% level.
* Relative error increases substantially for extremely rare events.
* Stopping often occurs after very few positive observations in rare-event scenarios.

## Reproducibility

All simulations are initialized with a fixed random seed:

```r
seed = 123
```

allowing exact replication of the reported results.

## Citation

If you use this code, please cite:

Anonymized 

## License

This repository is released for research and educational purposes.

Please refer to the LICENSE file for details.
