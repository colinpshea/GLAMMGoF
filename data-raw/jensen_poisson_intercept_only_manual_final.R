# ══════════════════════════════════════════════════════════════════════════════
# Jensen sigma-sweep simulation for the GLAMMGoF vignette figure
#
# This script:
#   1. Runs the sigma-sweep simulation across three scenarios
#   2. Saves the summarised results to inst/extdata/sim_jensen_sigma_sweep.rds
#
# Runtime: approximately X minutes at 100 datasets x 100 CV reps.
# The .rds file is loaded and plotted by the vignette; this script does not
# need to run at package install or CRAN check time.
#
# Run from the package root directory.
# ══════════════════════════════════════════════════════════════════════════════

# Install 'GLAMMGoF' in R:
#install.packages('GLAMMGoF', repos = c('https://colinpshea.r-universe.dev', 'https://cloud.r-project.org'))

library(glmmTMB)
library(GLAMMGoF)
library(dplyr)
library(ggplot2)

# ══════════════════════════════════════════════════════════════════════════════
# Jensen's inequality bias: intercept-only Poisson GLMM simulation
# Purpose: isolate contribution of fixed effect estimation to the gap between
# theoretical and observed RBIAS by removing all covariate complexity.
# With y ~ 1 + (1|site), only beta0 and sigma^2 are estimated.
#
# Outer loop over n_datasets independent datasets per sigma value eliminates
# single-dataset noise that caused wobbly lines at high N.
# ══════════════════════════════════════════════════════════════════════════════

set.seed(42)

# ── Design ────────────────────────────────────────────────────────────────────
sigma_grid  <- seq(0.10, 1.00, length.out = 12)
n_datasets  <- 100    # 100 outer loop: independent datasets per sigma value
n_reps_per  <- 100    # 100 inner loop: CV replicates per dataset
intercept   <- 1.5
run_manual <- TRUE
run_tmb <- FALSE
n_datasets_tmb <- 25 # 25
n_reps_tmb <- 25 # 25

sample_sizes <- list(
 "N = 20000" = list(n_sites = 50, n_years = 20, n_obs = 20)
)

scenarios <- c("glm_no_re", "glmm_one_re", "glmm_two_re")

safe_bind <- function(x) bind_rows(Filter(Negate(is.null), x))

# ── Data simulator (intercept only) ──────────────────────────────────────────
simulate_dat <- function(n_sites, n_years, n_obs, re_sd, scenario) {

  n_total <- n_sites * n_years * n_obs
  site_id <- rep(rep(seq_len(n_sites), each = n_years), each = n_obs)
  year_id <- rep(rep(seq_len(n_years), times = n_sites), each = n_obs)

  beta0 <- switch(scenario,
    glm_no_re    = intercept,
    glmm_one_re  = intercept - re_sd^2 / 2,
    glmm_two_re  = intercept - re_sd^2
  )

  # Intercept only -- no covariates
  log_mu <- rep(beta0, n_total)

  if (scenario == "glmm_one_re") {
    re_site <- rnorm(n_sites, 0, re_sd)
    log_mu  <- log_mu + re_site[site_id]
  } else if (scenario == "glmm_two_re") {
    re_site <- rnorm(n_sites, 0, re_sd)
    re_year <- rnorm(n_years, 0, re_sd)
    log_mu  <- log_mu + re_site[site_id] + re_year[year_id]
  }

  y <- rpois(n_total, lambda = exp(log_mu))

  data.frame(y = y,
             site = factor(site_id),
             year = factor(year_id))
}

# ── Model fitter (intercept only) ────────────────────────────────────────────
fit_model <- function(dat, scenario) {
  tryCatch({
    if (scenario == "glm_no_re") {
      glmmTMB(y ~ 1, data = dat, family = poisson)
    } else if (scenario == "glmm_one_re") {
      glmmTMB(y ~ 1 + (1 | site), data = dat, family = poisson)
    } else {
      glmmTMB(y ~ 1 + (1 | site) + (1 | year), data = dat, family = poisson)
    }
  }, error = function(e) NULL)
}

# ── Sweep function ────────────────────────────────────────────────────────────
run_sweep <- function(scenario, bias_adj, ss_name, ss_params) {

  # Use fewer datasets and reps for tmb to keep runtime manageable
  n_ds <- if (bias_adj == "tmb") n_datasets_tmb else n_datasets
  n_rep <- if (bias_adj == "tmb") n_reps_tmb else n_reps_per

  safe_bind(lapply(sigma_grid, function(re_sd) {

    dataset_results <- safe_bind(lapply(seq_len(n_ds), function(d) {

      dat <- simulate_dat(ss_params$n_sites, ss_params$n_years,
                          ss_params$n_obs, re_sd, scenario)
      m <- fit_model(dat, scenario)
      if (is.null(m)) return(NULL)

      vc <- VarCorr(m)$cond
      re_sd_est <- if (length(vc) > 0) sqrt(sum(sapply(vc, function(v) v[1,1]))) else 0

      cv_out <- tryCatch(
        bias_precision(nReps = n_rep,
                       testModel = m,
                       testData = dat,
                       method = "holdout",
                       propTrain = 0.8,
                       DHARMaPlot = FALSE,
                       bias_adjust = bias_adj, # ← use the argument!
                       verbose = FALSE),
        error = function(e) NULL
      )
      if (is.null(cv_out)) return(NULL)

      rbias_out <- cv_out$bias_precision_results |>
        filter(Metric == "RBIAS", Group == "Out-of-sample performance") |>
        pull(value)

      data.frame(
        dataset_id = d,
        re_sd_est = round(re_sd_est, 3),
        rbias_mean = mean(rbias_out, na.rm = TRUE)
      )
    }))

    if (is.null(dataset_results) || nrow(dataset_results) == 0) return(NULL)

    data.frame(
      sigma_true = re_sd,
      re_sd_est = mean(dataset_results$re_sd_est, na.rm = TRUE),
      rbias_mean = mean(dataset_results$rbias_mean, na.rm = TRUE),
      rbias_sd = sd(dataset_results$rbias_mean, na.rm = TRUE),
      scenario = scenario,
      bias_adjust = bias_adj, # ← store which method was used
      sample_size = ss_name,
      n_datasets = nrow(dataset_results)
    )
  }))
}

# ── Run ───────────────────────────────────────────────────────────────────────
all_results <- list()

for (ss_name in names(sample_sizes)) {
  ss <- sample_sizes[[ss_name]]
  message("\n── Sample size: ", ss_name, " ──")

  for (sc in scenarios) {
    message("  Scenario: ", sc)

    message("    bias_adjust = none")
    all_results[[paste(ss_name, sc, "none", sep = "_")]] <-
      run_sweep(sc, "none", ss_name, ss)

    if (run_manual && sc != "glm_no_re") {
      message("    bias_adjust = manual")
      all_results[[paste(ss_name, sc, "manual", sep = "_")]] <-
        run_sweep(sc, "manual", ss_name, ss)
    }

    if (run_tmb && sc != "glm_no_re") {
      message("    bias_adjust = tmb (slow...)")
      all_results[[paste(ss_name, sc, "tmb", sep = "_")]] <-
        run_sweep(sc, "tmb", ss_name, ss)
    }
  }
}

results <- safe_bind(all_results) |>
  mutate(
    lo = rbias_mean - 1.96 * rbias_sd / sqrt(n_datasets),
    hi = rbias_mean + 1.96 * rbias_sd / sqrt(n_datasets),
    scenario = factor(scenario,
      levels = c("glm_no_re", "glmm_one_re", "glmm_two_re"),
      labels = c("GLM (no RE)", "GLMM: site RE", "GLMM: site + year RE")),
    sample_size = factor(sample_size, levels = "N = 20000")
  )

# Dynamic factor levels and palette based on toggles
active_methods <- c("none",
                    if (run_manual) "manual",
                    if (run_tmb)    "tmb")
active_labels  <- c('No correction  ("none")',
                    if (run_manual) 'Analytical  ("manual")',
                    if (run_tmb)    'TMB  ("tmb")')

results <- results |>
  mutate(bias_adjust = factor(bias_adjust,
                              levels = active_methods,
                              labels = active_labels))

pal_full <- c(
  'No correction  ("none")'   = "#D85A30",
  'Analytical  ("manual")'    = "#185FA5",
  'TMB  ("tmb")'              = "#3DAA6A"
)
pal <- pal_full[active_labels]

# ── Summary table ─────────────────────────────────────────────────────────────
results |>
  select(sample_size, scenario, sigma_true, re_sd_est, rbias_mean) |>
  as_tibble() |>
  print(digits = 3, n = Inf)

# ── Figure ────────────────────────────────────────────────────────────────────
ggplot(results,
       aes(x = sigma_true, y = rbias_mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10,
              fill = "#D85A30", colour = NA) +
  geom_line(colour = "#D85A30", linewidth = 0.8) +
  geom_point(colour = "#D85A30", size = 1.8) +
  scale_x_continuous(
    name   = expression(paste("True RE ", sigma)),
    breaks = round(sigma_grid, 2)
  ) +
  labs(
    y        = "Mean out-of-sample RBIAS  [(predicted \u2212 observed) / observed]",
    title    = "Jensen's inequality bias: intercept-only Poisson GLMM",
    subtitle = paste0("Poisson  |  intercept only (y ~ 1)  |  ",
                      n_datasets, " datasets \u00d7 ", n_reps_per,
                      " CV reps  | N = 20,000")
  ) +
  facet_grid(sample_size ~ scenario) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(size = 9)
  )


# ── Figure ────────────────────────────────────────────────────────────────────
method_label <- paste(
  paste0(n_datasets, " datasets \u00d7 ", n_reps_per, " CV reps"), " | ",
  if (run_manual) '"manual" uses known true \u03c3' else NULL,
  if (run_tmb)    "tmb included" else NULL)

finplot <- ggplot(results,
       aes(x = sigma_true, y = rbias_mean,
           colour = bias_adjust, fill = bias_adjust)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = pal) +
  scale_fill_manual(values = pal) +
  scale_x_continuous(
    name   = expression(paste("True RE ", sigma)),
    breaks = round(sigma_grid, 2)
  ) +
  scale_y_continuous(limits = c(-65, 10), breaks = seq(-65, 10, 5), expand = expansion(add = c(0, 0))) +
  labs(
    y        = "Mean out-of-sample RBIAS  [(predicted \u2212 observed) / observed]",
    title    = "Jensen's inequality bias: intercept-only Poisson GLMM",
    subtitle = paste0("Poisson  |  intercept only (y ~ 1)  |  ",
                      n_datasets, " datasets \u00d7 ", n_reps_per,
                      " CV reps  | N = 20,000")
  ) +
  facet_grid(sample_size ~ scenario) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(colour = "grey40", size = 9),
    strip.background = element_rect(fill = "grey95"),
    strip.text       = element_text(size = 9)
  )

ggsave("jp_int_manual_072926.png", dpi = 600, height = 8, width = 12)

save.image("jensen_poisson_int_only_manual_final_072926.RData")

# ══════════════════════════════════════════════════════════════════════════════
# Prepare Jensen sigma-sweep simulation results for the vignette
# Source: jensen_poisson_int_only_manual_final_072926.RData
# ══════════════════════════════════════════════════════════════════════════════

# results is already summarised (one row per sigma x scenario x bias_adjust)
# so we just need to strip it down to what the figure needs and store metadata

sim_jensen_sigma_sweep <- list(
  results_df = results |>
    dplyr::select(sigma_true, scenario, bias_adjust,
                  rbias_mean, lo, hi, n_datasets, re_sd_est),
  metadata = list(
    design = "Intercept-only Poisson GLMM",
    n_datasets = 100,
    n_reps_per = 100,
    sigma_grid = seq(0.10, 1.00, length.out = 12),
    scenarios = c("GLM (no RE)", "GLMM: site RE", "GLMM: site + year RE"),
    n_total = 20000,
    seed = 42,
    session_info = sessionInfo()
  )
)

## save to inst/extdata
saveRDS(sim_jensen_sigma_sweep,
        "C:/Users/Colin.Shea/OneDrive - Florida Fish and Wildlife Conservation/R Packages/GLAMMGoF/inst/extdata/sim_jensen_sigma_sweep.rds", compress = "xz")
