# ══════════════════════════════════════════════════════════════════════════════
# GLAMMGoF Supplement Simulation: bootstrap vs holdout under a 2x2x2 design
# Factors:
#   1. Model specification: nbinom2 (correct) vs Poisson (misspecified)
#   2. Observations per group: 3 (sparse) vs 25 (dense)
#   3. Resampling method: holdout vs bootstrap
# Response: nbinom2 with size = 1.5 in all cases
# Outputs: RBIAS under bias_adjust = "manual"; per-refit RE variance
# ══════════════════════════════════════════════════════════════════════════════

library(glmmTMB)
library(GLAMMGoF)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(42)

# ── Design ────────────────────────────────────────────────────────────────────
n_sites     <- 200
re_sd       <- 0.7            # variance = 0.49
size_nb     <- 1.5            # nbinom2 dispersion in generative process
intercept   <- 1.5
beta_x      <- 0.3

n_datasets  <- 30             # outer loop; increase for final run
n_reps_per  <- 50             # inner CV reps per dataset
n_reps_diag <- 50             # RE-variance diagnostic refits per dataset

design_grid <- expand.grid(
  fit_family = c("nbinom2", "poisson"),  # nbinom2 = correctly specified
  n_obs      = c(3, 25),                  # sparse vs dense
  stringsAsFactors = FALSE
)

# ── Data simulator (matches your existing one, response only) ─────────────────
simulate_dat <- function(n_sites, n_obs, re_sd, size_nb) {
  n_total <- n_sites * n_obs
  site_id <- rep(seq_len(n_sites), each = n_obs)
  x       <- rnorm(n_total)
  re_site <- rnorm(n_sites, 0, re_sd)
  beta0   <- intercept - re_sd^2 / 2   # centering so E[y] ~ exp(intercept)
  log_mu  <- beta0 + beta_x * x + re_site[site_id]
  y       <- rnbinom(n_total, mu = exp(log_mu), size = size_nb)
  data.frame(y = y, x = x, site = factor(site_id))
}

# ── Helper: fit and extract per-refit RE variance ─────────────────────────────
get_re_var <- function(fam, dat) {
  fit <- tryCatch(
    glmmTMB(y ~ x + (1 | site), family = fam, data = dat),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  VarCorr(fit)$cond[[1]][1, 1]
}

# ── Main simulation loop ──────────────────────────────────────────────────────
results_rbias  <- list()
results_revar  <- list()

t0 <- Sys.time()

for (g in seq_len(nrow(design_grid))) {
  
  fit_family <- design_grid$fit_family[g]
  n_obs      <- design_grid$n_obs[g]
  
  cat(sprintf("\n[cell %d/%d] family=%s, n_obs=%d\n",
              g, nrow(design_grid), fit_family, n_obs))
  
  for (d in seq_len(n_datasets)) {
    
    dat <- simulate_dat(n_sites, n_obs, re_sd, size_nb)
    
    # Fit full-data model in this cell's family
    fam_obj <- if (fit_family == "nbinom2") glmmTMB::nbinom2 else stats::poisson
    m_full  <- tryCatch(
      glmmTMB(y ~ x + (1 | site), family = fam_obj, data = dat),
      error = function(e) NULL
    )
    if (is.null(m_full)) next
    
    # 1. bias_precision under holdout + bootstrap, both with manual correction
    bp_h <- tryCatch(
      bias_precision(nReps = n_reps_per, testModel = m_full, testData = dat,
                     method = "holdout",   bias_adjust = "manual",
                     DHARMaPlot = FALSE, verbose = FALSE),
      error = function(e) NULL
    )
    bp_b <- tryCatch(
      bias_precision(nReps = n_reps_per, testModel = m_full, testData = dat,
                     method = "bootstrap", bias_adjust = "manual",
                     DHARMaPlot = FALSE, verbose = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(bp_h)) {
      results_rbias[[length(results_rbias) + 1]] <- bp_h$bias_precision_results |>
        filter(Metric == "RBIAS") |>
        mutate(fit_family = fit_family, n_obs = n_obs,
               method = "holdout", dataset = d)
    }
    if (!is.null(bp_b)) {
      results_rbias[[length(results_rbias) + 1]] <- bp_b$bias_precision_results |>
        filter(Metric == "RBIAS") |>
        mutate(fit_family = fit_family, n_obs = n_obs,
               method = "bootstrap", dataset = d)
    }
    
    # 2. Per-refit RE variance diagnostic
    n_row <- nrow(dat)
    revar_h <- revar_b <- numeric(n_reps_diag)
    for (i in seq_len(n_reps_diag)) {
      tr_h <- dat[sample.int(n_row, round(0.8 * n_row)), ]
      tr_b <- dat[sample.int(n_row, n_row, replace = TRUE), ]
      revar_h[i] <- get_re_var(fam_obj, tr_h)
      revar_b[i] <- get_re_var(fam_obj, tr_b)
    }
    results_revar[[length(results_revar) + 1]] <- data.frame(
      fit_family = fit_family, n_obs = n_obs, dataset = d,
      full_data_re_var = VarCorr(m_full)$cond[[1]][1, 1],
      method = rep(c("holdout", "bootstrap"), each = n_reps_diag),
      refit_re_var = c(revar_h, revar_b)
    )
  }
  cat(sprintf("  elapsed: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat(sprintf("\nTotal elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ── Assemble output tibbles ───────────────────────────────────────────────────
rbias_df <- bind_rows(results_rbias) |>
  mutate(
    fit_family = factor(fit_family, levels = c("nbinom2", "poisson"),
                        labels = c("Correctly specified (nbinom2)",
                                   "Misspecified (Poisson)")),
    n_obs = factor(n_obs, levels = c(3, 25),
                   labels = c("Sparse (3 obs/site)", "Dense (25 obs/site)"))
  )

revar_df <- bind_rows(results_revar) |>
  mutate(
    fit_family = factor(fit_family, levels = c("nbinom2", "poisson"),
                        labels = c("Correctly specified (nbinom2)",
                                   "Misspecified (Poisson)")),
    n_obs = factor(n_obs, levels = c(3, 25),
                   labels = c("Sparse (3 obs/site)", "Dense (25 obs/site)"))
  )

# Save outputs
saveRDS(rbias_df, "sim_rbias.rds")
saveRDS(revar_df, "sim_revar.rds")

# ── Panel A: RBIAS (manual) by method x family x n_obs ────────────────────────
panelA <- rbias_df |>
  filter(Group == "Out-of-sample performance") |>
  ggplot(aes(x = method, y = value, fill = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.4, alpha = 0.7) +
  facet_grid(n_obs ~ fit_family) +
  labs(x = NULL, y = "Out-of-sample RBIAS (%)\nbias_adjust = 'manual'",
       title = "A. RBIAS after Jensen correction, by method") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# ── Panel B: per-refit RE variance by method x family x n_obs ─────────────────
panelB <- revar_df |>
  ggplot(aes(x = method, y = refit_re_var, fill = method)) +
  geom_hline(aes(yintercept = full_data_re_var),
             linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.size = 0.4, alpha = 0.7) +
  facet_grid(n_obs ~ fit_family, scales = "free_y") +
  labs(x = NULL, y = "Per-refit RE variance",
       title = "B. Per-refit random-effect variance, by method") +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

# Combine panels (if you have patchwork installed)
 library(patchwork); panelA / panelB

fin <- panelA / panelB

ggsave("ggof_comp.png", dpi = 600, height = 7, width = 5)

# after running the simulation
sim_bootstrap_holdout <- list(
  rbias_df = rbias_df |> select(fit_family, n_obs, method, Group, Metric, value),
  revar_df = revar_df |> select(fit_family, n_obs, method, full_data_re_var, refit_re_var),
  metadata = list(
    n_datasets = 30, n_reps_per = 50, n_reps_diag = 50,
    re_sd = 0.7, size_nb = 1.5,
    seed = 42, session_info = sessionInfo()
  )
)
saveRDS(sim_bootstrap_holdout, "C:/Users/Colin.Shea/OneDrive - Florida Fish and Wildlife Conservation/R Packages/GLAMMGoF/inst/extdata/sim_bootstrap_holdout.rds")


