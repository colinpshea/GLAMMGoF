# tests/testthat/test-lognormal.R
#
# Tests for the lognormal (residual + random-effect) bias correction added to
# jensen_correct() and bias_precision(). The critical properties under test:
#
#   1. A lognormal model (log(y) ~ .) recovers the arithmetic mean only after
#      the residual variance is included in the correction (the Monte Carlo test).
#   2. The residual term is added ONLY when the response is natural-log
#      transformed - NOT for log-link models and NOT for gaussian(link = "log"),
#      where the residual sits outside exp() (the "trap" regression tests).
#   3. The new error paths fire (identity-link ambiguity, non-natural log base,
#      tmb + lognormal).
#
# Simulators are self-contained so the tests do not depend on the shipped
# countData / countModel_* objects.

# ---------------------------------------------------------------------------
# Simulation helpers
# ---------------------------------------------------------------------------

# Lognormal LMM:  log(Y) = b0 + bx*x + b_site + eps
#   b_site ~ N(0, sd_re^2), eps ~ N(0, sd_resid^2)
#   => E[Y] = exp(b0 + bx*x + (sd_re^2 + sd_resid^2) / 2)
sim_lognormal <- function(seed = 123, n_site = 50, n_per = 20,
                          b0 = 2, bx = 0.3, sd_re = 0.5, sd_resid = 0.4) {
  set.seed(seed)
  N      <- n_site * n_per
  site   <- factor(rep(seq_len(n_site), each = n_per))
  b_site <- rnorm(n_site, 0, sd_re)
  x      <- rnorm(N)
  eta    <- b0 + bx * x + b_site[as.integer(site)]
  logy   <- eta + rnorm(N, 0, sd_resid)
  data.frame(y = exp(logy), x = x, site = site)
}

# Poisson GLMM (log link, no residual-inside-exp term)
sim_count <- function(seed = 202, n_site = 50, n_per = 20,
                      b0 = 1, bx = 0.4, sd_re = 0.5) {
  set.seed(seed)
  N      <- n_site * n_per
  site   <- factor(rep(seq_len(n_site), each = n_per))
  b_site <- rnorm(n_site, 0, sd_re)
  x      <- rnorm(N)
  eta    <- b0 + bx * x + b_site[as.integer(site)]
  data.frame(y = rpois(N, exp(eta)), x = x, site = site)
}

# gaussian(link = "log"): mean is exp(eta) but residual is additive on the
# RESPONSE scale, so it must NOT enter the correction.
sim_loglink_gaussian <- function(seed = 303, n_site = 50, n_per = 20,
                                 b0 = 2, bx = 0.3, sd_re = 0.5, sd_resid = 0.5) {
  set.seed(seed)
  N      <- n_site * n_per
  site   <- factor(rep(seq_len(n_site), each = n_per))
  b_site <- rnorm(n_site, 0, sd_re)
  x      <- rnorm(N)
  mu     <- exp(b0 + bx * x + b_site[as.integer(site)])
  y      <- mu + rnorm(N, 0, sd_resid)
  data.frame(y = pmax(y, 1e-3), x = x, site = site)
}


# ---------------------------------------------------------------------------
# 1. Monte Carlo recovery: the whole point of the residual term
# ---------------------------------------------------------------------------

test_that("lognormal LMM correction recovers the arithmetic mean; uncorrected underestimates", {
  skip_if_not_installed("glmmTMB")

  dat <- sim_lognormal(seed = 1, n_site = 60, n_per = 40)  # 2400 obs, tight MC error
  m   <- glmmTMB::glmmTMB(log(y) ~ x + (1 | site), family = gaussian, data = dat)

  factor_hat  <- jensen_correct(m)                              # exp((s2_re + s2_resid)/2)
  eta_hat     <- predict(m, type = "link", re.form = ~0)        # marginal log-scale mean
  uncorrected <- exp(eta_hat)                                   # geometric mean (biased low)
  corrected   <- uncorrected * factor_hat                       # arithmetic mean

  # Corrected sample mean matches the empirical arithmetic mean of Y
  expect_equal(mean(corrected), mean(dat$y), tolerance = 0.08)

  # Uncorrected is systematically below the arithmetic mean (Jensen's inequality)
  expect_lt(mean(uncorrected), mean(dat$y))

  # ... and by enough that the correction is doing real work here (~15-20%)
  expect_gt(mean(dat$y) / mean(uncorrected) - 1, 0.10)
})


# ---------------------------------------------------------------------------
# 2. jensen_correct(): residual included for lognormal, excluded otherwise
# ---------------------------------------------------------------------------

test_that("jensen_correct includes the residual variance for a lognormal LMM", {
  skip_if_not_installed("glmmTMB")

  dat <- sim_lognormal(seed = 2)
  m   <- glmmTMB::glmmTMB(log(y) ~ x + (1 | site), family = gaussian, data = dat)

  s2_re    <- glmmTMB::VarCorr(m)$cond$site[1, 1]
  s2_resid <- sigma(m)^2
  expect_equal(jensen_correct(m), exp((s2_re + s2_resid) / 2), tolerance = 1e-6)

  # Must differ from the RE-only factor by exactly the residual contribution
  expect_false(isTRUE(all.equal(jensen_correct(m), exp(s2_re / 2))))
})

test_that("jensen_correct on a lognormal LM (no RE) uses the residual term only", {
  dat <- sim_lognormal(seed = 3)
  m   <- lm(log(y) ~ x, data = dat)
  expect_equal(jensen_correct(m), exp(sigma(m)^2 / 2), tolerance = 1e-6)
})

test_that("log-link GLMM correction excludes the residual (trap regression)", {
  skip_if_not_installed("glmmTMB")

  dat <- sim_count(seed = 4)
  m   <- glmmTMB::glmmTMB(y ~ x + (1 | site), family = poisson, data = dat)

  s2_re <- glmmTMB::VarCorr(m)$cond$site[1, 1]
  expect_equal(jensen_correct(m), exp(s2_re / 2), tolerance = 1e-6)  # RE variance only
})

test_that("gaussian(link='log') gets RE-only correction, no residual (trap regression)", {
  skip_if_not_installed("glmmTMB")

  dat <- sim_loglink_gaussian(seed = 5)
  m   <- glmmTMB::glmmTMB(y ~ x + (1 | site),
                          family = gaussian(link = "log"), data = dat)

  s2_re <- glmmTMB::VarCorr(m)$cond$site[1, 1]
  expect_equal(jensen_correct(m), exp(s2_re / 2), tolerance = 1e-6)

  # The residual exists (sigma > 0) but must be excluded: the corrected factor
  # must NOT equal the residual-inclusive one.
  expect_false(isTRUE(all.equal(jensen_correct(m),
                                exp((s2_re + sigma(m)^2) / 2))))
})

test_that("lognormal factor agrees across glmmTMB and lme4 backends", {
  skip_if_not_installed("glmmTMB")
  skip_if_not_installed("lme4")

  dat   <- sim_lognormal(seed = 6)
  m_tmb <- glmmTMB::glmmTMB(log(y) ~ x + (1 | site), family = gaussian, data = dat)
  m_mer <- lme4::lmer(log(y) ~ x + (1 | site), data = dat, REML = FALSE)
  expect_equal(jensen_correct(m_tmb), jensen_correct(m_mer), tolerance = 0.02)
})


# ---------------------------------------------------------------------------
# 3. Error / edge paths
# ---------------------------------------------------------------------------

test_that("auto errors on an identity-link untransformed response", {
  skip_if_not_installed("lme4")

  dat  <- sim_lognormal(seed = 7)
  m_lm <- lm(y ~ x, data = dat)
  expect_error(jensen_correct(m_lm), "identity link")

  m_id <- lme4::lmer(y ~ x + (1 | site), data = dat)   # Gaussian identity LMM
  expect_error(jensen_correct(m_id), "identity link")
})

test_that("type='lognormal' forces the residual term on a pre-logged column", {
  dat      <- sim_lognormal(seed = 8)
  dat$logy <- log(dat$y)
  m_col    <- lm(logy ~ x, data = dat)              # transform invisible to formula

  expect_error(jensen_correct(m_col))               # auto cannot tell -> ambiguous
  expect_equal(jensen_correct(m_col, type = "lognormal"),
               exp(sigma(m_col)^2 / 2), tolerance = 1e-6)
})

test_that("non-natural log bases are rejected", {
  skip_if_not_installed("lme4")

  dat <- sim_lognormal(seed = 9)
  m10 <- lme4::lmer(log10(y) ~ x + (1 | site), data = dat)
  expect_error(jensen_correct(m10), "natural-log")

  m1p <- lme4::lmer(log1p(y) ~ x + (1 | site), data = dat)
  expect_error(jensen_correct(m1p), "natural-log")
})

test_that("log-link GLM without random effects returns factor 1 with a warning", {
  dat <- sim_count(seed = 10)
  mp  <- glm(y ~ x, family = poisson, data = dat)
  expect_warning(f <- jensen_correct(mp), "correction factor is 1")
  expect_equal(f, 1)
})

test_that("applied correction is multiplicative on both scales", {
  skip_if_not_installed("glmmTMB")

  dat <- sim_lognormal(seed = 11)
  m   <- glmmTMB::glmmTMB(log(y) ~ x + (1 | site), family = gaussian, data = dat)
  f   <- jensen_correct(m)
  p   <- c(1, 2, 5)

  resp <- jensen_correct(m, predictions = p, scale = "response")
  expect_equal(resp$predictions_adjusted, p * f)

  link <- jensen_correct(m, predictions = log(p), scale = "link")
  expect_equal(link$predictions_adjusted, p * f)
})


# ---------------------------------------------------------------------------
# 4. bias_precision(): lognormal support end to end
# ---------------------------------------------------------------------------

test_that("bias_precision runs on a lognormal model and returns finite metrics", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  dat <- sim_lognormal(seed = 12, n_site = 40, n_per = 20)
  m   <- glmmTMB::glmmTMB(log(y) ~ x + (1 | site), family = gaussian, data = dat)

  res <- bias_precision(nReps = 20, testModel = m, testData = dat,
                        method = "holdout", DHARMaPlot = FALSE,
                        seed = 1, verbose = FALSE)
  expect_true(all(is.finite(res$bias_precision_summary$mn)))

  # Sanity: RMedAE should be a sensible percentage, not exploded by a scale error
  rmedae_out <- res$bias_precision_summary$mn[
    res$bias_precision_summary$Metric == "RMedAE" &
    res$bias_precision_summary$Group  == "Out-of-sample performance"]
  expect_true(rmedae_out > 0 && rmedae_out < 500)
})

test_that("bias_adjust='manual' moves RBIAS toward zero for a lognormal model", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  dat <- sim_lognormal(seed = 13, n_site = 40, n_per = 20)
  m   <- glmmTMB::glmmTMB(log(y) ~ x + (1 | site), family = gaussian, data = dat)

  out_rbias <- function(res) {
    s <- res$bias_precision_summary
    s$mn[s$Metric == "RBIAS" & s$Group == "Out-of-sample performance"]
  }

  none <- bias_precision(nReps = 25, testModel = m, testData = dat,
                         DHARMaPlot = FALSE, bias_adjust = "none",
                         seed = 1, verbose = FALSE)
  man  <- bias_precision(nReps = 25, testModel = m, testData = dat,
                         DHARMaPlot = FALSE, bias_adjust = "manual",
                         seed = 1, verbose = FALSE)

  expect_lt(out_rbias(none), 0)                        # uncorrected underestimates
  expect_lt(abs(out_rbias(man)), abs(out_rbias(none))) # correction pulls toward 0
})

test_that("bias_adjust='tmb' errors for lognormal models", {
  skip_if_not_installed("glmmTMB")

  dat <- sim_lognormal(seed = 14, n_site = 30, n_per = 15)
  m   <- glmmTMB::glmmTMB(log(y) ~ x + (1 | site), family = gaussian, data = dat)
  expect_error(
    bias_precision(nReps = 5, testModel = m, testData = dat,
                   DHARMaPlot = FALSE, bias_adjust = "tmb", verbose = FALSE),
    "lognormal")
})

test_that("bias_precision rejects non-natural log responses", {
  skip_if_not_installed("glmmTMB")

  dat <- sim_lognormal(seed = 15, n_site = 30, n_per = 15)
  m   <- glmmTMB::glmmTMB(log10(y) ~ x + (1 | site), family = gaussian, data = dat)
  expect_error(
    bias_precision(nReps = 5, testModel = m, testData = dat,
                   DHARMaPlot = FALSE, verbose = FALSE),
    "natural-log")
})

test_that("Poisson GLMM metrics are unaffected by the lognormal changes (regression)", {
  skip_on_cran()
  skip_if_not_installed("glmmTMB")

  dat <- sim_count(seed = 16, n_site = 40, n_per = 20)
  m   <- glmmTMB::glmmTMB(y ~ x + (1 | site), family = poisson, data = dat)

  # Runs cleanly and is NOT treated as lognormal (no exp()-of-response rescaling)
  res <- bias_precision(nReps = 20, testModel = m, testData = dat,
                        DHARMaPlot = FALSE, bias_adjust = "none",
                        seed = 1, verbose = FALSE)
  expect_true(all(is.finite(res$bias_precision_summary$mn)))
})
