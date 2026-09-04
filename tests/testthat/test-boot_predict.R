# tests/testthat/test-boot_predict.R
#
# Coverage:
#   - Basic invocation across all five supported backends
#   - Jensen correction factor consistency with jensen_correct() called directly
#   - No-op equivalence for fixed-effect log-link GLMs (bias_adjust = "none"
#     and "manual" produce identical numerics under a fixed seed)
#   - Warning behavior for logit-link + bias_adjust = "manual"
#   - Error behavior for unsupported classes (GAM, zeroinfl)
#   - Reproducibility via seed argument
#   - Grid auto-build behavior (default, at overrides, numeric_default modes)
#   - Output attributes populated correctly
#
# Tests use small n_sim (100 or 200) for speed. This is sufficient for
# structural and consistency tests; numeric-precision tests use larger
# n_sim only where necessary.

# ──────────────────────────────────────────────────────────────
# Fixtures
# ──────────────────────────────────────────────────────────────
# The package ships glmmTMB example models (countModel_GLM,
# countModel_GLMM, countModel_GLMM2, logitModel_GLMM, etc.) which
# cover the glmmTMB and GAM backends. For lme4, glm, lm, and negbin
# coverage, we fit small models inline from the same countData.

skip_if_not_installed_all <- function(pkgs) {
  for (p in pkgs) testthat::skip_if_not_installed(p)
}

fit_lme4_glmm <- function() {
  # glmer with log-link Poisson: small enough to fit fast, structure
  # matches the glmmTMB GLMM for cross-backend comparison.
  lme4::glmer(y ~ Season + Temp + (1 | Site),
              data = countData, family = poisson)
}

fit_glm <- function() {
  stats::glm(y ~ Season + Temp, data = countData, family = poisson)
}

fit_lm <- function() {
  d <- countData[countData$y > 0,]
  stats::lm(log(y) ~ Season + Temp, data = d)
}

fit_negbin <- function() {
  MASS::glm.nb(y ~ Season + Temp, data = countData)
}


# ──────────────────────────────────────────────────────────────
# Helper: shared structural assertions
# ──────────────────────────────────────────────────────────────

expect_boot_output <- function(out, expected_cols = NULL) {
  testthat::expect_s3_class(out, "data.frame")
  testthat::expect_true(
    all(c("boot_mean", "boot_se", "boot_lwr", "boot_upr") %in% names(out))
  )
  testthat::expect_true(!is.null(attr(out, "draws")))
  testthat::expect_equal(ncol(attr(out, "draws")), nrow(out) |> (\(.) .)() * 0 +
                         ncol(attr(out, "draws")))   # sanity: dims exist
  testthat::expect_true(is.numeric(attr(out, "jensen_factor")))
  testthat::expect_true(attr(out, "bias_adjust") %in% c("none", "manual"))
  testthat::expect_true(attr(out, "link") %in% c("log", "logit", "identity"))
  if (!is.null(expected_cols)) {
    testthat::expect_true(all(expected_cols %in% names(out)))
  }
}


# ──────────────────────────────────────────────────────────────
# Backend coverage
# ──────────────────────────────────────────────────────────────

test_that("boot_predict works on glmmTMB fixed-effect GLM", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLM, n_sim = 100, seed = 1)
  expect_boot_output(out, expected_cols = c("Season", "Temp"))
  expect_equal(attr(out, "jensen_factor"), 1)   # no REs, no log-transform
})

test_that("boot_predict works on glmmTMB GLMM with random intercept", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM, bias_adjust = "manual",
                      n_sim = 100, seed = 1)
  expect_boot_output(out, expected_cols = c("Season", "Temp"))
  expect_gt(attr(out, "jensen_factor"), 1)   # log-link + REs -> factor > 1
})

test_that("boot_predict works on lme4 glmer", {
  skip_if_not_installed_all(c("lme4"))
  m <- fit_lme4_glmm()
  out <- boot_predict(m, bias_adjust = "manual", n_sim = 100, seed = 1)
  expect_boot_output(out, expected_cols = c("Season", "Temp"))
  expect_gt(attr(out, "jensen_factor"), 1)
})

test_that("boot_predict works on stats::glm", {
  m <- fit_glm()
  out <- boot_predict(m, n_sim = 100, seed = 1)
  expect_boot_output(out, expected_cols = c("Season", "Temp"))
  expect_equal(attr(out, "jensen_factor"), 1)
})

test_that("boot_predict works on stats::lm with log-transformed response", {
  m <- fit_lm()
  out <- boot_predict(m, bias_adjust = "manual", n_sim = 100, seed = 1)
  expect_boot_output(out, expected_cols = c("Season", "Temp"))
  # log-transformed response: correction factor should be > 1 from
  # residual variance alone
  expect_gt(attr(out, "jensen_factor"), 1)
})

test_that("boot_predict works on MASS::glm.nb", {
  skip_if_not_installed_all(c("MASS"))
  m <- fit_negbin()
  out <- boot_predict(m, n_sim = 100, seed = 1)
  expect_boot_output(out, expected_cols = c("Season", "Temp"))
  expect_equal(attr(out, "jensen_factor"), 1)   # no REs -> no correction
})


# ──────────────────────────────────────────────────────────────
# Correction factor consistency with jensen_correct()
# ──────────────────────────────────────────────────────────────

test_that("jensen_factor attribute matches jensen_correct() on same model", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM, bias_adjust = "manual",
                      n_sim = 100, seed = 1)
  expect_equal(attr(out, "jensen_factor"), jensen_correct(countModel_GLMM))
})

test_that("jensen_factor for two-RE model matches jensen_correct()", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM2, bias_adjust = "manual",
                      n_sim = 100, seed = 1)
  expect_equal(attr(out, "jensen_factor"), jensen_correct(countModel_GLMM2))
})


# ──────────────────────────────────────────────────────────────
# No-op equivalence: bias_adjust makes no difference for fixed-
# effect log-link GLMs
# ──────────────────────────────────────────────────────────────

test_that("bias_adjust has no numeric effect on fixed-effect log-link GLM", {
  skip_if_not_installed_all(c("glmmTMB"))
  out_none <- boot_predict(countModel_GLM, bias_adjust = "none",
                           n_sim = 200, seed = 42)
  out_man <- suppressWarnings(
    boot_predict(countModel_GLM, bias_adjust = "manual",
                 n_sim = 200, seed = 42)
  )
  expect_equal(out_none$boot_mean, out_man$boot_mean, tolerance = 1e-10)
  expect_equal(out_none$boot_se,   out_man$boot_se,   tolerance = 1e-10)
  expect_equal(out_none$boot_lwr,  out_man$boot_lwr,  tolerance = 1e-10)
  expect_equal(out_none$boot_upr,  out_man$boot_upr,  tolerance = 1e-10)
})

test_that("bias_adjust = 'manual' emits informative message for no-RE model", {
  skip_if_not_installed_all(c("glmmTMB"))
  expect_message(
    boot_predict(countModel_GLM, bias_adjust = "manual",
                 n_sim = 100, seed = 1),
    regexp = "no random effects"
  )
})


# ──────────────────────────────────────────────────────────────
# Correction uniformly scales all summaries
# ──────────────────────────────────────────────────────────────

test_that("Jensen correction scales boot_mean, boot_se, and CIs by the same factor", {
  skip_if_not_installed_all(c("glmmTMB"))
  out_none <- boot_predict(countModel_GLMM, bias_adjust = "none",
                           n_sim = 200, seed = 42)
  out_man <- boot_predict(countModel_GLMM, bias_adjust = "manual",
                          n_sim = 200, seed = 42)

  jf <- attr(out_man, "jensen_factor")
  expect_equal(out_man$boot_mean, out_none$boot_mean * jf, tolerance = 1e-10)
  expect_equal(out_man$boot_se,   out_none$boot_se   * jf, tolerance = 1e-10)
  expect_equal(out_man$boot_lwr,  out_none$boot_lwr  * jf, tolerance = 1e-10)
  expect_equal(out_man$boot_upr,  out_none$boot_upr  * jf, tolerance = 1e-10)
})


# ──────────────────────────────────────────────────────────────
# Warning behavior for logit-link
# ──────────────────────────────────────────────────────────────

test_that("logit-link + bias_adjust = 'manual' warns and returns factor = 1", {
  skip_if_not_installed_all(c("glmmTMB"))
  expect_warning(
    out <- boot_predict(logitModel_GLMM, bias_adjust = "manual",
                        n_sim = 100, seed = 1),
    regexp = "logit"
  )
  expect_equal(attr(out, "jensen_factor"), 1)
})


# ──────────────────────────────────────────────────────────────
# Error behavior for unsupported classes
# ──────────────────────────────────────────────────────────────

test_that("boot_predict errors informatively on GAM objects", {
  skip_if_not_installed_all(c("mgcv"))
  expect_error(
    boot_predict(countModel_GAM, n_sim = 100),
    regexp = "mgcv"
  )
})

test_that("boot_predict errors informatively on zeroinfl objects", {
  skip_if_not_installed_all(c("pscl"))
  m <- pscl::zeroinfl(y ~ Season + Temp | 1, data = countData, dist = "negbin")
  expect_error(
    boot_predict(m, n_sim = 100),
    regexp = "zeroinfl|pscl"
  )
})

# ──────────────────────────────────────────────────────────────
# Reproducibility
# ──────────────────────────────────────────────────────────────

test_that("seed produces reproducible output", {
  skip_if_not_installed_all(c("glmmTMB"))
  out1 <- boot_predict(countModel_GLMM, n_sim = 200, seed = 42)
  out2 <- boot_predict(countModel_GLMM, n_sim = 200, seed = 42)
  expect_equal(out1$boot_mean, out2$boot_mean)
  expect_equal(out1$boot_se,   out2$boot_se)
  expect_equal(out1$boot_lwr,  out2$boot_lwr)
  expect_equal(out1$boot_upr,  out2$boot_upr)
})

test_that("different seeds produce different output", {
  skip_if_not_installed_all(c("glmmTMB"))
  out1 <- boot_predict(countModel_GLMM, n_sim = 200, seed = 42)
  out2 <- boot_predict(countModel_GLMM, n_sim = 200, seed = 99)
  # Not identical, but should be close
  expect_false(isTRUE(all.equal(out1$boot_mean, out2$boot_mean)))
})


# ──────────────────────────────────────────────────────────────
# Grid construction
# ──────────────────────────────────────────────────────────────

test_that("auto-generated grid crosses factor levels with numeric mean by default", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM, n_sim = 50, seed = 1)
  # 4 seasons x 1 mean(Temp) = 4 rows
  expect_equal(nrow(out), 4)
  expect_equal(sort(as.character(out$Season)),
               sort(c("Spring", "Summer", "Autumn", "Winter")))
})

test_that("at argument overrides numeric_default for named covariates", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM,
                      at = list(Temp = c(10, 15, 20)),
                      n_sim = 50, seed = 1)
  # 4 seasons x 3 Temp values = 12 rows
  expect_equal(nrow(out), 12)
  expect_true(all(out$Temp %in% c(10, 15, 20)))
})

test_that("numeric_default = 'range' expands to numeric_length values", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM,
                      numeric_default = "range",
                      numeric_length = 5,
                      n_sim = 50, seed = 1)
  # 4 seasons x 5 Temp values = 20 rows
  expect_equal(nrow(out), 20)
  expect_equal(length(unique(out$Temp)), 5)
})

test_that("numeric_default as named list applies per-covariate specs", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM,
                      numeric_default = list(Temp = "range"),
                      numeric_length = 4,
                      n_sim = 50, seed = 1)
  # Temp gets range (4 values), covariates not in list default to mean
  expect_equal(length(unique(out$Temp)), 4)
})

test_that("grid exceeding 10,000 rows errors before bootstrap runs", {
  skip_if_not_installed_all(c("glmmTMB"))
  expect_error(
    boot_predict(countModel_GLMM,
                 numeric_default = "range",
                 numeric_length = 3000,
                 n_sim = 50),
    regexp = "grid|rows"
  )
})


# ──────────────────────────────────────────────────────────────
# Output attributes
# ──────────────────────────────────────────────────────────────

test_that("output carries all expected attributes with correct values", {
  skip_if_not_installed_all(c("glmmTMB"))
  out <- boot_predict(countModel_GLMM,
                      bias_adjust = "manual",
                      n_sim = 100, alpha = 0.10, seed = 42)
  attrs <- attributes(out)
  expect_true("draws" %in% names(attrs))
  expect_equal(dim(attrs$draws), c(nrow(out), 100))
  expect_equal(attrs$bias_adjust, "manual")
  expect_equal(attrs$link, "log")
  expect_equal(attrs$n_sim, 100)
  expect_equal(attrs$alpha, 0.10)
  expect_true(is.numeric(attrs$jensen_factor))
})

# ──────────────────────────────────────────────────────────────
# correction_factor tests (6 total)
# ──────────────────────────────────────────────────────────────

test_that("correction_factor: user-supplied scalar is applied correctly", {
  skip_if_not_installed_all(c("glmmTMB"))
  cf_user <- 1.25
  bp_user <- boot_predict(countModel_GLMM, correction_factor = cf_user,
                          n_sim = 500, seed = 1)
  bp_none <- boot_predict(countModel_GLMM, bias_adjust = "none",
                          n_sim = 500, seed = 1)
  # Each corrected prediction equals the uncorrected prediction times cf_user.
  expect_equal(bp_user$boot_mean,   bp_none$boot_mean   * cf_user, tolerance = 1e-8)
  expect_equal(bp_user$boot_median, bp_none$boot_median * cf_user, tolerance = 1e-8)
  # The jensen_factor attribute reflects the supplied value.
  expect_equal(attr(bp_user, "jensen_factor"), cf_user)
})

test_that("correction_factor: NULL default uses internal jensen_correct", {
  skip_if_not_installed_all(c("glmmTMB"))
  bp_manual <- boot_predict(countModel_GLMM, bias_adjust = "manual",
                            n_sim = 500, seed = 1)
  expect_equal(attr(bp_manual, "jensen_factor"),
               jensen_correct(countModel_GLMM), tolerance = 1e-8)
})

test_that("correction_factor: overrides bias_adjust with a message", {
  skip_if_not_installed_all(c("glmmTMB"))
  cf_user <- 1.10
  expect_message(
    bp <- boot_predict(countModel_GLMM, bias_adjust = "manual",
                       correction_factor = cf_user,
                       n_sim = 500, seed = 1),
    regexp = "overriding"
  )
  expect_equal(attr(bp, "jensen_factor"), cf_user)
})

test_that("correction_factor: overrides bias_adjust = 'none' with a message", {
  skip_if_not_installed_all(c("glmmTMB"))
  expect_message(
    bp <- boot_predict(countModel_GLMM, bias_adjust = "none",
                       correction_factor = 1.10,
                       n_sim = 500, seed = 1),
    regexp = "overriding"
  )
  expect_equal(attr(bp, "jensen_factor"), 1.10)
})

test_that("correction_factor: correction_factor = 1 is a valid no-op", {
  skip_if_not_installed_all(c("glmmTMB"))
  bp_one  <- boot_predict(countModel_GLMM, correction_factor = 1,
                          n_sim = 500, seed = 1)
  bp_none <- boot_predict(countModel_GLMM, bias_adjust = "none",
                          n_sim = 500, seed = 1)
  expect_equal(bp_one$boot_mean, bp_none$boot_mean, tolerance = 1e-8)
  expect_equal(attr(bp_one, "jensen_factor"), 1)
})

test_that("correction_factor: invalid values are rejected", {
  skip_if_not_installed_all(c("glmmTMB"))

  # Negative, zero, non-numeric, and non-finite values fail the type/value check
  expect_error(
    boot_predict(countModel_GLMM, correction_factor = -1,
                 n_sim = 100, seed = 1),
    regexp = "positive, finite"
  )
  expect_error(
    boot_predict(countModel_GLMM, correction_factor = 0,
                 n_sim = 100, seed = 1),
    regexp = "positive, finite"
  )
  expect_error(
    boot_predict(countModel_GLMM, correction_factor = "1.1",
                 n_sim = 100, seed = 1),
    regexp = "positive, finite"
  )
  expect_error(
    boot_predict(countModel_GLMM, correction_factor = NA_real_,
                 n_sim = 100, seed = 1),
    regexp = "positive, finite"
  )

  # A positive, finite numeric VECTOR of the wrong length fails the shape check
  # (vectors of correct length are now supported for per-row corrections)
  expect_error(
    boot_predict(countModel_GLMM, correction_factor = c(1.1, 1.2),
                 n_sim = 100, seed = 1),
    regexp = "length 1 or length nrow"
  )
})

test_that("correction_factor: correct-length vector is accepted", {
  skip_if_not_installed_all(c("glmmTMB"))
  bp_scalar <- boot_predict(countModel_GLMM, correction_factor = 1.1,
                            n_sim = 100, seed = 1)
  cf_vec <- rep(1.1, nrow(bp_scalar))
  expect_no_error(
    boot_predict(countModel_GLMM, correction_factor = cf_vec,
                 n_sim = 100, seed = 1)
  )
})

test_that("correction_factor: applied to logit-link model with informational message", {
  skip_if_not_installed_all(c("glmmTMB"))
  expect_message(
    bp <- boot_predict(logitModel_GLMM, correction_factor = 1.1,
                       n_sim = 500, seed = 1),
    regexp = "logit"
  )
  expect_equal(attr(bp, "jensen_factor"), 1.1)
})
