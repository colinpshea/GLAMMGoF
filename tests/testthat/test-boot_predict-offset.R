# Tests for boot_predict() offset handling and median-across-draws column.
# These tests cover the offset fix (formula-side, argument-side, user-supplied)
# and the boot_median companion to boot_mean.

# ---- Shared test data ------------------------------------------------------

make_offset_data <- function(n = 300, seed = 42) {
  set.seed(seed)
  effort <- runif(n, 0.5, 50)
  x      <- runif(n, -1, 1)
  mu     <- effort * exp(0.5 + 1.2 * x)
  y      <- rpois(n, mu)
  data.frame(y = y, x = x, effort = effort)
}

# ---- Group 1: Offset resolution paths --------------------------------------

test_that("boot_predict handles formula-side offset (glm)", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  gold <- stats::predict(mod, newdata = nd, type = "response")
  bp   <- boot_predict(mod, newdata = nd, n_sim = 3000, seed = 1)

  expect_equal(bp$boot_median, unname(gold), tolerance = 0.02)
  expect_equal(bp$boot_mean,   unname(gold), tolerance = 0.02)
})

test_that("boot_predict handles argument-side offset (glm)", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x, offset = log(effort), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  gold <- stats::predict(mod, newdata = nd, type = "response")
  bp   <- boot_predict(mod, newdata = nd, n_sim = 3000, seed = 1)

  expect_equal(bp$boot_median, unname(gold), tolerance = 0.02)
})

test_that("boot_predict handles no-offset model (backward compatibility)", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x, family = poisson, data = dat)  # no offset
  nd  <- data.frame(x = c(-0.5, 0, 0.5))

  gold <- stats::predict(mod, newdata = nd, type = "response")
  bp   <- boot_predict(mod, newdata = nd, n_sim = 3000, seed = 1)

  expect_equal(bp$boot_median, unname(gold), tolerance = 0.02)
  # Attribute records offset vector of zeros
  expect_equal(attr(bp, "offset"), rep(0, nrow(nd)))
})

# ---- Group 2: User-supplied offset argument --------------------------------

test_that("user-supplied `offset` argument overrides model offset", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  # Predict at unit exposure regardless of what newdata says
  bp   <- boot_predict(mod, newdata = nd, offset = log(1),
                       n_sim = 3000, seed = 1)
  gold <- exp(stats::coef(mod)[1] + stats::coef(mod)[2] * nd$x)

  expect_equal(bp$boot_median, unname(gold), tolerance = 0.02)
  expect_equal(attr(bp, "offset"), rep(0, nrow(nd)))
})

test_that("user-supplied `offset` accepts a per-row vector", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  offsets <- log(c(1, 5, 20))
  bp <- boot_predict(mod, newdata = nd, offset = offsets,
                     n_sim = 3000, seed = 1)
  gold <- exp(stats::coef(mod)[1] + stats::coef(mod)[2] * nd$x) * c(1, 5, 20)

  expect_equal(bp$boot_median, unname(gold), tolerance = 0.02)
  expect_equal(attr(bp, "offset"), offsets)
})

test_that("user-supplied `offset` warns when model has no offset", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x, family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5))

  expect_warning(
    boot_predict(mod, newdata = nd, offset = log(2), n_sim = 500, seed = 1),
    "model has no offset"
  )
})

test_that("user-supplied `offset` rejects non-numeric input", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  expect_error(
    boot_predict(mod, newdata = nd, offset = "log(1)", n_sim = 500),
    "must be numeric"
  )
})

# ---- Group 3: Auto-grid behavior with offsets ------------------------------

test_that("auto-grid excludes offset variables from the prediction grid", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)

  bp <- suppressWarnings(
    boot_predict(mod, at = list(x = c(-0.5, 0, 0.5)),
                 offset = log(10), n_sim = 500, seed = 1)
  )

  expect_false("effort" %in% names(bp))
  expect_true(all(c("x", "boot_mean", "boot_median") %in% names(bp)))
  expect_equal(nrow(bp), 3)
})

test_that("boot_predict warns when offset var missing from newdata", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5))  # no effort column

  expect_warning(
    boot_predict(mod, newdata = nd, n_sim = 500, seed = 1),
    "offset variable"
  )

  # And when it warns, predictions should equal exp(Xβ) at offset = 0
  bp <- suppressWarnings(
    boot_predict(mod, newdata = nd, n_sim = 3000, seed = 1)
  )
  gold_zero <- exp(stats::coef(mod)[1] + stats::coef(mod)[2] * nd$x)
  expect_equal(bp$boot_median, unname(gold_zero), tolerance = 0.02)
})

test_that("boot_predict does not silently fall through to globalenv for offset var", {
  # Regression test: earlier draft's eval() could pick up a stale 'effort'
  # in the caller's global env if newdata omitted the column. Verify that
  # the current implementation warns + defaults to 0 instead.
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5))

  # 'effort' exists in this test's environment — must not be used
  local_effort <- 999
  effort <- local_effort

  bp <- suppressWarnings(
    boot_predict(mod, newdata = nd, n_sim = 3000, seed = 1)
  )
  gold_zero <- exp(stats::coef(mod)[1] + stats::coef(mod)[2] * nd$x)
  expect_equal(bp$boot_median, unname(gold_zero), tolerance = 0.02)
})

# ---- Group 4: boot_median column ------------------------------------------

test_that("boot_predict output contains both boot_mean and boot_median", {
  dat <- make_offset_data()
  mod <- stats::glm(y ~ x, family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5))

  bp <- boot_predict(mod, newdata = nd, n_sim = 500, seed = 1)

  expect_true(all(c("boot_mean", "boot_median", "boot_se",
                    "boot_lwr", "boot_upr") %in% names(bp)))
  expect_type(bp$boot_median, "double")
  expect_length(bp$boot_median, nrow(nd))
})

test_that("boot_mean and boot_median agree closely on well-identified models", {
  dat <- make_offset_data(n = 1000)
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  bp <- boot_predict(mod, newdata = nd, n_sim = 3000, seed = 1)

  # Sampling-uncertainty Jensen effect is tiny at n=1000
  ratio <- bp$boot_mean / bp$boot_median
  expect_true(all(ratio > 0.99 & ratio < 1.01),
              info = paste("ratio range:", paste(round(range(ratio), 4), collapse = "-")))
})

test_that("boot_median tracks the plug-in point prediction", {
  dat <- make_offset_data(n = 1000)
  mod <- stats::glm(y ~ x + offset(log(effort)), family = poisson, data = dat)
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  gold_plugin <- stats::predict(mod, newdata = nd, type = "response")
  bp <- boot_predict(mod, newdata = nd, n_sim = 3000, seed = 1)

  # boot_median approximates exp(Xβhat) closely
  expect_equal(bp$boot_median, unname(gold_plugin), tolerance = 0.01)
})

# ---- Group 5: glmmTMB backend ---------------------------------------------

test_that("boot_predict handles offset with glmmTMB nbinom2 backend", {
  skip_if_not_installed("glmmTMB")

  dat <- make_offset_data(n = 500)
  mod <- glmmTMB::glmmTMB(
    y ~ x + offset(log(effort)),
    family = glmmTMB::nbinom2, data = dat
  )
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  gold <- stats::predict(mod, newdata = nd, type = "response")
  bp   <- boot_predict(mod, newdata = nd, n_sim = 2000, seed = 1)

  expect_equal(bp$boot_median, unname(gold), tolerance = 0.03)
  expect_true(all(c("boot_mean", "boot_median") %in% names(bp)))
})

test_that("boot_predict user-offset argument works with glmmTMB", {
  skip_if_not_installed("glmmTMB")

  dat <- make_offset_data(n = 500)
  mod <- glmmTMB::glmmTMB(
    y ~ x + offset(log(effort)),
    family = glmmTMB::nbinom2, data = dat
  )
  nd  <- data.frame(x = c(-0.5, 0, 0.5), effort = 10)

  bp <- boot_predict(mod, newdata = nd, offset = log(1),
                     n_sim = 2000, seed = 1)
  # Should NOT equal predict(mod, newdata=nd) — that's at effort=10, not 1
  gold_at_10 <- stats::predict(mod, newdata = nd, type = "response")
  # Should equal predictions at effort=1 (i.e., gold_at_10 / 10)
  expect_equal(bp$boot_median, unname(gold_at_10) / 10, tolerance = 0.03)
})
