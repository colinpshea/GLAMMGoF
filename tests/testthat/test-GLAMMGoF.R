# ══════════════════════════════════════════════════════════════════════════════
# testthat tests for GLAMMGoF
# Tests are designed to be fast -- heavy resampling functions use nReps = 3
# to verify structural correctness rather than statistical performance.
# Package data objects (countData, logitData, countModel_GLMM etc.) are used
# throughout to avoid model fitting overhead.
# ══════════════════════════════════════════════════════════════════════════════

library(testthat)
library(GLAMMGoF)

# ── Load package data ─────────────────────────────────────────────────────────
data("countData")
data("logitData")
data("countModel_GLMM")       # glmmTMB, 1 RE
data("countModel_GLMM2")      # glmmTMB, 2 RE
data("countModel_GLM")        # glmmTMB, no RE
data("logitModel_GLMM")       # glmmTMB, binomial, 1 RE

# ══════════════════════════════════════════════════════════════════════════════
# jensen_correct() tests
# ══════════════════════════════════════════════════════════════════════════════

test_that("jensen_correct returns scalar > 1 for glmmTMB with RE", {
  cf <- jensen_correct(countModel_GLMM)
  expect_type(cf, "double")
  expect_length(cf, 1)
  expect_gt(cf, 1)
})

test_that("jensen_correct returns 1 with warning for model with no REs", {
  expect_warning(cf <- jensen_correct(countModel_GLM))
  expect_equal(cf, 1)
})

test_that("jensen_correct 2 RE correction > 1 RE correction", {
  cf_1re <- jensen_correct(countModel_GLMM)
  cf_2re <- jensen_correct(countModel_GLMM2)
  expect_gt(cf_2re, cf_1re)
})

test_that("jensen_correct returns list when predictions supplied on response scale", {
  preds  <- predict(countModel_GLMM, type = "response", re.form = ~0)
  result <- jensen_correct(countModel_GLMM, predictions = preds,
                           scale = "response")
  expect_type(result, "list")
  expect_named(result, c("correction", "predictions_adjusted"))
  expect_equal(length(result$predictions_adjusted), length(preds))
  expect_true(all(result$predictions_adjusted >= preds))
})

test_that("jensen_correct se_adjusted = se * correction", {
  preds  <- predict(countModel_GLMM, type = "response", re.form = ~0)
  se     <- rep(0.5, length(preds))
  result <- jensen_correct(countModel_GLMM, predictions = preds,
                           se = se, scale = "response")
  expect_true("se_adjusted" %in% names(result))
  expect_equal(result$se_adjusted, se * result$correction)
})

test_that("jensen_correct lwr and upr adjusted correctly", {
  preds  <- predict(countModel_GLMM, type = "response", re.form = ~0)
  lwr    <- preds * 0.8
  upr    <- preds * 1.2
  result <- jensen_correct(countModel_GLMM, predictions = preds,
                           lwr = lwr, upr = upr, scale = "response")
  expect_true(all(c("lwr_adjusted", "upr_adjusted") %in% names(result)))
  expect_equal(result$lwr_adjusted, lwr * result$correction)
  expect_equal(result$upr_adjusted, upr * result$correction)
})

test_that("jensen_correct link scale backtransforms correctly", {
  preds_link <- predict(countModel_GLMM, type = "link", re.form = ~0)
  result     <- jensen_correct(countModel_GLMM, predictions = preds_link,
                               scale = "link")
  cf         <- jensen_correct(countModel_GLMM)
  expect_equal(result$predictions_adjusted, exp(preds_link) * cf)
})

test_that("jensen_correct errors for unsupported model class", {
  m_lm <- lm(y ~ Season, data = countData)
  expect_error(jensen_correct(m_lm))
})

# ══════════════════════════════════════════════════════════════════════════════
# bias_precision() tests
# ══════════════════════════════════════════════════════════════════════════════

test_that("bias_precision returns correct list structure", {
  out <- bias_precision(nReps = 3, testModel = countModel_GLMM,
                        testData = countData, DHARMaPlot = FALSE,
                        verbose = FALSE, seed = 42)
  expect_type(out, "list")
  expect_named(out, c("bias_precision_results",
                      "bias_precision_hist",
                      "bias_precision_summary"))
})

test_that("bias_precision results tibble has correct columns", {
  out <- bias_precision(nReps = 3, testModel = countModel_GLMM,
                        testData = countData, DHARMaPlot = FALSE,
                        verbose = FALSE, seed = 42)
  expect_named(out$bias_precision_results,
               c("simRep", "Group", "Metric", "value"))
})

test_that("bias_precision summary has correct metrics", {
  out <- bias_precision(nReps = 3, testModel = countModel_GLMM,
                        testData = countData, DHARMaPlot = FALSE,
                        verbose = FALSE, seed = 42)
  metrics <- unique(as.character(out$bias_precision_summary$Metric))
  expect_true(all(c("RRMSE", "RMAE", "RMedAE", "RBIAS") %in% metrics))
})

test_that("bias_precision summary has in-sample and out-of-sample groups", {
  out <- bias_precision(nReps = 3, testModel = countModel_GLMM,
                        testData = countData, DHARMaPlot = FALSE,
                        verbose = FALSE, seed = 42)
  groups <- unique(as.character(out$bias_precision_summary$Group))
  expect_true("In-sample performance" %in% groups)
  expect_true("Out-of-sample performance" %in% groups)
})

test_that("bias_precision works with bias_adjust = 'manual'", {
  expect_no_error(
    bias_precision(nReps = 3, testModel = countModel_GLMM,
                   testData = countData, DHARMaPlot = FALSE,
                   bias_adjust = "manual", verbose = FALSE, seed = 42)
  )
})

test_that("bias_precision conditional_predictions = TRUE fires message", {
  expect_message(
    bias_precision(nReps = 3, testModel = countModel_GLMM,
                   testData = countData, DHARMaPlot = FALSE,
                   conditional_predictions = TRUE, verbose = FALSE, seed = 42)
  )
})

test_that("bias_precision errors for bias_adjust = 'manual' + conditional_predictions", {
  expect_error(
    bias_precision(nReps = 3, testModel = countModel_GLMM,
                   testData = countData, DHARMaPlot = FALSE,
                   bias_adjust = "manual", conditional_predictions = TRUE,
                   verbose = FALSE, seed = 42)
  )
})

test_that("bias_precision works with bootstrap method", {
  expect_no_error(
    bias_precision(nReps = 3, testModel = countModel_GLMM,
                   testData = countData, DHARMaPlot = FALSE,
                   method = "bootstrap", verbose = FALSE, seed = 42)
  )
})

test_that("bias_precision manual correction moves RBIAS toward zero", {
  out_none <- bias_precision(nReps = 25, testModel = countModel_GLMM2,
                              testData = countData, DHARMaPlot = FALSE,
                              bias_adjust = "none", verbose = FALSE, seed = 42)
  out_manual <- bias_precision(nReps = 25, testModel = countModel_GLMM2,
                                testData = countData, DHARMaPlot = FALSE,
                                bias_adjust = "manual", verbose = FALSE, seed = 42)
  rbias_none <- out_none$bias_precision_summary$mn[
    out_none$bias_precision_summary$Metric == "RBIAS" &
    out_none$bias_precision_summary$Group == "In-sample performance"]
  rbias_manual <- out_manual$bias_precision_summary$mn[
    out_manual$bias_precision_summary$Metric == "RBIAS" &
    out_manual$bias_precision_summary$Group == "In-sample performance"]
  expect_lt(abs(rbias_manual), abs(rbias_none))
})

# ══════════════════════════════════════════════════════════════════════════════
# brier_auc() tests
# ══════════════════════════════════════════════════════════════════════════════

test_that("brier_auc returns correct list structure", {
  out <- brier_auc(nReps = 3, testModel = logitModel_GLMM,
                   testData = logitData, DHARMaPlot = FALSE, seed = 42)
  expect_type(out, "list")
  expect_true(length(out) > 0)
})

test_that("brier_auc summary contains AUC, Brier, and LogLoss", {
  out  <- brier_auc(nReps = 3, testModel = logitModel_GLMM,
                    testData = logitData, DHARMaPlot = FALSE, seed = 42)
  summ <- out[[grep("summary", names(out))]]
  metrics <- unique(as.character(summ$Metric))
  expect_true(any(c("AUC statistic", "Brier score", "Log loss") %in% metrics))
})

test_that("brier_auc errors for bias_adjust = 'manual'", {
  expect_error(
    brier_auc(nReps = 3, testModel = logitModel_GLMM,
              testData = logitData, DHARMaPlot = FALSE,
              bias_adjust = "manual", seed = 42)
  )
})

test_that("brier_auc conditional_predictions = TRUE fires message", {
  expect_message(
    brier_auc(nReps = 3, testModel = logitModel_GLMM,
              testData = logitData, DHARMaPlot = FALSE,
              conditional_predictions = TRUE, seed = 42)
  )
})

test_that("brier_auc works with bootstrap method", {
  expect_no_error(
    brier_auc(nReps = 3, testModel = logitModel_GLMM,
              testData = logitData, DHARMaPlot = FALSE,
              method = "bootstrap", seed = 42)
  )
})

test_that("brier_auc errors for non-binary response", {
  expect_error(
    brier_auc(nReps = 3, testModel = countModel_GLMM,
              testData = countData, DHARMaPlot = FALSE, seed = 42)
  )
})
