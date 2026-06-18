#' Bootstrap or Monte Carlo assessment of RRMSE, RMAE, RMedAE, and RBIAS predictive performance statistics
#'
#' @description Assess in- and out-of-sample predictive performance of generalized linear and generalized additive models with continuous or integer response variables and with or without random effects and zero-inflation, using either repeated random holdout (Monte Carlo cross-validation) or bootstrap resampling with out-of-bag evaluation. Four performance statistics are reported: relative root mean squared error (RRMSE), relative Mean Absolute Error (RMAE), relative median absolute error (RMedAE), and relative bias (RBIAS). Note that all performance measures are based on population-level predictions (i.e., random effects are ignored). For models with a zero-inflation component, predictions account for zero-inflation (e.g., for glmmTMB, the predicted value represents the product of the mean_count and (1 - prob_zero)).
#'
#' All three accuracy statistics (`RRMSE`, `RMAE`, and `RMedAE`) express prediction error as a percentage of the mean observed value, which makes them interpretable on a common scale regardless of the units or magnitude of the response. To recover the corresponding raw error in the original units of the response, multiply any of these values by the observed mean and divide by 100. For example, an RRMSE of 18% with a sample mean of 50 implies a root mean squared error of 9 in the original units.
#'
#' `RRMSE (Relative Root Mean Squared Error)`: The square root of average squared prediction errors, expressed as a percentage of the mean observed value. Because errors are squared before averaging, `RRMSE` penalizes large individual errors more heavily than small ones, making it sensitive to cases where the model produces occasional extreme mispredictions.
#'
#'  `RRMSE = sqrt(mean((predicted - observed)^2))/mean(observed)*100`
#'
#' `RMAE (Relative Mean Absolute Error)`: The average absolute prediction error expressed as a percentage of the mean observed value. `RMAE` treats all errors proportionally to their size without the extra weight that squaring applies, and is often the most intuitive summary of typical prediction accuracy.
#'
#'  `RMAE = mean(abs((predicted - observed)))/mean(observed)*100`
#'
#' `RMedAE (Relative Median Absolute Error)`: The median absolute prediction error expressed as a percentage of the mean observed value. Because it is based on the median rather than the mean of the error distribution, `RMedAE` is the most robust of the three to outlying predictions and gives the best picture of accuracy for a typical observation when the error distribution is skewed.
#'
#'  `RMedAE = median(abs((predicted - observed)))/mean(observed)*100`
#'
#' Taken together, these three statistics tell a more complete story than any one alone. If RRMSE is notably larger than `RMAE`, that signals a handful of high-leverage mispredictions are inflating the squared-error average, which is a pattern `RMedAE` will often fail to reflect if the bulk of predictions are accurate. Conversely, close agreement among all three suggests errors are roughly symmetric and there are no extreme outliers driving the summary.
#'
#' `RBIAS (Relative Bias)`: The mean signed prediction error expressed as a percentage of the mean observed value, where positive values indicate systematic over-prediction and negative values indicate systematic under-prediction. RBIAS is independent of the accuracy metrics above: a model can be nearly unbiased on average yet still produce large errors, or it can be highly biased while still ranking observations correctly. Reporting RBIAS alongside RRMSE and RMAE therefore distinguishes random prediction noise from systematic directional error.
#'
#'  `RBIAS = mean((predicted - observed))/mean(observed)*100`
#'
#' @param nReps Desired number of bootstrap or Monte Carlo replicates. The default value is 100, but this number should be at least 1000 in practice.
#' @param testModel A regression model fit to testData in glmmTMB (with or without random effects), glmer/glmer.nb/lmer (with random effects), glm/glm.nb/lm (without random effects), or gam (with or without random effects). The response variable can be continuous or an integer, and possible statistical distributions include Poisson, negative binomial, gamma, tweedie, and gaussian.
#' @param testData A data frame with a continuous or integer response variable and continuous and/or categorical predictors.
#' @param propTrain The proportion of `testData` used for model-fitting and in-sample predictive performance when method = `holdout` (the default value is 0.8). The remaining proportion is used to assess out-of-sample predictive performance. This argument is ignored when method = `bootstrap`.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the simulateResiduals() function of the DHARMa package? The default is TRUE.
#' @param DHARMaReps You can also specify DHARMaReps if you want something other than the default of 1000 simulation replicates.
#' @param testZI Logical. If `TRUE` and `DHARMaPlot = TRUE`, runs `testZeroInflation` on the simulated residuals. Most relevant for count models (Poisson, negative binomial, ZIP, hurdle). Default is `TRUE`.
#' @param seed Optional integer seed for reproducibility. If NULL (the default), no seed is set and results will differ across runs.
#' @param method The resampling method to use. The default, `holdout`, repeatedly splits the data into random training and testing data sets (Monte Carlo cross-validation), whereas `bootstrap` samples the training data with replacement and evaluates in-sample performance on the bootstrap sample and out-of-sample performance on the out-of-bag observations not selected in the bootstrap sample (approximately 36.8% of observations on average). For well-behaved models and reasonably sized datasets, both methods should produce similar results; differences are most likely to emerge with small datasets, highly overdispersed data, or poorly specified models.
#' @param bias_adjust Logical. If `TRUE` and the model is a `glmmTMB` fit, applies a bias correction to marginal predictions on the response scale via `do.bias.correct = TRUE` in `predict.glmmTMB`. This corrects for the systematic underestimation of the arithmetic mean that arises when backtransforming population-level predictions through a nonlinear link function in the presence of random effects (Jensen's inequality). The default is `FALSE`, which preserves the uncorrected marginal prediction behavior and allows `RBIAS` to reflect the full magnitude of RE-driven prediction bias. This argument is silently ignored for non-`glmmTMB` models. See the note below and Thorson & Kristensen (2016) for details.
#' @note This function does not currently support binomial models with cbind() or proportion responses, and for binary 0/1 responses, use brier_auc(). This function also supports models with spatial random effects (e.g, in glmmTMB), but it is much slower than for more conventional GLM(M)s and GAM(M)s.
#'
#' **Random effects and Jensen's inequality:** All predictions are population-level (i.e., random effects are set to zero via `re.form = ~0`). For models with a nonlinear link function (e.g., log, logit) and random effects, backtransforming the linear predictor to the response scale introduces a systematic negative bias in the predicted arithmetic mean. This occurs because Jensen's inequality implies that `E[exp(eta)] > exp(E[eta])` for any random variable `eta`, so `exp(beta_0)` underestimates the true mean by a factor of approximately `exp(sigma^2 / 2)`, where `sigma^2` is the total random effect variance (summed across all RE terms). The magnitude of this bias grows rapidly with RE variance: a model with two random effects of modest size (e.g., SD = 0.3 and 0.4) can produce marginal predictions that underestimate observed values by 10% or more. Consistent negative `RBIAS` in the output of this function -- particularly when both in-sample and out-of-sample values are negative -- may therefore reflect this structural property of the model rather than misspecification. To obtain bias-corrected marginal predictions for `glmmTMB` models, set `bias_adjust = TRUE`. Note that this correction is not available for `lme4` or `mgcv` models; however, `lme4` and `mgcv` GAMMs already use marginal predictions that exclude random effects, and the Jensen bias is present in those models regardless.
#' @references Hyndman, R.J. and Koehler, A.B. (2006) Another look at measures of forecast accuracy. \emph{International Journal of Forecasting}, 22, 679--688.
#' @references Thorson, J.T. and Kristensen, K. (2016) Implementing a generic method for bias correction in statistical models using random effects, with spatial and population dynamics examples. \emph{Fisheries Research}, 175, 66--74.
#' @return This function returns either three, four, or five objects depending on the values of `DHARMaPlot` and `testZI`: a data frame with all bootstrapping or Monte Carlo resampling results (i.e., all `nReps` values for each performance statistic), a data frame with a summary (mean and 95% confidence intervals) of all replicates for each performance statistic, and a histogram of values for each performance statistic. If `DHARMaPlot = TRUE`, a fourth object is also returned: a goodness-of-fit plot based on scaled residuals from `simulateResiduals()`. If `testZI = TRUE` and `DHARMaPlot = TRUE`, a fifth object `dharmaZI` is also returned containing the result of `testZeroInflation()`. In the histogram, a blue dotted vertical line indicates the mean across replicates.
#'
#' This package contains an example data set for fitting a negative binomial or Poisson regression called countData. Six example negative binomial regression model objects are also included: countModel_GLM is a GLM with no random effects; countModel_GLMM is a GLMM with one random effect; countModel_GLMM2 is a GLMM with two random effects; countModel_GAM is a GAM with no random effects; countModel_GAMM is a GAMM with one random effect; and countModel_GAMM2 is a GAMM with two random effects. GLMs and GLMMs were fitted using glmmTMB, whereas GAMs and GAMMs were fitted using mgcv:
#'
#' countModel_GLM <- glmmTMB(y ~ Season + Temp, family = nbinom2, data = countData)
#'
#' countModel_GLMM <- glmmTMB(y ~ Season + Temp + (1|Site), family = nbinom2, data = countData)
#'
#' countModel_GLMM2 <- glmmTMB(y ~ Season + Temp + (1|Site) + (1|Year), family = nbinom2, data = countData)
#'
#' countModel_GAM <- gam(y ~ Season + s(Temp), family = nb, data = countData)
#'
#' countModel_GAMM <- gam(y ~ Season + s(Temp) + s(Site, bs = "re"), family = nb, data = countData)
#'
#' countModel_GAMM2 <- gam(y ~ Season + s(Temp) + s(Site, bs = "re") + s(Year, bs = "re"), family = nb, data = countData)
#'
#' Bootstrapping or Monte Carlo resampling of the performance statistics requires specifying the data and model being tested, the desired number of replicates (the default is 100 but should be at least 1000 in practice),  the proportion of data used for training when `method = "holdout"` (the default is 0.8), whether to use `DHARMa` residual diagnostics (the default is `TRUE`), whether to use `DHARMa` to test for zero-inflation (the default is `TRUE`), the number of `DHARMa` simulation replicates (the default is 1000), and an optional integer seed for reproducibility, and the resampling method `holdout` or `bootstrap`:
#'
#' bias_precision(nReps = 100, testModel = countModel_GLMM, testData = countData, propTrain = 0.8, DHARMaPlot = TRUE, testZI = TRUE, DHARMaReps = 1000, seed = 123, method = "holdout")
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate bind_rows
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_grid theme_bw theme element_blank element_line element_text labs unit scale_y_continuous
#' @importFrom DHARMa simulateResiduals testZeroInflation
#' @importFrom glmmTMB ranef glmmTMB
#' @importFrom lme4 glmer lmer glmer.nb
#' @importFrom MASS glm.nb
#' @importFrom mgcv gam predict.gam
#' @importFrom stats complete.cases formula
#' @export
bias_precision <- function(nReps = 100, testModel = NULL, testData = NULL,
                           propTrain = 0.8, DHARMaPlot = TRUE, testZI = TRUE, DHARMaReps = 1000,
                           seed = NULL, method = c("holdout", "bootstrap"),
                           bias_adjust = FALSE) {

  # --- specify bootstrapping method
  method = match.arg(method)

  # --- Optional seed ---
  if (!is.null(seed)) set.seed(seed)

  # --- Cost functions ---
  fit_cost_rrmse <- function(y, yhat) sqrt(mean((yhat - y)^2)) / mean(y) * 100
  fit_cost_rmedae  <- function(y, yhat) median(abs(yhat - y)) / mean(y) * 100
  fit_cost_rmae  <- function(y, yhat) mean(abs(yhat - y)) / mean(y) * 100
  fit_cost_rbias <- function(y, yhat) mean(yhat - y) / mean(y) * 100

  # --- Validate inputs ---
  stopifnot("testModel cannot be NULL" = !is.null(testModel))
  stopifnot("testData cannot be NULL"  = !is.null(testData))
  stopifnot("propTrain must be between 0 and 1" = propTrain > 0 && propTrain < 1)

  resp_var  <- all.vars(formula(testModel))[1]
  is_binary <- function(x) length(unique(x)) == 2 && all(x %in% c(0, 1))
  stopifnot("Response variable is binary! Use brier_auc() instead" =
              !is_binary(testData[[resp_var]]))

  # Check for unsupported binomial response types
  resp_expr <- formula(testModel)[[2]]
  is_cbind  <- is.call(resp_expr) && deparse(resp_expr[[1]]) == "cbind"
  is_prop   <- is.call(resp_expr) && grepl("/", deparse(resp_expr))
  stopifnot(
    "cbind() and proportion binomial responses are not supported. See ?bias_precision for details." =
      !is_cbind && !is_prop
  )

  # Remove rows with NA in any model variable (response or covariates)
  model_vars    <- all.vars(formula(testModel))
  model_vars    <- model_vars[model_vars %in% names(testData)]
  n_before      <- nrow(testData)
  complete_rows <- complete.cases(testData[, model_vars, drop = FALSE])
  testData      <- testData[complete_rows, ]
  n_dropped     <- n_before - nrow(testData)
  if (n_dropped == 1) warning(n_dropped, " row with NA values in model variables (response or covariates) was removed before resampling.")
  if (n_dropped > 1) warning(n_dropped, " rows with NA values in model variables (response or covariates) were removed before resampling.")

  # --- Pre-compute model class flags (once, outside loop) ---
  mc         <- class(testModel)
  is_glmmTMB <- "glmmTMB"  %in% mc
  is_gam     <- "gam"      %in% mc
  is_glmer   <- "glmerMod" %in% mc
  is_lmer    <- "lmerMod"  %in% mc
  is_negbin  <- "negbin"   %in% mc
  is_glm     <- "glm"      %in% mc && !is_gam
  is_lm      <- "lm"       %in% mc && !is_gam && !is_glm

  # --- Pre-compute GAM RE metadata (once, outside loop) ---
  gam_re_labels <- NULL
  gam_re_terms  <- NULL
  if (is_gam) {
    re_smooths <- Filter(function(s) isTRUE(s$random), testModel$smooth)
    if (length(re_smooths) > 0) {
      gam_re_labels <- sapply(re_smooths, function(s) s$label)
      gam_re_terms  <- sapply(re_smooths, function(s) s$term)
    }
  }

  # --- Fit helper: dispatch on model class, return NULL on failure ---
  fit_model <- function(train) {
    tryCatch({
      if (is_glmmTMB) {
        glmmTMB(formula(testModel, component = "cond"),
                family      = family(testModel),
                dispformula = formula(testModel, component = "disp"),
                ziformula   = formula(testModel, component = "zi"),
                data        = train)
      } else if (is_gam) {
        gam(formula(testModel), family = family(testModel), data = train)
      } else if (is_glmer) {
        if (grepl("Negative Binomial", family(testModel)$family))
          glmer.nb(formula(testModel), data = train)
        else
          glmer(formula(testModel), family = family(testModel), data = train)
      } else if (is_lmer) {
        lmer(formula(testModel), data = train)
      } else if (is_negbin) {
        glm.nb(formula(testModel), data = train)
      } else if (is_glm) {
        glm(formula(testModel), family = family(testModel), data = train)
      } else if (is_lm) {
        lm(formula(testModel), data = train)
      }
    }, error = function(e) {
      message("Model failed on replicate: ", e$message)
      NULL
    })
  }

  # --- Predict helper ---
  get_preds <- function(m, newdata) {
    if (is_glmmTMB) {
      preds <- predict(m, type = "response", newdata = newdata, re.form = ~0,
                       allow.new.levels = TRUE, do.bias.correct = bias_adjust)
      # do.bias.correct = TRUE returns a matrix with columns Estimate, Std. Error,
      # Est. (bias.correct), and Std. (bias.correct); extract the bias-corrected point estimate
      if (is.matrix(preds)) preds[, "Est. (bias.correct)"] else preds
    } else if (is_gam) {
      if (!is.null(gam_re_labels)) {
        predict(m, type = "response",
                exclude            = gam_re_labels,
                newdata            = newdata,
                newdata.guaranteed = TRUE)
      } else {
        predict(m, type = "response", newdata = newdata)
      }
    } else if (is_glmer || is_lmer) {
      predict(m, type = "response", re.form = ~0,
              allow.new.levels = TRUE, newdata = newdata)
    } else {
      predict(m, type = "response", newdata = newdata)
    }
  }

  # --- Bootstrap loop ---
  results <- vector("list", nReps)

  for (j in seq_len(nReps)) {
    if (method == "holdout"){
      train_idx <- sample(seq_len(nrow(testData)), size = floor(propTrain * nrow(testData)))
      train <- testData[ train_idx, ]
      test  <- testData[-train_idx, ]
    } else {
      train_idx <- sample(seq_len(nrow(testData)), size = nrow(testData), replace = TRUE)
      train <- testData[ train_idx, ]
      test  <- testData[setdiff(seq_len(nrow(testData)), unique(train_idx)), ]
    }

    m_train <- fit_model(train)
    if (is.null(m_train)) next

    y_train    <- train[[resp_var]]
    y_test     <- test[[resp_var]]
    yhat_train <- get_preds(m_train, train)
    yhat_test  <- get_preds(m_train, test)

    results[[j]] <- data.frame(
      train_RRMSE  = fit_cost_rrmse( y_train, yhat_train),
      test_RRMSE   = fit_cost_rrmse( y_test,  yhat_test),
      train_RMedAE = fit_cost_rmedae(y_train, yhat_train),
      test_RMedAE  = fit_cost_rmedae(y_test,  yhat_test),
      train_RMAE   = fit_cost_rmae(  y_train, yhat_train),
      test_RMAE    = fit_cost_rmae(  y_test,  yhat_test),
      train_RBIAS  = fit_cost_rbias( y_train, yhat_train),
      test_RBIAS   = fit_cost_rbias( y_test,  yhat_test)
    )
  }

  # Guard against all replicates failing
  results_clean <- results[!vapply(results, is.null, logical(1))]
  if (length(results_clean) == 0)
    stop("All model fits failed - no results to summarize.")

  n_failed <- nReps - length(results_clean)
  if (n_failed > 0)
    message(n_failed, " of ", nReps, " bootstrap replicates failed (", round(n_failed / nReps*100, 1), "%). If this percentage is high, consider reviewing your model structure or increasing propTrain to use a larger proportion of data for model-fitting.")

  # --- Tidy results ---
  results_df <- bind_rows(results_clean, .id = "simRep") %>%
    pivot_longer(cols = -simRep, names_to = "metric", values_to = "value") %>%
    separate(metric, into = c("Group", "Metric")) %>%
    mutate(
      Group  = factor(Group, levels = c("train", "test"),
                      labels = c("In-sample performance", "Out-of-sample performance")),
      Metric = factor(Metric, levels = c("RRMSE", "RMAE", "RMedAE", "RBIAS"))
    )

  # ---- count up, report, and omit results with NA
  n_na <- sum(is.na(results_df$value) | is.infinite(results_df$value))

  if (n_na > 0) cat(n_na, " NA or Inf values removed from ", nrow(results_df)," total bootstrap observations. ", "This may indicate model instability or sparse data.")

  results_df <- results_df[!is.na(results_df$value) & !is.infinite(results_df$value), ]

  results_summary <- results_df %>%
    group_by(Group, Metric) %>%
    summarise(mn    = mean(value),
              lwr95 = quantile(value, 0.025),
              upr95 = quantile(value, 0.975),
              .groups = "drop")

  # --- Jensen's inequality diagnostic message ---
  # Triggered when both in- and out-of-sample RBIAS are consistently negative,
  # which is the signature of marginal prediction bias from strong random effects
  # on a nonlinear link scale. Only fires when bias_adjust = FALSE to avoid
  # misleading the user when the correction has already been applied.
  if (!bias_adjust && is_glmmTMB) {
    rbias_summary <- results_summary[results_summary$Metric == "RBIAS", ]
    both_negative <- all(rbias_summary$mn < -10)
    if (both_negative) {
      message(
        "Note: Both in-sample and out-of-sample RBIAS are consistently negative (< -10%). ",
        "This may indicate that strong random effects are causing marginal predictions to ",
        "underestimate the arithmetic mean on the response scale (Jensen's inequality): ",
        "exp(beta) underestimates the true mean when random effect variance is large. ",
        "Consider re-running with bias_adjust = TRUE to apply a lognormal bias correction. ",
        "See ?bias_precision for details."
      )
    }
  }

  results_plot <- ggplot(results_df, aes(x = value)) +
    geom_histogram(color = "black", fill = "grey") +
    facet_grid(Group ~ Metric, scales = "free") +
    geom_vline(data = results_summary, aes(xintercept = mn), color = "blue", linetype = "dotted", linewidth = 0.8) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey90", linetype = "solid"),
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
          axis.text          = element_text(colour = "black"),
          panel.spacing      = unit(1.5, "lines")) +
    labs(x = "% relative to true mean", y = "Frequency") +
    scale_y_continuous(expand = expansion(mult = c(0,0.01)))

  if (DHARMaPlot) {
    dharmaPlot <- tryCatch(
      withCallingHandlers(
        simulateResiduals(n = DHARMaReps, testModel, plot = TRUE),
        warning = function(w) {
          message("DHARMa warning (", class(testModel)[1], "): ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        message("DHARMa failed (", class(testModel)[1], "): ", conditionMessage(e))
        NULL
      }
    )
    dharmaZI <- if (testZI && !is.null(dharmaPlot)) {
      tryCatch(
        withCallingHandlers(
          testZeroInflation(dharmaPlot),
          warning = function(w) {
            message("DHARMa ZI warning (", class(testModel)[1], "): ", conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(e) {
          message("DHARMa ZI test failed (", class(testModel)[1], "): ", conditionMessage(e))
          NULL
        }
      )
    } else NULL
    return(list(bias_precision_results  = results_df,
                bias_precision_hist     = results_plot,
                bias_precision_summary  = results_summary,
                dharmaPlot             = dharmaPlot,
                dharmaZI               = dharmaZI))
  } else {
    return(list(bias_precision_results  = results_df,
                bias_precision_hist     = results_plot,
                bias_precision_summary  = results_summary))
  }
}
