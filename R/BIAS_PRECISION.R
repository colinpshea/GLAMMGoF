#' Bootstrap or Monte Carlo assessment of RRMSE, RMAE, RMedAE, and RBIAS predictive performance statistics
#'
#' @description Assess in- and out-of-sample predictive performance of generalized linear and generalized additive models with continuous or integer response variables and with or without random effects and zero-inflation, using either repeated random holdout (Monte Carlo cross-validation) or bootstrap resampling with out-of-bag evaluation. Three performance statistics are reported: relative root mean squared error (RRMSE), calculated as sqrt(mean((observed - predicted)^2))/mean(observed)*100; relative mean absolute error (RMAE), calculated as mean(abs((observed - predicted)))/mean(observed)*100; relative median absolute error (RMedAE), calculated as median(abs((observed - predicted)))/mean(observed)*100; and relative bias (RBIAS), calculated as mean((observed - predicted))/mean(observed)*100. Note that all performance measures are based on population-level predictions (i.e., random effects are ignored). For models with a zero-inflation component, predictions account for zero-inflation (e.g., for glmmTMB, the predicted value represents the product of the mean_count and (1 - prob_zero)).
#'
#' All three accuracy statistics (`RRMSE`, `RMAE`, and `RMedAE`) express prediction error as a percentage of the mean observed value, which makes them interpretable on a common scale regardless of the units or magnitude of the response. To recover the corresponding raw error in the original units of the response, multiply any of these values by the observed mean and divide by 100. For example, an RRMSE of 18% with a sample mean of 50 implies a root mean squared error of 9 in the original units.
#'
#' `RRMSE (Relative Root Mean Squared Error)`: The square root of average squared prediction errors, expressed as a percentage of the mean observed value. Because errors are squared before averaging, `RRMSE` penalizes large individual errors more heavily than small ones, making it sensitive to cases where the model produces occasional extreme mispredictions.
#'
#' `RMAE (Relative Mean Absolute Error)`: The average absolute prediction error expressed as a percentage of the mean observed value. `RMAE` treats all errors proportionally to their size without the extra weight that squaring applies, and is often the most intuitive summary of typical prediction accuracy.
#'
#' `RMedAE (Relative Median Absolute Error)`: The median absolute prediction error expressed as a percentage of the mean observed value. Because it is based on the median rather than the mean of the error distribution, `RMedAE` is the most robust of the three to outlying predictions and gives the best picture of accuracy for a typical observation when the error distribution is skewed.
#'
#' Taken together, these three statistics tell a more complete story than any one alone. If `RRMSE` is notably larger than `RMAE`, that signals a handful of high-leverage mispredictions are inflating the squared-error average, which is a pattern `RMedAE` will often fail to reflect if the bulk of predictions are accurate. Conversely, close agreement among all three suggests errors are roughly symmetric and there are no extreme outliers driving the summary.
#'
#' `RBIAS (Relative Bias)`: The mean signed prediction error expressed as a percentage of the mean observed value, where positive values indicate systematic under-prediction and negative values indicate systematic over-prediction. `RBIAS` is independent of the accuracy metrics above: a model can be nearly unbiased on average yet still produce large errors, or it can be highly biased while still ranking observations correctly. Reporting `RBIAS` alongside `RRMSE` and `RMAE` therefore distinguishes random prediction noise from systematic directional error.
#'
#' @param nReps Desired number of bootstrap or Monte Carlo replicates. The default value is 100, but this number should be at least 1000 in practice.
#' @param testModel A regression model fit to testData in `glmmTMB` (with or without random effects), `glmer`/`glmer.nb`/`lmer' (with random effects), `glm`/`glm.nb`/`lm` (without random effects), or `gam` (with or without random effects). The response variable can be continuous or an integer, and possible statistical distributions include Poisson, negative binomial, gamma, tweedie, and gaussian.
#' @param testData A data frame with a continuous or integer response variable and continuous and/or categorical predictors.
#' @param propTrain The proportion of `testData` used for model-fitting and in-sample predictive performance when method = `holdout` (the default value is 0.8). The remaining proportion is used to assess out-of-sample predictive performance. This argument is ignored when method = `bootstrap`.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the `simulateResiduals()` function of the `DHARMa` package? The default is `TRUE`.
#' @param DHARMaReps You can also specify DHARMaReps if you want something other than the default of 1000 simulation replicates.
#' @param seed Optional integer seed for reproducibility. If `NULL` (the default), no seed is set and results will differ across runs.
#' @param method The resampling method to use. The default, `holdout`, repeatedly splits the data into random training and testing data sets (Monte Carlo cross-validation), whereas `bootstrap` samples the training data with replacement and evaluates in-sample performance on the bootstrap sample and out-of-sample performance on the out-of-bag observations not selected in the bootstrap sample (approximately 36.8% of observations on average). For well-behaved models and reasonably sized datasets, both methods should produce similar results; differences are most likely to emerge with small datasets, highly overdispersed data, or poorly specified models.
#' @note This function does not currently support binomial models with cbind() or proportion responses, and for binary 0/1 responses, use BRIER_AUC(). This function also supports models with spatial random effects (e.g, in glmmTMB), but it is much slower than for more conventional GLM(M)s and GAM(M)s.
#' @return This function returns four objects: a data frame with all of the bootstrapping or Monte Carlo resampling results (i.e., all `nReps` values for each performance statistic), a data frame with a summary (mean and 95% confidence intervals) of all replicates for each performance statistic, a histogram of values for each performance statistic, and a goodness-of-fit plot based on scaled residuals from the `simulateResiduals()` function of the `DHARMa` package. If `DHARMaPlot = FALSE`, then `simulateResiduals` is not used to assess the model residuals and only three of the four objects are returned.
#'
#' This package contains an example data set for fitting a negative binomial or Poisson regression called countData. Six example negative binomial regression model objects are also included: countModel_GLM is a GLM with no random effects; countModel_GLMM is a GLMM with one random effect; countModel_GLMM2 is a GLMM with two random effects; countModel_GAM is a GAM with no random effects; countModel_GAMM is a GAMM with one random effect; and countModel_GAMM2 is a GAMM with two random effects. GLMs and GLMMs were fitted using `glmmTMB`, whereas GAMs and GAMMs were fitted using `mgcv`:
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
#' Bootstrapping or Monte Carlo resampling of the performance statistics requires specifying the data and model being tested, the desired number of replicates (the default is 100 but should be at least 1000 in practice), the resampling method `holdout` or `bootstrap`, the proportion of data used for training when `method = "holdout"` (the default is 0.8), whether to use `DHARMa` residual diagnostics (the default is `TRUE`), the number of `DHARMa` simulation replicates (the default is 1000), and an optional integer seed for reproducibility:
#'
#' BIAS_PRECISION(nReps = 100, testModel = countModel_GLMM, testData = countData, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000, seed = 123, method = "holdout")
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate bind_rows
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_grid theme_bw theme element_blank element_line element_text labs unit scale_y_continuous
#' @importFrom DHARMa simulateResiduals
#' @importFrom glmmTMB ranef glmmTMB
#' @importFrom lme4 glmer lmer glmer.nb
#' @importFrom MASS glm.nb
#' @importFrom mgcv gam predict.gam
#' @export
BIAS_PRECISION <- function(nReps = 100, testModel = NULL, testData = NULL,
                             propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000,
                             seed = NULL, method = c("holdout", "bootstrap")) {

  # --- specify bootstrapping method
  method = match.arg(method)

  # --- Optional seed ---
  if (!is.null(seed)) set.seed(seed)

  # --- Cost functions ---
  fit_cost_rrmse <- function(y, yhat) sqrt(mean((y - yhat)^2)) / mean(y) * 100
  fit_cost_rmedae  <- function(y, yhat) median(abs(y - yhat)) / mean(y) * 100
  fit_cost_rmae  <- function(y, yhat) mean(abs(y - yhat)) / mean(y) * 100
  fit_cost_rbias <- function(y, yhat) mean(y - yhat) / mean(y) * 100

  # --- Validate inputs ---
  stopifnot("testModel cannot be NULL" = !is.null(testModel))
  stopifnot("testData cannot be NULL"  = !is.null(testData))
  stopifnot("propTrain must be between 0 and 1" = propTrain > 0 && propTrain < 1)

  resp_var  <- all.vars(formula(testModel))[1]
  is_binary <- function(x) length(unique(x)) == 2 && all(x %in% c(0, 1))
  stopifnot("Response variable is binary! Use BRIER_AUC() instead" =
              !is_binary(testData[[resp_var]]))

  # Check for unsupported binomial response types
  resp_expr <- formula(testModel)[[2]]
  is_cbind  <- is.call(resp_expr) && deparse(resp_expr[[1]]) == "cbind"
  is_prop   <- is.call(resp_expr) && grepl("/", deparse(resp_expr))
  stopifnot(
    "cbind() and proportion binomial responses are not supported. See ?BIAS_PRECISION for details." =
      !is_cbind && !is_prop
  )

  # --- Pre-compute model class flags (once, outside loop) ---
  mc         <- class(testModel)
  is_glmmTMB <- "glmmTMB"  %in% mc
  is_gam     <- "gam"      %in% mc
  is_glmer   <- "glmerMod" %in% mc
  is_lmer    <- "lmerMod"  %in% mc
  is_negbin  <- "negbin"   %in% mc
  is_glm     <- "glm"      %in% mc && !is_gam
  is_lm      <- "lm"       %in% mc && !is_gam && !is_glm  # fixed: was !is_glm && !is_glm

  # --- Pre-compute GAM RE metadata (once, outside loop) ---
  # smooth$label passed to exclude= (e.g. "s(Site)"); smooth$term is the column name.
  # Using smooth object fields directly avoids any regex on the label string.
  gam_re_labels <- NULL
  gam_re_terms  <- NULL
  if (is_gam) {
    re_smooths <- Filter(function(s) isTRUE(s$random), testModel$smooth)
    if (length(re_smooths) > 0) {
      gam_re_labels <- sapply(re_smooths, function(s) s$label)  # vector: one per RE smooth
      gam_re_terms  <- sapply(re_smooths, function(s) s$term)   # retained for metadata
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

  # --- Predict helper: marginal/population-level predictions, RE excluded ---
  # glmmTMB: allow.new.levels = TRUE handles RE groups absent from training data
  # GAM:     exclude= + newdata.guaranteed=TRUE drops all RE smooths cleanly;
  #          no column deletion needed, full newdata passed untouched
  # lme4:    re.form = ~0 suppresses all random effects
  get_preds <- function(m, newdata) {
    if (is_glmmTMB) {
      predict(m, type = "response", newdata = newdata,
              allow.new.levels = TRUE)
    } else if (is_gam) {
      if (!is.null(gam_re_labels)) {
        predict(m, type = "response",
                exclude            = gam_re_labels,  # vector: all RE smooth labels excluded
                newdata            = newdata,
                newdata.guaranteed = TRUE)            # mgcv won't look for RE columns in newdata
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
    if (is.null(m_train)) next  # skip failed fits cleanly; no stale model carried forward

    y_train    <- train[[resp_var]]
    y_test     <- test[[resp_var]]
    yhat_train <- get_preds(m_train, train)
    yhat_test  <- get_preds(m_train, test)

    results[[j]] <- data.frame(
      train_RRMSE = fit_cost_rrmse(y_train, yhat_train),
      test_RRMSE  = fit_cost_rrmse(y_test,  yhat_test),
      train_RMedAE  = fit_cost_rmedae(y_train,  yhat_train),
      test_RMedAE   = fit_cost_rmedae(y_test,   yhat_test),
      train_RMAE  = fit_cost_rmae(y_train,  yhat_train),
      test_RMAE   = fit_cost_rmae(y_test,   yhat_test),
      train_RBIAS = fit_cost_rbias(y_train, yhat_train),
      test_RBIAS  = fit_cost_rbias(y_test,  yhat_test)
    )
  }

  # Guard against all replicates failing
  results_clean <- results[!vapply(results, is.null, logical(1))]
  if (length(results_clean) == 0)
    stop("All model fits failed - no results to summarise.")

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
    return(list(bias_precision_results  = results_df,
                bias_precision_hist     = results_plot,
                bias_precision_summary  = results_summary,
                dharmaPlot          = dharmaPlot))
  }

  list(bias_precision_results = results_df,
       bias_precision_hist    = results_plot,
       bias_precision_summary = results_summary)
}
