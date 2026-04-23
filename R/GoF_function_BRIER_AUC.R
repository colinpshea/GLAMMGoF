#' Bootstrap or Monte Carlo assessment of RRMSE, RMAD, and RBIAS predictive performance statistics
#'
#' @description Assess in- and out-of-sample predictive performance of generalized linear and generalized additive models with binary response variables and with or without random effects, using either repeated random holdout (Monte Carlo cross-validation) or bootstrap resampling with out-of-bag evaluation. Two performance statistics are reported: Brier scores (see the `rms` package documentation for details), which range from 0 to 1 with values closer to 0 indicating a better-predicting model and where sqrt(Brier score) is the average difference between the predicted probability and the observed value (0 or 1); and AUC, an aggregated metric that evaluates how well a model classifies positive and negative outcomes at all possible probability cutoffs, ranging from 0 to 1 with values closer to 1 indicating a better classifier and where an AUC of 0.5 suggests performance no better than random guessing. Note that all performance measures are based on population-level predictions (i.e., random effects are ignored when present).
#' @param nReps Desired number of bootstrap or Monte Carlo replicates. The default value is 100, but this number should be at least 1000 in practice.
#' @param testModel A logistic regression model fitted to testData using `glmmTMB` (with or without random effects), `glmer` (with random effects), `glm` (without random effects), or `gam` (with or without random effects).
#' @param testData A data frame with a binary response variable and continuous and/or categorical predictor variables.
#' @param propTrain The proportion of `testData` used for model-fitting and in-sample predictive performance when method = `holdout` (the default value is 0.8). The remaining proportion is used to assess out-of-sample predictive performance. This argument is ignored when method = `bootstrap`.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the `simulateResiduals()` function of the `DHARMa` package? The default is `TRUE`.
#' @param DHARMaReps If DHARMaPlot is `TRUE`, you can also specify DHARMaReps if you want something other than the default of 1000 simulation replicates.
#' @param seed Optional integer seed for reproducibility. If `NULL` (the default), no seed is set and results will differ across runs.
#' @param method The resampling method to use. The default, `holdout`, repeatedly splits the data into random training and testing data sets (Monte Carlo cross-validation), whereas `bootstrap` samples the training data with replacement and evaluates in-sample performance on the bootstrap sample and out-of-sample performance on the out-of-bag observations not selected in the bootstrap sample (approximately 36.8% of observations on average). For well-behaved models and reasonably sized datasets, both methods should produce similar results; differences are most likely to emerge with small datasets, highly overdispersed data, or poorly specified models.
#' @note This function only supports binary 0/1 responses and does not currently support binomial models with cbind() or proportion responses. This function also supports models with spatial random effects  (e.g, in glmmTMB), but it is much slower than for more conventional GLM(M)s and GAM(M)s
#' @return This function returns four objects: a data frame with all of the bootstrapping or Monte Carlo resampling results (i.e., all `nReps` values for each performance statistic), a data frame with a summary (mean and 95% confidence intervals) of all replicates for each performance statistic, a histogram of values for each performance statistic, and a goodness-of-fit plot based on scaled residuals from the `simulateResiduals()` function of the `DHARMa` package. If `DHARMaPlot = FALSE`, then `simulateResiduals` is not used to assess the model residuals and only three of the four objects are returned.
#'
#' This package contains an example data set to fit a logistic regression called logitData. Six example logistic regression model objects are also included: logitModel_GLM is a GLM with no random effects; logitModel_GLMM is a GLMM with one random effect; logitModel_GLMM2 is a GLMM with two random effects; logitModel_GAM is a GAM with no random effects; logitModel_GAMM is a GAMM with one random effect; and logitModel_GAMM2 is a GAMM with two random effects. GLMs and GLMMs were fitted using `glmmTMB`, whereas GAMs and GAMMs were fitted using `mgcv`:
#'
#' logitModel_GLM <- glmmTMB(y ~ Season + Temp, family = binomial, data = logitData)
#'
#' logitModel_GLMM <- glmmTMB(y ~ Season + Temp + (1|Site), family = binomial, data = logitData)
#'
#' logitModel_GLMM2 <- glmmTMB(y ~ Season + Temp + (1|Site) + (1|Year), family = binomial, data = logitData)
#'
#' logitModel_GAM <- gam(y ~ Season + s(Temp), family = binomial, data = logitData)
#'
#' logitModel_GAMM <- gam(y ~ Season + s(Temp) + s(Site, bs = "re"), family = binomial, data = logitData)
#'
#' logitModel_GAMM2 <- gam(y ~ Season + s(Temp) + s(Site, bs = "re") + s(Year, bs = "re"), family = binomial, data = logitData)
#'
#' Bootstrapping or Monte Carlo resampling of the performance statistics requires specifying the data and model being tested, the desired number of replicates (the default is 100 but should be at least 1000 in practice), the resampling method `holdout` or `bootstrap`, the proportion of data used for training when `method = "holdout"` (the default is 0.8), whether to use `DHARMa` residual diagnostics (the default is `TRUE`), the number of `DHARMa` simulation replicates (the default is 1000), and an optional integer seed for reproducibility:
#'
#' BRIER_AUC(nReps = 100, testModel = logitModel_GLMM, testData = logitData, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000, seed = 42, method = "holdout")
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate bind_rows
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_histogram facet_grid theme_bw theme element_blank element_line element_text labs unit scale_x_continuous expansion
#' @importFrom DHARMa simulateResiduals
#' @importFrom rms val.prob
#' @importFrom glmmTMB glmmTMB
#' @importFrom lme4 glmer
#' @importFrom mgcv gam predict.gam
#' @export
BRIER_AUC <- function(nReps = 100, testModel = NULL, testData = NULL,
                      propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000,
                      seed = NULL, method = c("holdout", "bootstrap")) {

  # --- specify bootstrapping method
  method = match.arg(method)

  # --- Optional seed ---
  if (!is.null(seed)) set.seed(seed)

  # --- Validate inputs ---
  stopifnot("testModel cannot be NULL" = !is.null(testModel))
  stopifnot("testData cannot be NULL"  = !is.null(testData))
  stopifnot("propTrain must be between 0 and 1" = propTrain > 0 && propTrain < 1)

  resp_var  <- all.vars(formula(testModel))[1]
  is_binary <- function(x) length(unique(x)) == 2 && all(x %in% c(0, 1))
  stopifnot("Response variable is not binary! Use RRMSE_RMAD_RBIAS() instead" =
              is_binary(testData[[resp_var]]))

  # --- Pre-compute model class flags (once, outside loop) ---
  mc         <- class(testModel)
  is_glmmTMB <- "glmmTMB"  %in% mc
  is_gam     <- "gam"      %in% mc
  is_glmer   <- "glmerMod" %in% mc
  is_glm     <- "glm"      %in% mc && !is_gam
  # Note: lmer and negbin not included as they are not applicable to binary outcomes

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
        glmer(formula(testModel), family = family(testModel), data = train)
      } else if (is_glm) {
        glm(formula(testModel), family = family(testModel), data = train)
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
    } else if (is_glmer) {
      predict(m, type = "response", re.form = ~0,
              allow.new.levels = TRUE, newdata = newdata)
    } else {
      predict(m, type = "response", newdata = newdata)
    }
  }

  # --- val.prob helper: call once per dataset, extract both metrics ---
  # val.prob() computes AUC and Brier simultaneously; calling it twice
  # (once per metric) would be redundant and twice as slow.
  get_stats <- function(phat, y) {
    vp <- val.prob(p = phat, y = y, smooth = FALSE, pl = FALSE)
    c(auc = unname(vp["C (ROC)"]), brier = unname(vp["Brier"]))
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

    yhat_train <- get_preds(m_train, train)
    yhat_test  <- get_preds(m_train, test)

    stats_train <- get_stats(yhat_train, train[[resp_var]])
    stats_test  <- get_stats(yhat_test,  test[[resp_var]])

    results[[j]] <- data.frame(
      auc_train   = stats_train["auc"],
      brier_train = stats_train["brier"],
      auc_test    = stats_test["auc"],
      brier_test  = stats_test["brier"]
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
  # Note: separate() splits on "_" giving Metric (auc/brier) then Group (train/test),
  # which is the reverse order from RRMSE_RMAD_RBIAS where Group comes first.
  results_df <- bind_rows(results_clean, .id = "simRep") %>%
    pivot_longer(cols = -simRep, names_to = "metric", values_to = "value") %>%
    separate(metric, into = c("Metric", "Group")) %>%
    mutate(
      Group  = factor(Group, levels = c("train", "test"),
                      labels = c("In-sample performance", "Out-of-sample performance")),
      Metric = factor(Metric, levels = c("auc", "brier"),
                      labels = c("AUC statistic", "Brier score"))
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
    facet_grid(Group ~ Metric) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey90", linetype = "solid"),
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
          axis.text          = element_text(colour = "black"),
          panel.spacing      = unit(1.5, "lines")) +
    labs(x = "Value", y = "Frequency") +
    scale_x_continuous(breaks = seq(0, 1, 0.2),
                       expand = expansion(add = c(0.05, 0.05))) +
    coord_cartesian(xlim = c(0, 1))

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
    return(list(brier_auc_results  = results_df,
                brier_auc_hist     = results_plot,
                brier_auc_summary  = results_summary,
                dharmaPlot         = dharmaPlot))
  }

  list(brier_auc_results = results_df,
       brier_auc_hist    = results_plot,
       brier_auc_summary = results_summary)
}
