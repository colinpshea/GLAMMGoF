#' Bootstrap Brier score and AUC fit statistics
#'
#' @description Bootstrap Brier score and AUC fit statistics (see `rms` package documentation for details) for generalized linear and generalized additive models with binary response variables and with or without random effects. Brier scores range from 0 to 1, with values closer to 0 indicating a better-predicting model, and where sqrt(Brier score) is the average difference, across all observations, between the predicted probability and the observed value (0 or 1). The AUC score is an aggregated metric that evaluates how well a logistic regression model classifies positive and negative outcomes at all possible cutoffs (i.e., probability cutoffs for defining positive (1) vs negative (0) cases). AUC statistics range from 0 to 1, where values closer to 1 indicate a better-predicting model (i.e., a better classifier), and where an AUC score of 0.5 suggests a model performs no better than random guessing. Note that although this function works with GLMMs and GAMMs fitted using a variety of packages/functions, all performance measures are based on model predictions that ignore random effects when present in a model (i.e., performance statistics are based on marginal or population-level model predictions).
#' @param nReps Desired number of bootstrap replicates. The default value is 100, but this number should be at least 1000 in practice.
#' @param testModel A logistic regression model fitted to testData using `glmmTMB` (with or without random effects), `glmer` (with random effects), `glm` (without random effects), or `gam` (with or without random effects).
#' @param testData A data frame with a binary response variable and continuous and/or categorical predictor variables.
#' @param propTrain The proportion of `testData` that is used for model-fitting and in-sample predictive performance (the default value is 0.8). The remaining % is used to assess out-of-sample predictive performance.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the `simulateResiduals()` function of the `DHARMa` package? The default is `TRUE`.
#' @param DHARMaReps If DHARMaPlot is `TRUE`, you can also specify DHARMaReps if you want something other than the default of 1000 simulation replicates.
#' @param seed Optional integer seed for reproducibility. If `NULL` (the default), no seed is set and results will differ across runs.
#' @return This function returns four objects: a data frame with all of the bootstrapping results (i.e., all `nReps` bootstrapped values for each performance statistic), a data frame with a summary (mean and 95% CLs) of all bootstrap replicates for each performance statistic, a histogram of values for each performance statistic, and a goodness-of-fit plot based on scaled residuals from the `simulateResiduals()` function of the `DHARMa` package. If DHARMaPlot = `FALSE`, then `simulateResiduals()` isn't used to assess the model's residuals and only three of the four objects are returned.
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
#' Bootstrapping the performance statistics requires specifying the data and model being tested, the desired number of bootstrap replicates (the default is 100 but it should be higher in practice), the proportion of data used in the training (in-sample performance) data set, whether you want to use DHARMa to assess the residuals (the default is TRUE), how many simulation replicates you want to use in DHARMa's `simulateResiduals()` function (the default is 1000), and an optional integer seed for reproducibility:
#'
#' BRIER_AUC(nReps = 100, testModel = logitModel_GLMM, testData = logitData, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000, seed = 42)
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
                      seed = NULL) {

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
  # Note: lmer/negbin not included — not applicable to binary outcomes

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
    train_idx <- sample(seq_len(nrow(testData)), size = floor(propTrain * nrow(testData)))
    train <- testData[ train_idx, ]
    test  <- testData[-train_idx, ]

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
    stop("All model fits failed — no results to summarise.")

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
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                       expand = expansion(add = c(0.05, 0.05)))

  if (DHARMaPlot) {
    dharmaPlot <- tryCatch(
      simulateResiduals(n = DHARMaReps, testModel, plot = TRUE),
      error = function(e) {
        warning("DHARMa plot failed: ", e$message)
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
