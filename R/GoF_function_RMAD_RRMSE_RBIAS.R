#' Bootstrap RRMSE, RMAD, and RBIAS predictive performance statistics
#'
#' @description Bootstrap RRMSE, RMAD, and RBIAS predictive performance statistics for generalized linear and generalized additive models with continuous or integer response variables and with or without random effects and zero-inflation. The performance statistics are relative root mean squared error (RRMSE), calculated as sqrt(mean((observed - predicted)^2))/mean(observed)*100, relative median absolute deviation (RMAD), calculated as median(abs((observed - predicted)))/mean(observed)*100, and relative bias (RBIAS), calculated as mean((observed - predicted))/mean(observed)*100. Note that although this function works with GLMMs and GAMMs fitted using a variety of packages/functions, all performance measures are based on model predictions that ignore random effects when present in a model (i.e., performance statistics are based on marginal or population-level model predictions). For models with a zero-inflation component, the performance statistics are based on predictions that account for zero-inflation (e.g., for glmmTMB and model objects, the predicted value represents the product of the mean_count and (1 - prob_zero)).
#' @param nReps Desired number of bootstrap replicates. The default value is 100, but this number should be at least 1000 in practice.
#' @param testModel A regression model fit to testData in `glmmTMB` (with or without random effects), `glmer` (with random effects), `glm`/`lm` (without random effects), or `gam` (with or without random effects). The response variable can be continuous or an integer, and possible statistical distributions include Poisson, negative binomial, gamma, tweedie, and gaussian.
#' @param testData A data frame with a continuous or integer response variable and continuous and/or categorical predictors.
#' @param propTrain Proportion of `testData` that is used for model-fitting and in-sample predictive performance (the remaining % is used to assess out-of-sample predictive performance). The default value is 0.8.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the `simulateResiduals()` function of the `DHARMa` package? The default is `TRUE`.
#' @param DHARMaReps You can also specify DHARMaReps if you want something other than the default of 1000 simulation replicates.
#' @param seed Optional integer seed for reproducibility. If `NULL` (the default), no seed is set and results will differ across runs.
#' @note This function does not currently support binomial models with cbind() or proportion responses, and for binary 0/1 responses, use BRIER_AUC(). This function also supports models with spatial random effects (e.g, in glmmTMB), but it is much slower than for more conventional GLM(M)s and GAM(M)s.
#' @return This function returns four objects: a data frame with all of the bootstrapping results (i.e., all nReps bootstrapped values for each performance statistic), a data frame with a summary (mean and 95% CLs) of all bootstrap replicates for each performance statistic, a histogram of values for each performance statistic, and a goodness-of-fit plot based on scaled residuals from the `simulateResiduals()` function of the `DHARMa` package. If DHARMaPlot = `FALSE`, then `simulateResiduals()` isn't used to assess the model's residuals, and only three of the four objects are returned.
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
#' Bootstrapping the performance statistics requires specifying the data and model being tested, the desired number of bootstrap replicates (the default is 100 but it should be higher in practice), the proportion of data used in the training (in-sample performance) data set, whether you want to use DHARMa to assess the residuals (the default is `TRUE`), how many simulation replicates you want to use in DHARMa's simulateResiduals() function (the default is 1000), and an optional integer seed for reproducibility:
#'
#' RRMSE_RMAD_RBIAS(nReps = 100, testModel = countModel_GLMM, testData = countData, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000, seed = 42)
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate bind_rows
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_histogram facet_grid theme_bw theme element_blank element_line element_text labs unit
#' @importFrom DHARMa simulateResiduals
#' @importFrom glmmTMB ranef glmmTMB
#' @importFrom lme4 glmer lmer glmer.nb
#' @importFrom MASS glm.nb
#' @importFrom mgcv gam predict.gam
#' @export
RRMSE_RMAD_RBIAS <- function(nReps = 100, testModel = NULL, testData = NULL,
                             propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000,
                             seed = NULL) {

  # --- Optional seed ---
  if (!is.null(seed)) set.seed(seed)

  # --- Cost functions ---
  fit_cost_rrmse <- function(y, yhat) sqrt(mean((y - yhat)^2)) / mean(y) * 100
  fit_cost_rmad  <- function(y, yhat) median(abs(y - yhat)) / mean(y) * 100
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
    "cbind() and proportion binomial responses are not supported. See ?RRMSE_RMAD_RBIAS for details." =
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
    train_idx <- sample(seq_len(nrow(testData)), size = floor(propTrain * nrow(testData)))
    train <- testData[ train_idx, ]
    test  <- testData[-train_idx, ]

    m_train <- fit_model(train)
    if (is.null(m_train)) next  # skip failed fits cleanly; no stale model carried forward

    y_train    <- train[[resp_var]]
    y_test     <- test[[resp_var]]
    yhat_train <- get_preds(m_train, train)
    yhat_test  <- get_preds(m_train, test)

    results[[j]] <- data.frame(
      train_RRMSE = fit_cost_rrmse(y_train, yhat_train),
      test_RRMSE  = fit_cost_rrmse(y_test,  yhat_test),
      train_RMAD  = fit_cost_rmad(y_train,  yhat_train),
      test_RMAD   = fit_cost_rmad(y_test,   yhat_test),
      train_RBIAS = fit_cost_rbias(y_train, yhat_train),
      test_RBIAS  = fit_cost_rbias(y_test,  yhat_test)
    )
  }

  # Guard against all replicates failing
  results_clean <- results[!vapply(results, is.null, logical(1))]
  if (length(results_clean) == 0)
    stop("All model fits failed — no results to summarise.")

  # --- Tidy results ---
  results_df <- bind_rows(results_clean, .id = "simRep") %>%
    pivot_longer(cols = -simRep, names_to = "metric", values_to = "value") %>%
    separate(metric, into = c("Group", "Metric")) %>%
    mutate(
      Group  = factor(Group, levels = c("train", "test"),
                      labels = c("In-sample performance", "Out-of-sample performance")),
      Metric = factor(Metric, levels = c("RRMSE", "RMAD", "RBIAS"))
    )

  results_summary <- results_df %>%
    group_by(Group, Metric) %>%
    summarise(mn    = mean(value),
              lwr95 = quantile(value, 0.025),
              upr95 = quantile(value, 0.975),
              .groups = "drop")

  results_plot <- ggplot(results_df, aes(x = value)) +
    geom_histogram(color = "black", fill = "grey") +
    facet_grid(Group ~ Metric, scales = "free") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey90", linetype = "solid"),
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"),
          axis.text          = element_text(colour = "black"),
          panel.spacing      = unit(1.5, "lines")) +
    labs(x = "% relative to true mean", y = "Frequency")

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
    return(list(rrmse_rmad_results  = results_df,
                rrmse_rmad_hist     = results_plot,
                rrmse_rmad_summary  = results_summary,
                dharmaPlot          = dharmaPlot))
  }

  list(rrmse_rmad_results = results_df,
       rrmse_rmad_hist    = results_plot,
       rrmse_rmad_summary = results_summary)
}
