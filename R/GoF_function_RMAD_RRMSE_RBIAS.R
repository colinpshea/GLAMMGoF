#' Bootstrap RRMSE, RMAD, and RBIAS predictive performance statistics
#'
#' @description Bootstrap RRMSE, RMAD, and RBIAS predictive performance statistics for generalized linear and generalized additive models with continuous or integer response variables and with or without random effects and zero-inflation. The performance statistics are relative root mean squared error (RRMSE), calculated as sqrt(mean((observed - predicted)^2))/mean(observed)*100, relative median absolute deviation (RMAD), calculated as median(abs((observed - predicted)))/mean(observed)*100, and relative bias (RBIAS), calculated as mean((observed - predicted))/mean(observed)*100. Note that although this function works with GLMMs and GAMMs fitted using a variety of packages/functions, all performance measures are based on model predictions that ignore random effects when present in a model (i.e., performance statistics are based on marginal or population-level model predictions). For models with a zero-inflation component, the performance statistics are based on predictions that account for zero-inflation (e.g., for glmmTMB and pscl model objects, the predicted value represents the product of the mean_count and (1 - prob_zero)).
#' @param nReps Desired number of bootstrap replicates. The default value is 100, but this number should be at least 1000 in practice.
#' @param testModel A regression model fit to testData in `glmmTMB` (with or without random effects), `glmer` (with random effects), `glm`/`lm` (without random effects), or `gam` (with or without random effects). The response variable can be continuous or an integer, and possible statistical distributions include Poisson, negative binomial, gamma, tweedie, and gaussian.
#' @param testData A data frame with a continuous or integer response variable and continuous and/or categorical predictors.
#' @param propTrain Proportion of `testData` that is used for model-fitting and in-sample predictive performance (the remaining % is used to assess out-of-sample predictive performance). The default value is 0.8.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the `simulateResiduals()` function of the `DHARMa` package? The default is `TRUE`.
#' @param DHARMaReps You can also specify DHARMaReps if you want something other than the default of 1000 simulation replicates.
#' @param seed Optional integer seed for reproducibility. If `NULL` (the default), no seed is set and results will differ across runs.
#' @return This function returns four objects: a data frame with all of the bootstrapping results (i.e., all nReps bootstrapped values for each performance statistic), a data frame with a summary (mean and 95% CLs) of all bootstrap replicates for each performance statistic, a histogram of values for each performance statistic, and a goodness-of-fit plot based on scaled residuals from the `simulateResiduals()` function of the `DHARMa` package. If DHARMaPlot = `FALSE`, then `simulateResiduals()` isn't used to assess the model's residuals, and only three of the four objects are returned.
#'
#' This package contains an example data set for fitting a negative binomial or Poisson regression called countData. This data set has an integer response variable, but data with a continuous response variable could also be used. Four example negative binomial regression model objects are also included: countModel1GLM is a GLM that includes a random effect; countModel2GLM is a GLM that does not include a random effect; countModel1GAM is a GAM that includes a random effect; and countModel2GAM is a GAM that does not include a random effect. Both GLMs were fitted using `glmmTMB`, whereas the GAMs were fitted using `mgcv`. countModel1 could also be a `glmer` (fitted using `glmer` or `glmer.nb` from `lme4`) or `gam` (fitted using `mgcv`) model object, and countModel2 could also be a `glm.nb` (from the `MASS` package) or `gam` model object:
#'
#' countModel1GLM <- glmmTMB(y ~ Season + River + Temp + Snags + Year + AvgDepth + (1|RiverSeasonYear), family = nbinom2, data = countData)
#'
#' countModel2GLM <- glmmTMB(y ~ Season + River + Temp + Snags + Year + AvgDepth, family = nbinom2, data = countData)
#'
#' countModel1GAM <- gam(y ~ Season + River + s(Temp) + Snags + Year + s(AvgDepth) + s(RiverSeasonYear, bs = "re"), family = nb, data = countData)
#'
#' countModel2GAM <- gam(y ~ Season + River + s(Temp) + Snags + Year + s(AvgDepth), family = nb, data = countData)
#'
#' Bootstrapping the performance statistics requires specifying the data and model being tested, the desired number of bootstrap replicates (the default is 100 but it should be higher in practice), the proportion of data used in the training (in-sample performance) data set, whether you want use DHARMa to assess the residuals (the default is `TRUE`), how many simulation replicates you want to use in DHARMa's simulateResiduals() function (the default is 1000), and an optional integer seed for reproducibility:
#'
#' RRMSE_RMAD_RBIAS(nReps = 100, testModel = countModel1, testData = countData, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000, seed = 42)
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
RRMSE_RMAD_RBIAS <- function(
    nReps = 100,
    testModel = NULL,
    testData = NULL,
    propTrain = 0.8,
    DHARMaPlot = TRUE,
    DHARMaReps = 1000,
    seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  # ----------------------------
  # Metrics
  # ----------------------------
  fit_cost_rrmse <- function(y, yhat)
    sqrt(mean((y - yhat)^2)) / mean(y) * 100

  fit_cost_rmad <- function(y, yhat)
    median(abs(y - yhat)) / mean(y) * 100

  fit_cost_rbias <- function(y, yhat)
    mean(y - yhat) / mean(y) * 100

  # ----------------------------
  # Checks
  # ----------------------------
  stopifnot(!is.null(testModel))
  stopifnot(!is.null(testData))
  stopifnot(propTrain > 0 && propTrain < 1)

  resp_var <- all.vars(formula(testModel))[1]

  is_binary <- function(x)
    length(unique(x)) == 2 && all(x %in% c(0, 1))

  stopifnot(
    "Response variable is binary! Use BRIER_AUC() instead" =
      !is_binary(testData[[resp_var]])
  )

  mc <- class(testModel)

  is_glmmTMB <- "glmmTMB" %in% mc
  is_gam     <- "gam" %in% mc
  is_glmer   <- "glmerMod" %in% mc
  is_lmer    <- "lmerMod" %in% mc
  is_negbin  <- "negbin" %in% mc
  is_glm     <- "glm" %in% mc && !is_gam
  is_lm      <- "lm" %in% mc && !is_gam && !is_glm

  # ----------------------------
  # BULLETPROOF GLMM DETECTION
  # ----------------------------
  has_re_syntax <- grepl("\\|", deparse(formula(testModel)))

  is_glmm_like <- has_re_syntax && (
    is_glmmTMB || is_glmer || is_lmer
  )

  # ----------------------------
  # GAM RE metadata
  # ----------------------------
  gam_re_labels <- NULL
  gam_re_terms  <- NULL

  if (is_gam) {
    re_smooths <- Filter(function(s) isTRUE(s$random), testModel$smooth)
    if (length(re_smooths) > 0) {
      gam_re_labels <- sapply(re_smooths, function(s) s$label)
      gam_re_terms  <- sapply(re_smooths, function(s) s$term)
    }
  }

  # ----------------------------
  # Fit model
  # ----------------------------
  fit_model <- function(train) {

    tryCatch({

      # =========================
      # GLMM (SAFE ROUTING)
      # =========================
      if (is_glmm_like) {

        if (is_glmmTMB) {

          glmmTMB(
            formula(testModel, component = "cond"),
            family      = family(testModel),
            dispformula = formula(testModel, component = "disp"),
            ziformula   = formula(testModel, component = "zi"),
            data        = train
          )

        } else {

          fam <- family(testModel)

          if (grepl("Negative Binomial", fam$family)) {
            lme4::glmer.nb(formula(testModel), data = train)
          } else {
            lme4::glmer(formula(testModel),
                        family = fam,
                        data   = train)
          }
        }

        # =========================
        # GAM
        # =========================
      } else if (is_gam) {

        mgcv::gam(formula(testModel),
                  family = family(testModel),
                  data   = train)

        # =========================
        # GLM / LM
        # =========================
      } else if (is_glm) {

        glm(formula(testModel),
            family = family(testModel),
            data   = train)

      } else if (is_lm) {

        lm(formula(testModel), data = train)

        # =========================
        # fallback GLMM objects
        # =========================
      } else if (is_glmer) {

        if (grepl("Negative Binomial", family(testModel)$family))
          lme4::glmer.nb(formula(testModel), data = train)
        else
          lme4::glmer(formula(testModel),
                      family = family(testModel),
                      data   = train)

      } else if (is_lmer) {

        lme4::lmer(formula(testModel), data = train)

      } else if (is_negbin) {

        MASS::glm.nb(formula(testModel), data = train)
      }

    }, error = function(e) {
      message("Model failed: ", e$message)
      NULL
    })
  }

  # ----------------------------
  # Prediction helper
  # ----------------------------
  get_preds <- function(m, newdata) {

    if (is_glmmTMB) {

      predict(m, type = "response", newdata = newdata)

    } else if (is_gam) {

      if (!is.null(gam_re_labels)) {
        predict(m,
                type = "response",
                exclude = gam_re_labels,
                newdata = newdata,
                newdata.guaranteed = TRUE)
      } else {
        predict(m, type = "response", newdata = newdata)
      }

    } else if (is_glmer || is_lmer) {

      predict(m, type = "response", re.form = ~0, newdata = newdata)

    } else {

      predict(m, type = "response", newdata = newdata)
    }
  }

  # ----------------------------
  # Loop
  # ----------------------------
  results <- vector("list", nReps)

  for (j in seq_len(nReps)) {

    train_idx <- sample(seq_len(nrow(testData)),
                        size = floor(propTrain * nrow(testData)))

    train <- testData[train_idx, ]
    test  <- testData[-train_idx, ]

    m_train <- fit_model(train)
    if (is.null(m_train)) next

    y_train <- train[[resp_var]]
    y_test  <- test[[resp_var]]

    yhat_train <- get_preds(m_train, train)
    yhat_test  <- get_preds(m_train, test)

    results[[j]] <- data.frame(
      train_RRMSE = fit_cost_rrmse(y_train, yhat_train),
      test_RRMSE  = fit_cost_rrmse(y_test, yhat_test),
      train_RMAD  = fit_cost_rmad(y_train, yhat_train),
      test_RMAD   = fit_cost_rmad(y_test, yhat_test),
      train_RBIAS = fit_cost_rbias(y_train, yhat_train),
      test_RBIAS  = fit_cost_rbias(y_test, yhat_test)
    )
  }

  # ----------------------------
  # SAFE bind_rows (fixes simRep crash)
  # ----------------------------
  results_clean <- results[!vapply(results, is.null, logical(1))]

  if (length(results_clean) == 0) {
    stop("All model fits failed — no results to summarise.")
  }

  results_df <- bind_rows(results_clean, .id = "simRep") %>%
    tidyr::pivot_longer(cols = -simRep,
                        names_to = "metric",
                        values_to = "value") %>%
    tidyr::separate(metric,
                    into = c("Group", "Metric")) %>%
    dplyr::mutate(
      Group = factor(Group,
                     levels = c("train", "test"),
                     labels = c("In-sample performance",
                                "Out-of-sample performance")),
      Metric = factor(Metric,
                      levels = c("RRMSE", "RMAD", "RBIAS"))
    )

  results_summary <- results_df %>%
    dplyr::group_by(Group, Metric) %>%
    dplyr::summarise(
      mn    = mean(value),
      lwr95 = quantile(value, 0.025),
      upr95 = quantile(value, 0.975),
      .groups = "drop"
    )

  results_plot <- ggplot2::ggplot(results_df,
                                  ggplot2::aes(x = value)) +
    ggplot2::geom_histogram(color = "black", fill = "grey") +
    ggplot2::facet_grid(Group ~ Metric) +
    ggplot2::theme_bw()

  if (DHARMaPlot) {
    dharmaPlot <- DHARMa::simulateResiduals(
      n = DHARMaReps,
      testModel,
      plot = TRUE
    )

    return(list(
      rrmse_rmad_results = results_df,
      rrmse_rmad_hist    = results_plot,
      rrmse_rmad_summary = results_summary,
      dharmaPlot         = dharmaPlot
    ))
  }

  list(
    rrmse_rmad_results = results_df,
    rrmse_rmad_hist    = results_plot,
    rrmse_rmad_summary = results_summary
  )
}
