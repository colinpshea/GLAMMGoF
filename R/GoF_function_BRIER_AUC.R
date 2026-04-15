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
#' This package contains an example data set to fit a logistic regression called logitData, and four example logistic regression model objects are also included: logitModel1GLM is a GLM that includes a random effect; logitModel2GLM is a GLM that does not include a random effect; logitModel1GAM is a GAM that includes a random effect; and logitModel2GAM is a GAM that does not include a random effect. Both GLMs were fitted using `glmmTMB`, whereas the GAMs were fitted using `mgcv`. logitModel1GLM could also be a `glmer` (from `lme4`) or `gam` (from `mgcv`) model object, and logitModel2GLM could also be a `glm` or `gam` (from `mgcv`) model object:
#'
#' logitModel1GLM <- glmmTMB(y ~ totalLengthcm + Zone + (1|Year), family = binomial, data = logitData)
#'
#' logitModel2GLM <- glmmTMB(y ~ totalLengthcm + Zone, family = binomial, data = logitData)
#'
#' logitModel1GAM <- gam(y ~ s(totalLengthcm) + Zone + s(Year, bs = "re"), family = binomial, data = logitData)
#'
#' logitModel2GAM <- gam(y ~ s(totalLengthcm) + Zone, family = binomial, data = logitData)
#'
#' Bootstrapping the performance statistics requires specifying the data and model being tested, the desired number of bootstrap replicates (the default is 100 but it should be higher in practice), the proportion of data used in the training (in-sample performance) data set, whether you want to use DHARMa to assess the residuals (the default is TRUE), how many simulation replicates you want to use in DHARMa's `simulateResiduals()` function (the default is 1000), and an optional integer seed for reproducibility:
#'
#' BRIER_AUC(nReps = 100, testModel = logitModel1GLM, testData = logitData, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000, seed = 42)
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate bind_rows
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_histogram facet_grid theme_bw theme element_blank element_line element_text labs unit scale_x_continuous expansion
#' @importFrom DHARMa simulateResiduals
#' @importFrom rms val.prob
#' @importFrom glmmTMB ranef glmmTMB
#' @importFrom lme4 glmer
#' @importFrom mgcv gam predict.gam
#' @export
BRIER_AUC <- function(
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
  # Checks
  # ----------------------------
  stopifnot(!is.null(testModel))
  stopifnot(!is.null(testData))
  stopifnot(propTrain > 0 && propTrain < 1)

  resp_var <- all.vars(formula(testModel))[1]

  is_binary <- function(x)
    length(unique(x)) == 2 && all(x %in% c(0, 1))

  stopifnot(
    "Response variable is not binary! Use RRMSE_RMAD_RBIAS() instead" =
      is_binary(testData[[resp_var]])
  )

  mc <- class(testModel)

  # ----------------------------
  # CLASS-BASED MODEL DETECTION
  # ----------------------------
  is_glmmTMB <- "glmmTMB" %in% mc
  is_gam     <- "gam" %in% mc
  is_glmer   <- "glmerMod" %in% mc
  is_lmer    <- "lmerMod" %in% mc
  is_negbin  <- "negbin" %in% mc
  is_glm     <- "glm" %in% mc && !is_gam
  is_lm      <- "lm" %in% mc && !is_glm && !is_glm

  is_glmm_like <- is_glmmTMB || is_glmer || is_lmer

  # ----------------------------
  # GAM RE metadata
  # ----------------------------
  gam_re_labels <- NULL

  if (is_gam) {
    re_smooths <- Filter(function(s) isTRUE(s$random), testModel$smooth)
    if (length(re_smooths) > 0) {
      gam_re_labels <- sapply(re_smooths, function(s) s$label)
    }
  }

  # ----------------------------
  # FIT MODEL
  # ----------------------------
  fit_model <- function(train) {

    tryCatch({

      # GLMM
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

        # GAM
      } else if (is_gam) {

        mgcv::gam(formula(testModel),
                  family = family(testModel),
                  data   = train)

        # GLM
      } else if (is_glm) {

        glm(formula(testModel),
            family = family(testModel),
            data   = train)

        # LM
      } else if (is_lm) {

        lm(formula(testModel), data = train)

        # MASS negbin
      } else if (is_negbin) {

        MASS::glm.nb(formula(testModel), data = train)
      }

    }, error = function(e) {
      message("Model failed: ", e$message)
      NULL
    })
  }

  # ----------------------------
  # PREDICTION
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
  # METRICS
  # ----------------------------
  get_stats <- function(phat, y) {
    vp <- val.prob(p = phat, y = y, smooth = FALSE, pl = FALSE)
    c(auc = unname(vp["C (ROC)"]), brier = unname(vp["Brier"]))
  }

  # ----------------------------
  # LOOP
  # ----------------------------
  results <- vector("list", nReps)

  for (j in seq_len(nReps)) {

    train_idx <- sample(seq_len(nrow(testData)),
                        size = floor(propTrain * nrow(testData)))

    train <- testData[train_idx, ]
    test  <- testData[-train_idx, ]

    m_train <- fit_model(train)
    if (is.null(m_train)) next

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

  results_clean <- results[!vapply(results, is.null, logical(1))]

  if (length(results_clean) == 0) {
    stop("All model fits failed — no results.")
  }

  results_df <- bind_rows(results_clean, .id = "simRep") %>%
    tidyr::pivot_longer(cols = -simRep,
                        names_to = "metric",
                        values_to = "value") %>%
    tidyr::separate(metric, into = c("Metric", "Group")) %>%
    dplyr::mutate(
      Group = factor(Group,
                     levels = c("train", "test"),
                     labels = c("In-sample performance",
                                "Out-of-sample performance")),
      Metric = factor(Metric,
                      levels = c("auc", "brier"),
                      labels = c("AUC statistic", "Brier score"))
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

  # ----------------------------
  # DHARMa (guarded)
  # ----------------------------
  if (DHARMaPlot) {
    dharmaPlot <- tryCatch({
      DHARMa::simulateResiduals(n = DHARMaReps,
                                testModel, plot = TRUE)
    }, error = function(e) {
      warning("DHARMa plot failed: ", e$message)
      NULL
    })

    return(list(
      brier_auc_results = results_df,
      brier_auc_hist    = results_plot,
      brier_auc_summary = results_summary,
      dharmaPlot        = dharmaPlot
    ))
  }

  list(
    brier_auc_results = results_df,
    brier_auc_hist    = results_plot,
    brier_auc_summary = results_summary
  )
}
