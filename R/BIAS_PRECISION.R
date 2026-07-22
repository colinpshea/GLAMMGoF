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
#' @param bias_adjust Character string specifying the bias adjustment method for marginal predictions in `glmmTMB` and `lme4` models. One of `"none"` (the default), `"manual"`, or `"tmb"`. `"none"` uses standard population-level predictions (`re.form = ~0`) with no correction, preserving the full RBIAS signal driven by Jensen's inequality and allowing it to be used diagnostically. `"manual"` applies an analytical correction to marginal predictions: population-level predictions (`re.form = ~0`) are multiplied by `exp(V / 2)`, where `V` is the sum of the variance components that are additive on the log scale, i.e. that sit inside the exponential back-transformation. For log-link models `V` is the total random effect variance summed across all RE terms extracted from `VarCorr()`. For lognormal models - those fit to a natural-log-transformed response, `log(y) ~ .` - the residual variance is also additive on the log scale, so `V` additionally includes `sigma(testModel)^2`; such models are subject to retransformation bias from both sources simultaneously and a correction omitting the residual term will under-correct, often substantially, since the residual variance is frequently the larger component. Whether the residual term is included is determined by inspecting the response transformation in the model formula, not the family or link, so a Gaussian model with a log *link* on an untransformed response (`family = gaussian("log")`, where the residual is additive on the response scale and outside the exponential) correctly receives the random-effect-only correction. This option is supported for `glmmTMB` and `lme4` model objects, and for `lm`/`glm` objects with a log-transformed response. `"tmb"` uses TMB's built-in bias correction (`do.bias.correct = TRUE`), which requires `re.form = NULL` as a computational necessity so that TMB's automatic differentiation can access the RE variance structure; this is not a choice to make group-specific conditional predictions, and the output is still a corrected population-level marginal mean equivalent to `exp(X*beta) * exp(sigma^2/2)` for standard random intercept models; note that `"tmb"` is only available for `glmmTMB` models and will throw an informative error for `lme4` objects, and is additionally unavailable for lognormal models, where TMB's correction integrates over the random effects but omits the residual retransformation term and would therefore under-correct. For most diagnostic purposes `"none"` is recommended; for calibrated response-scale predictions `"manual"` and `"tmb"` should give similar results and can be compared as a consistency check. This argument is silently ignored for `mgcv` GAM/GAMM models. See the note below and Thorson & Kristensen (2016) for details.
#' @note This function does not currently support binomial models with cbind() or proportion responses, and for binary 0/1 responses, use brier_auc(). This function also supports models with spatial random effects (e.g, in glmmTMB), but it is much slower than for more conventional GLM(M)s and GAM(M)s.
#'
#' **Random effects and Jensen's inequality:** By default, all predictions are population-level (i.e., random effects are set to zero via `re.form = ~0`). For models with a nonlinear link function (e.g., log, logit) and random effects, back-transforming the linear predictor to the response scale introduces a systematic negative bias in the predicted arithmetic mean. This occurs because Jensen's inequality implies that `E[exp(eta)] > exp(E[eta])` for any random variable `eta`, so `exp(beta_0)` underestimates the true mean by a factor of approximately `exp(sigma^2 / 2)`, where `sigma^2` is the total random effect variance (summed across all RE terms). It follows that log-normal models with random effects are subject to both sources of bias simultaneously (residual variance and random effect variance), requiring a combined correction of `exp((sigma^2_residual + sigma^2_random_effect) / 2)`, where both variance components are additive on the natural log scale. The magnitude of this bias grows rapidly with RE variance: a model with two random effects of modest size (e.g., SD = 0.3 and 0.4) can produce marginal predictions that underestimate observed values by 10% or more. Consistent negative `RBIAS` in the output of this function - particularly when both in-sample and out-of-sample values are negative - may therefore reflect this structural property of the model rather than misspecification. To obtain bias-corrected marginal predictions, set `bias_adjust = "manual"` (available for `glmmTMB` and `lme4` models, and for `lm`/`glm` fits to a log-transformed response) or `bias_adjust = "tmb"` (`glmmTMB` only, and not applicable to lognormal models). For lognormal models the correction applied is the combined `exp((sigma^2_residual + sum(sigma^2_RE)) / 2)`. Note that no correction is available for `mgcv` models.
#'
#' **Random slopes:** `bias_adjust = "manual"` assumes a random intercept-only structure. For models with random slopes, the variance of the linear predictor varies by observation as a function of the covariate associated with the random slope, so a single scalar correction `exp(sigma^2 / 2)` is not appropriate. In this case, use `bias_adjust = "tmb"`, which applies TMB's automatic differentiation over the full random effect covariance structure and correctly accounts for random slope variance. Applying `bias_adjust = "manual"` to a random slopes model will produce incorrect corrections and is not recommended; a warning is now issued when random slopes are detected. In such cases the correction is computed from the random intercept variances only, and slope variances and covariances are ignored.
#'
#' **Conditional predictions:** By default, both in-sample and out-of-sample predictions are made marginally (`re.form = ~0`) to assess population-level generalization. Setting `conditional_predictions = TRUE` switches both to conditional predictions (`re.form = NULL`, `allow.new.levels = TRUE`), which use the random effect estimates from the refitted training model for all predictions. Since GLAMMGoF's holdout CV splits rows within groups rather than groups themselves, the training model has RE estimates for all groups represented in the test set, making conditional out-of-sample predictions well-defined and directly comparable to conditional in-sample predictions. The gap between in-sample and out-of-sample metrics then reflects genuine row-level overfitting rather than any marginal vs conditional distinction. This is most appropriate when the study design involves repeated measurements within observed groups (sites, years, individuals) and the primary interest is within-group predictive accuracy rather than generalization to new groups. Note that `conditional_predictions = TRUE` is incompatible with `bias_adjust = "manual"` and will throw an informative error - conditional predictions already absorb Jensen's inequality bias via estimated random effects, so applying the lognormal correction on top would double-correct. Use `bias_adjust = "tmb"` with `conditional_predictions = TRUE` if a bias correction is desired alongside conditional predictions, or use `bias_adjust = "manual"` with `conditional_predictions = FALSE` for population-level generalization with Jensen bias correction. For lognormal models the reasoning behind this restriction differs slightly: conditional predictions absorb the random effect contribution to retransformation bias, but not the residual contribution, since the residual is always marginalized. A conditional lognormal prediction therefore still requires a residual-only correction of `exp(sigma^2_residual / 2)`, which is not currently implemented; the incompatibility error applies to lognormal models as well.
#'
#' **Correction factor and out-of-sample leakage:** For `bias_adjust = "manual"`, the correction factor `exp(V / 2)` is computed once from the full-data `testModel` via `VarCorr()` for the random effect variances, plus `sigma()` for the residual variance in lognormal models,  and applied uniformly across all resampling replicates rather than being recomputed from each refitted training model. This introduces a minor form of information leakage into out-of-sample performance metrics, since the correction factor incorporates variance information from the full dataset including held-out observations. However, this approach is intentional: recomputing the correction from each training resample produces highly unstable correction factors at high RE variance, where small-sample RE variance estimates are noisy and `exp(sigma^2/2)` is sensitive to overestimation. The stability of the bias-corrected metrics is considered a greater practical benefit than strict out-of-sample purity of the correction factor, particularly given that the primary purpose of `bias_adjust = "manual"` is diagnostic confirmation that Jensen's inequality is the source of observed negative RBIAS rather than model misspecification. While `bias_adjust = "manual"` is an approximation rather than an exact correction - unlike `bias_adjust = "tmb"` which uses TMB's automatic differentiation over the full RE covariance structure - it is a close and reliable approximation for models with random intercepts, and is recommended as a fast and interpretable diagnostic tool: if switching from `bias_adjust = "none"` to `bias_adjust = "manual"` moves RBIAS toward zero, this is strong evidence that the observed negative bias is attributable to RE variance and Jensen's inequality rather than model misspecification or poor predictive generalization.
#' **Recommended workflow:** The following stepwise workflow is recommended for assessing predictive performance of a glmmTMB GLMM with a log link. See the package vignette for a fully worked example.
#'
#' Step 1 - Run with defaults (`bias_adjust = "none"`, `conditional_predictions = FALSE`): assess population-level predictive generalization. Inspect RBIAS for both in-sample and out-of-sample performance. A message will fire if both are consistently negative (< -10%), suggesting Jensen's inequality may be contributing.
#'
#' Step 2 - If negative RBIAS is detected, re-run with `bias_adjust = "manual"`: if RBIAS moves toward zero, Jensen's inequality is confirmed as the source of the negative bias rather than model misspecification or poor generalization. The analytical correction `exp(sigma^2/2)` should then be applied to any response-scale predictions used for management or reporting purposes.
#'
#' Step 3 (optional) - Re-run with `bias_adjust = "tmb"` to validate the manual correction using TMB's automatic differentiation. If manual and TMB agree closely, the analytical correction is reliable for your RE structure. Note that this step is considerably slower than Step 2, particularly for models with spatial random effects.
#'
#' Step 4 (optional) - Re-run with `conditional_predictions = TRUE` to assess within-group predictive accuracy using estimated random effects. This measures how well the model interpolates within observed groups rather than generalizing to new ones. The in-sample vs out-of-sample gap now reflects genuine row-level overfitting rather than any marginal vs conditional distinction.
#'
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
#' Bootstrapping or Monte Carlo resampling of the performance statistics requires specifying the data and model being tested, the desired number of replicates (the default is 100 but should be at least 1000 in practice),  the proportion of data used for training when `method = "holdout"` (the default is 0.8), whether to use `DHARMa` residual diagnostics (the default is `TRUE`), whether to use `DHARMa` to test for zero-inflation (the default is `TRUE`), the number of `DHARMa` simulation replicates (the default is 1000), and an optional integer seed for reproducibility, and the resampling method `holdout` or `bootstrap`. Standard usage: bias_precision(nReps = 100, testModel = countModel_GLMM, testData = countData, propTrain = 0.8, DHARMaPlot = TRUE, testZI = TRUE, DHARMaReps = 1000, seed = 123, method = "holdout").
#'
#' To diagnose Jensen's inequality bias, first run with bias_adjust = "none" and inspect RBIAS, then re-run with bias_adjust = "manual". If RBIAS moves toward zero, Jensen's inequality is the source of the negative bias rather than model misspecification: bias_precision(nReps = 100, testModel = countModel_GLMM, testData = countData, method = "holdout", DHARMaPlot = FALSE, bias_adjust = "manual", seed = 123).
#'
#' To validate the manual correction against TMB's AD-based method (slower): bias_precision(nReps = 100, testModel = countModel_GLMM, testData = countData, method = "holdout", DHARMaPlot = FALSE, bias_adjust = "tmb", seed = 123).
#'
#' To use conditional predictions for both in-sample and out-of-sample performance; this is appropriate when repeated measurements within observed groups are present and within-group predictive accuracy is the primary interest, and both in-sample and out-of-sample RBIAS will reflect how well the model predicts held-out observations from groups it has already seen, and the in/out gap reflects genuine row-level overfitting: bias_precision(nReps = 100, testModel = countModel_GLMM, testData = countData, method = "holdout", DHARMaPlot = FALSE, conditional_predictions = TRUE, seed = 123).
#'
#' To apply the bias correction to new predictions outside of GLAMMGoF after diagnosing Jensen's inequality, use jensen_correct(), which resolves the appropriate correction automatically from the model, or extract the variance components manually (the total RE variance from VarCorr() for log-link models, plus the residual variance from sigma() for models fit to a log-transformed response), compute the correction factor as exp(total_variance / 2), and multiply marginal predictions by this factor (for random intercept models); or use predict() with do.bias.correct = TRUE for models with random slopes or complex RE structures. See the package vignette for worked examples.
#'
#' @param verbose Logical. If `TRUE` (the default), prints a diagnostic message when substantial negative RBIAS is detected in a `glmmTMB` model or a lognormal model with `bias_adjust = "none"`, suggesting the user consider applying a bias correction. Set to `FALSE` to suppress this message, which is useful when calling `bias_precision()` repeatedly in simulation or sweep contexts.
#' @param conditional_predictions Logical. If `TRUE` and the model is a `glmmTMB` or `lme4` fit with random effects, both in-sample and out-of-sample predictions are made conditionally on the estimated random effects (`re.form = NULL`, `allow.new.levels = TRUE`) rather than marginally (`re.form = ~0`). Since GLAMMGoF's holdout CV splits rows within groups rather than groups themselves, random effect estimates from the training model are valid for held-out observations from the same groups, making conditional out-of-sample predictions well-defined and meaningful. This produces performance metrics that reflect within-group predictive accuracy - how well the model predicts held-out observations from groups it has already seen - which is the most relevant question when the study design involves repeated measurements within sites, years, or individuals. In-sample and out-of-sample metrics remain directly comparable since both use the same conditional prediction strategy; the gap between them reflects genuine overfitting to training rows rather than any marginal vs conditional distinction. The default is `FALSE`, which uses marginal predictions for both in-sample and out-of-sample performance when `bias_adjust` is `"none"` or `"manual"`, assessing population-level generalization rather than within-group interpolation. Note that `bias_adjust = "tmb"` always uses conditional predictions internally regardless of this setting, since TMB's bias correction requires random effects to be present in the linear predictor; `conditional_predictions` therefore has no additional effect when `bias_adjust = "tmb"`. This argument is supported for `glmmTMB` and `lme4` model objects and is silently ignored for `mgcv` GAM/GAMM models.
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by summarise mutate bind_rows
#' @importFrom tidyr pivot_longer separate
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline facet_grid theme_bw theme element_blank element_line element_text labs unit scale_y_continuous
#' @param correction_factor Optional numeric scalar. If supplied and `bias_adjust = "manual"`, this value is used directly as the lognormal bias correction factor instead of computing it from `VarCorr(testModel)`. This is primarily useful in simulation contexts where the true random effect variance is known, allowing the correction to be anchored to the true `sigma^2` rather than the potentially noisy estimated variance from the fitted model. For example, in a simulation sweeping over sigma values, pass `exp(sigma^2 / 2)` directly to avoid instability at high sigma where `VarCorr()` estimates from training subsets are unreliable. In applied use with real data, leave this as `NULL` (the default) so the correction is computed from the full-data `testModel` via `VarCorr()`. Ignored when `bias_adjust` is not `"manual"`. For lognormal models, note that a supplied `correction_factor` replaces the entire correction, so it must include the residual variance as well as the random effect variance, i.e. `exp((sigma^2_residual + sigma^2_RE) / 2)`; supplying `exp(sigma^2_RE / 2)` will silently under-correct.
#' @importFrom DHARMa simulateResiduals testZeroInflation
#' @importFrom glmmTMB ranef glmmTMB
#' @importFrom nlme VarCorr
#' @importFrom lme4 glmer lmer glmer.nb
#' @importFrom MASS glm.nb
#' @importFrom mgcv gam predict.gam
#' @importFrom stats complete.cases formula predict sigma
#' @export
bias_precision <- function(nReps = 100, testModel = NULL, testData = NULL,
                           propTrain = 0.8, DHARMaPlot = TRUE, testZI = TRUE, DHARMaReps = 1000,
                           seed = NULL, method = c("holdout", "bootstrap"),
                           bias_adjust = c("none", "manual", "tmb"),
                           verbose = TRUE,
                           conditional_predictions = FALSE,
                           correction_factor = NULL) {

  # --- specify bootstrapping method and bias adjustment
  method      <- match.arg(method)
  bias_adjust <- match.arg(bias_adjust)

  # --- Optional seed ---
  if (!is.null(seed)) set.seed(seed)

  # --- Validate argument types ---
  if (inherits(testModel, "data.frame"))
    stop("testModel appears to be a data frame rather than a fitted model object. ",
         "Did you accidentally swap the testModel and testData arguments? ",
         "testModel should be a fitted model object (e.g., from glmmTMB, lme4, or mgcv) ",
         "and testData should be the data frame used to fit the model.")

  if (!inherits(testData, "data.frame"))
    stop("testData should be a data frame but received an object of class '",
         paste(class(testData), collapse = ", "), "'. ",
         "Did you accidentally swap the testModel and testData arguments?")

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

  # --- Detect a natural-log-transformed response (lognormal models) ---
  # A lognormal model (log(y) ~ .) is Gaussian on the log scale, so predict()
  # returns the mean on the LOG scale: predictions must be exponentiated before
  # comparison with the original-scale response, and the correction adds the
  # residual variance to the RE variance. Detection keys on the formula LHS,
  # never on family/link, so gaussian(link = "log") (bare LHS, untransformed
  # response) is NOT treated as lognormal and correctly receives no residual term.
  logresp      <- .detect_log_response(testModel)
  is_lognormal <- isTRUE(logresp$natural)
  if (logresp$logged && !is_lognormal)
    stop("bias_precision() supports natural-log response transforms only; the ",
         "formula uses '", logresp$base, "'. The back-transformation and its ",
         "variance scaling differ for other bases - refit on the natural-log ",
         "scale (log(y) ~ .).")

  # --- conditional_predictions message (after model class flags defined) ---
  if (conditional_predictions && (is_glmmTMB || is_glmer || is_lmer)) {
    if (bias_adjust == "manual")
      stop("bias_adjust = 'manual' is incompatible with conditional_predictions = TRUE. ",
           "Conditional predictions already absorb Jensen's inequality bias via estimated ",
           "random effects - applying the lognormal correction on top would double-correct. ",
           "Use bias_adjust = 'none' with conditional_predictions = TRUE to assess ",
           "within-group predictive accuracy, or bias_adjust = 'manual' with ",
           "conditional_predictions = FALSE to assess population-level generalization ",
           "with Jensen bias correction. See ?bias_precision for details.")
    message(
      "Note: conditional_predictions = TRUE assesses within-group predictive ",
      "accuracy using estimated random effects (re.form = NULL) for both ",
      "in-sample and out-of-sample predictions. This is appropriate when ",
      "repeated measurements within observed groups are present and within-group ",
      "interpolation is the primary interest. For population-level generalization ",
      "assessment - including detection of Jensen's inequality bias - use the ",
      "default conditional_predictions = FALSE."
    )
  }

  # --- Validate bias_adjust = "tmb" for lme4 models ---
  if (bias_adjust == "tmb" && (is_glmer || is_lmer))
    stop("bias_adjust = 'tmb' is not supported for lme4 models since TMB's automatic ",
         "differentiation bias correction is only available in glmmTMB. ",
         "Use bias_adjust = 'manual' instead, or refit your model in glmmTMB to ",
         "access the TMB correction. See ?bias_precision for details.")

  # --- Validate bias_adjust = "tmb" for lognormal models ---
  if (bias_adjust == "tmb" && is_lognormal)
    stop("bias_adjust = 'tmb' does not apply to lognormal (log-transformed ",
         "response) models: TMB's bias correction integrates over the random ",
         "effects but does not include the residual retransformation term, so it ",
         "would under-correct. Use bias_adjust = 'manual', which applies ",
         "exp((sigma^2_resid + sum(sigma^2_RE)) / 2).")

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

  # --- Pre-compute manual bias correction factor from full-data model (once, outside loop) ---
  # If correction_factor is supplied by the user (e.g. in simulation contexts where the
  # true sigma^2 is known), use it directly. Otherwise compute it once from testModel
  # via the shared .re_variance_sum() helper (and sigma() for the lognormal residual).
  # Using testModel rather than each resample's refitted model avoids unstable
  # correction factors at high RE variance where training subsets produce noisy
  # variance estimates that inflate exp(V/2) dramatically.
  #
  # V = sum of RE intercept variances (all RE models) + residual variance (lognormal
  # only). The residual term is what makes a log(y) ~ . model doubly biased; it is
  # omitted for log-link models, where the residual sits outside exp(). This branch
  # also covers no-RE lognormal fits (lm / glm), for which re_var = 0 and the
  # correction reduces to exp(sigma^2_resid / 2).
  cf_internal <- if (!is.null(correction_factor) && bias_adjust == "manual") {
    correction_factor   # use supplied value directly
  } else if (bias_adjust == "manual") {
    rev_full <- .re_variance_sum(testModel)
    if (rev_full$slopes)
      warning("Random slopes detected; bias_adjust = 'manual' uses intercept ",
              "variances only and is approximate. For random-slope models use ",
              "bias_adjust = 'tmb' (glmmTMB). See ?bias_precision.")
    resid_var <- if (is_lognormal) sigma(testModel)^2 else 0
    exp((rev_full$re_var + resid_var) / 2)
  } else {
    1
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

  # --- Predict helpers ---
  # Both train and test use the same prediction strategy via get_preds_base.
  # When conditional_predictions = TRUE, both use re.form = NULL since
  # GLAMMGoF's holdout CV splits rows within groups rather than groups
  # themselves - RE estimates from the training model are valid for
  # held-out observations from the same groups, making conditional
  # out-of-sample predictions well-defined and directly comparable to
  # conditional in-sample predictions.
  # When conditional_predictions = FALSE (default), both use marginal
  # predictions for population-level generalization assessment.

  get_preds_base <- function(m, newdata) {
    # --- Lognormal (log(y) ~ .): predict() returns the mean on the LOG scale ---
    # The response is Gaussian on the log scale (identity link), so predict()
    # gives the log-scale mean regardless of backend. Back-transform with exp()
    # to the original response scale, then apply the correction factor
    # (residual + RE variance via cf_internal when bias_adjust = "manual"; 1
    # otherwise). This runs before the branches below so the log-scale-vs-
    # original-scale mismatch is resolved for every lognormal backend, including
    # no-RE lm / glm fits.
    if (is_lognormal) {
      reform <- if (conditional_predictions) NULL else ~0
      if (is_glmmTMB) {
        eta <- predict(m, type = "response", newdata = newdata,
                       re.form = reform, allow.new.levels = TRUE)
      } else if (is_glmer || is_lmer) {
        eta <- predict(m, newdata = newdata, re.form = reform,
                       allow.new.levels = TRUE)
      } else {
        eta <- predict(m, newdata = newdata)          # lm / glm (no random effects)
      }
      return(exp(eta) * cf_internal)
    }

    if (is_glmmTMB) {
      if (conditional_predictions) {
        predict(m, type = "response", newdata = newdata,
                re.form = NULL, allow.new.levels = TRUE)
      } else if (bias_adjust == "manual") {
        predict(m, type = "response", newdata = newdata,
                re.form = ~0, allow.new.levels = TRUE) * cf_internal
      } else if (bias_adjust == "tmb") {
        preds <- predict(m, type = "response", newdata = newdata,
                         re.form = NULL, allow.new.levels = TRUE,
                         do.bias.correct = TRUE)
        preds[, "Est. (bias.correct)"]
      } else {
        predict(m, type = "response", newdata = newdata,
                re.form = ~0, allow.new.levels = TRUE)
      }
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
      if (conditional_predictions) {
        predict(m, type = "response", re.form = NULL,
                allow.new.levels = TRUE, newdata = newdata)
      } else if (bias_adjust == "manual") {
        predict(m, type = "response", re.form = ~0,
                allow.new.levels = TRUE, newdata = newdata) * cf_internal
      } else {
        predict(m, type = "response", re.form = ~0,
                allow.new.levels = TRUE, newdata = newdata)
      }
    } else {
      predict(m, type = "response", newdata = newdata)
    }
  }

  get_preds_train <- get_preds_base
  get_preds_test  <- get_preds_base

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
    yhat_train <- get_preds_train(m_train, train)
    yhat_test  <- get_preds_test(m_train, test)

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
  # on a nonlinear link scale. Only fires when bias_adjust = "none" and verbose = TRUE.
  if (verbose && bias_adjust == "none" && (is_glmmTMB || is_lognormal)) {
    rbias_summary <- results_summary[results_summary$Metric == "RBIAS", ]
    both_negative <- all(rbias_summary$mn < -10)
    if (both_negative) {
      message(
        "Note: Both in-sample and out-of-sample RBIAS are consistently negative (< -10%). ",
        "This may indicate that strong random effects are causing marginal predictions to ",
        "underestimate the arithmetic mean on the response scale (Jensen's inequality): ",
        "exp(beta) underestimates the true mean when random effect variance is large. ",
        "Consider re-running with bias_adjust = 'manual' or bias_adjust = 'tmb' to apply ",
        "a lognormal bias correction and confirm the source of the bias. ",
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
