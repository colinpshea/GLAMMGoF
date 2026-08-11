#' Parametric bootstrap predictions with optional Jensen's inequality correction
#'
#' Produces population-level (marginal) predictions from a fitted model,
#' with confidence intervals derived from a parametric bootstrap over the
#' joint asymptotic distribution of the fixed-effect coefficients. When
#' the fitted model has variance components on the log scale, the
#' analytical Jensen's inequality correction can be applied by setting
#' `bias_adjust = "manual"`, delivering response-scale predictions that
#' are calibrated to the arithmetic mean of the response rather than the
#' geometric mean.
#'
#' @param mod A fitted model object. Supported classes:
#'   \code{\link[glmmTMB]{glmmTMB}}, \code{\link[lme4]{merMod}} (including
#'   \code{lmerMod}, \code{glmerMod}, and negative binomial fits from
#'   \code{glmer.nb}), \code{\link[MASS]{glm.nb}} (\code{negbin} objects),
#'   \code{\link[stats]{glm}}, and \code{\link[stats]{lm}}. \code{mgcv}
#'   GAMs and \code{pscl::zeroinfl} models are not currently supported;
#'   passing one produces an informative error.
#' @param newdata Optional data frame of covariate values at which to
#'   predict. If \code{NULL} (the default), a grid is auto-generated
#'   from the model's covariates using the rules described in Details.
#' @param at Optional named list of explicit covariate values to use in
#'   the auto-generated grid, e.g. \code{at = list(x = c(0, 1, 2))}.
#'   Values passed via \code{at} override the defaults derived from
#'   \code{numeric_default}. Ignored when \code{newdata} is supplied.
#' @param numeric_default Controls how continuous covariates are expanded
#'   in the auto-generated grid. One of \code{"mean"} (single value at
#'   the mean; default), \code{"median"}, \code{"range"}
#'   (\code{numeric_length} equally spaced values from min to max),
#'   or \code{"quartiles"} (Q1, median, Q3). A user-supplied function
#'   taking a numeric vector and returning a vector may also be passed.
#'   A named list of specifications may be used for per-covariate
#'   control, e.g. \code{numeric_default = list(x = "range", z = "mean")};
#'   any covariate not named in the list falls back to \code{"mean"}.
#'   Ignored when \code{newdata} is supplied.
#' @param numeric_length Integer giving the length of the sequence used
#'   when \code{numeric_default = "range"}. Defaults to 20.
#' @param n_sim Integer giving the number of parametric bootstrap draws.
#'   Defaults to 5000. Larger values reduce Monte Carlo noise in the
#'   summary statistics at the cost of runtime; 1000 is usually
#'   sufficient for stable point estimates, 5000 or more for stable
#'   tail quantiles.
#' @param alpha Numeric in (0, 1) giving the two-sided confidence level.
#'   Defaults to 0.05, producing 95\% intervals via the
#'   \code{alpha / 2} and \code{1 - alpha / 2} quantiles of the
#'   bootstrap distribution.
#' @param bias_adjust Character. One of \code{"none"} (default) or
#'   \code{"manual"}. When \code{"manual"}, the analytical Jensen's
#'   inequality correction from \code{\link{jensen_correct}} is applied
#'   to the bootstrap predictions. See Details.
#' @param seed Optional integer passed to \code{set.seed()} for
#'   reproducibility. When \code{NULL} (the default), the current RNG
#'   state is used.
#'
#' @details
#' \strong{Prediction target.} Predictions are population-level
#' (marginal), corresponding to \code{re.form = ~0} in the standard
#' \code{predict()} interface: random effects are set to zero and the
#' output represents the expected response for a new observation whose
#' group-level context is unknown. This is the appropriate target for
#' predicting to new sites, new years, or new individuals, and it is
#' the target to which the Jensen correction applies.
#'
#' \strong{Parametric bootstrap.} Fixed-effect coefficient uncertainty
#' is propagated by drawing \code{n_sim} samples from the joint
#' asymptotic distribution \code{MvNormal(coef(mod), vcov(mod))} on the
#' link scale, transforming to the response scale via the model's
#' inverse-link function, and summarizing across draws. This captures
#' full covariance among fixed effects, including correlations that
#' delta-method approximations may not represent accurately for
#' non-linear back-transformations.
#'
#' \strong{Variance-component uncertainty is treated as plug-in.} The
#' Jensen correction factor \eqn{\exp(\Sigma \sigma^2 / 2)} uses
#' \eqn{\hat{\sigma}^2} at its fitted point estimate. Uncertainty in
#' \eqn{\hat{\sigma}^2} is not propagated through the bootstrap. This
#' matches the plug-in convention used elsewhere in the ecosystem
#' (\code{\link{jensen_correct}}, \code{ggeffects::predict_response}
#' with \code{bias_correction = TRUE}) and is a reasonable
#' approximation at ecological RE variances. For a discussion of when
#' plug-in behavior may under- or over-correct, see the tutorial
#' vignette section on residual bias sources.
#'
#' \strong{Jensen correction scope.} Setting \code{bias_adjust = "manual"}
#' invokes \code{\link{jensen_correct}}, which resolves the appropriate
#' correction factor from the fitted model. For log-link mixed models
#' (\code{glmmTMB} or \code{lme4} with a log link and random effects),
#' this is \eqn{\exp(\Sigma \sigma^2_u / 2)}. For log-transformed
#' responses (\code{lm(log(y) ~ x)} or the mixed-model analogue), this
#' is \eqn{\exp((\sigma^2_\varepsilon + \Sigma \sigma^2_u) / 2)}. For
#' fixed-effect log-link GLMs and identity-response models with no log
#' transform, no correction is needed and the factor is 1; a warning
#' is emitted in that case to alert the user that \code{bias_adjust}
#' had no effect. For \code{logit}-link models, no closed-form scalar
#' correction exists and a warning is emitted.
#'
#' \strong{Zero-inflated models.} Zero-inflation is supported for
#' \code{glmmTMB} models with either intercept-only or covariate-
#' dependent \code{ziformula}. Predictions are on the marginal scale
#' \eqn{E[Y] = (1 - p_{\mathrm{zi}}) \cdot E[Y | Y > 0]}, and the
#' Jensen correction (when requested) is applied to the whole
#' marginal expectation, which is correct when random effects appear
#' only in the conditional (count) submodel. Random effects on the
#' zero-inflation submodel are not currently supported.
#'
#' \strong{Grid auto-build.} When \code{newdata} is not supplied, a
#' prediction grid is constructed by crossing the levels of each
#' categorical predictor with the \code{numeric_default} expansion of
#' each continuous predictor. Grids exceeding 10,000 rows trigger an
#' error to prevent unintentionally massive bootstrap runs; supply
#' \code{newdata} explicitly, narrow the grid via \code{at}, or reduce
#' \code{numeric_length} to proceed.
#'
#' @return A data frame containing the columns of \code{newdata} (or
#'   the auto-generated grid) augmented with:
#'   \describe{
#'     \item{\code{boot_mean}}{Mean of the bootstrap distribution at
#'       each prediction row.}
#'     \item{\code{boot_se}}{Standard deviation of the bootstrap
#'       distribution at each prediction row.}
#'     \item{\code{boot_lwr}}{Lower confidence bound
#'       (\code{alpha / 2} quantile).}
#'     \item{\code{boot_upr}}{Upper confidence bound
#'       (\code{1 - alpha / 2} quantile).}
#'   }
#'   The returned object also carries the following attributes:
#'   \code{draws} (the full \code{n_rows} \eqn{\times} \code{n_sim}
#'   bootstrap matrix), \code{jensen_factor} (the scalar correction
#'   applied, or 1 if none), \code{bias_adjust}, \code{link},
#'   \code{n_sim}, and \code{alpha}.
#'
#' @section Reproducibility:
#' Setting \code{seed} produces fully reproducible output. Two calls
#' with the same \code{mod}, \code{newdata} (or grid inputs), and
#' \code{seed} return identical predictions.
#'
#' @seealso
#' \code{\link{jensen_correct}} for the correction-factor calculation
#' used internally, and for applying the correction to predictions
#' produced by other packages. \code{\link{bias_precision}} for
#' resampling-based validation that diagnoses Jensen bias in
#' out-of-sample predictions. The tutorial vignette
#' (\code{vignette("GLAMMGoF_tutorial")}) works through the full
#' diagnose-correct-predict workflow.
#'
#' @note
#' Support is currently limited to fixed effects in the parametric
#' bootstrap. GAMs (\code{mgcv}) are not supported because variance
#' components from smooth terms require a different extraction
#' pathway (\code{gam.vcomp()}) that is not yet integrated with
#' \code{\link{jensen_correct}}. Bayesian model objects
#' (\code{brmsfit}, \code{stanreg}) are not supported; use the
#' package's native posterior-prediction functions and apply
#' \code{\link{jensen_correct}} to the response-scale draws.
#'
#' @examples
#' \dontrun{
#' # Fixed-effect log-link GLM: no correction needed, factor = 1
#' fit_glm <- glmmTMB::glmmTMB(y ~ Season + Temp,
#'                             family = glmmTMB::nbinom2,
#'                             data = countData)
#' boot_predict(fit_glm, bias_adjust = "none", seed = 123)
#'
#' # Log-link GLMM with one random effect: correction applied
#' fit_glmm <- glmmTMB::glmmTMB(y ~ Season + Temp + (1 | Site),
#'                              family = glmmTMB::nbinom2,
#'                              data = countData)
#' out <- boot_predict(fit_glmm, bias_adjust = "manual", seed = 123)
#' attr(out, "jensen_factor")   # correction factor applied
#'
#' # Custom grid: predict at specific Temp values across all Seasons
#' boot_predict(fit_glmm,
#'              at = list(Temp = c(10, 15, 20)),
#'              bias_adjust = "manual",
#'              seed = 123)
#'
#' # Full workflow: diagnose bias, then generate corrected predictions
#' bias_precision(fit_glmm, testData = countData,
#'                bias_adjust = "none")     # RBIAS negative
#' bias_precision(fit_glmm, testData = countData,
#'                bias_adjust = "manual")   # RBIAS near zero
#' boot_predict(fit_glmm, bias_adjust = "manual", seed = 123)
#' }
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats coef vcov model.frame model.matrix reformulate terms formula plogis median quantile sd family
#' @importFrom reformulas nobars
#' @export
boot_predict <- function(mod,
                         newdata          = NULL,
                         at               = NULL,
                         numeric_default  = "mean",
                         numeric_length   = 20,
                         n_sim            = 5000,
                         alpha            = 0.05,
                         bias_adjust      = c("none", "manual"),
                         seed             = NULL) {

  bias_adjust <- match.arg(bias_adjust)

  # Fail early on unsupported model classes with an informative
  # message that points to the right alternative. See dispatcher
  # file for the full list of rejected classes.
  .bp_reject_unsupported(mod)

  if (!is.null(seed)) set.seed(seed)

  # ── Extract backend-agnostic quantities via dispatchers ────
  bhat     <- .bp_fixef(mod)          # named numeric vector
  Vhat     <- .bp_vcov_fixef(mod)     # square matrix, same order
  link     <- .bp_link(mod)           # "log" | "logit" | "identity"
  cond_frm <- .bp_formula(mod)        # RHS-usable formula, no bars
  zi       <- .bp_zi_info(mod)        # list: mode, fixef, vcov, formula

  if (!link %in% c("log", "logit", "identity")) {
    stop("boot_predict() currently supports log, logit, and identity ",
         "links only. Got link = '", link, "'.", call. = FALSE)
  }

  # ── Helper: expand a single numeric covariate to a vector ──
  # Kept internal to boot_predict so the numeric_default API is
  # self-contained; can be extracted if reused elsewhere.
  expand_numeric <- function(x, spec, len) {
    if (is.function(spec)) return(spec(x))
    switch(
      spec,
      mean      = mean(x, na.rm = TRUE),
      median    = stats::median(x, na.rm = TRUE),
      range     = seq(min(x, na.rm = TRUE),
                      max(x, na.rm = TRUE),
                      length.out = len),
      quartiles = as.numeric(stats::quantile(
        x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)),
      stop("Unknown numeric_default spec: '", spec, "'", call. = FALSE)
    )
  }

  # ── Auto-build grid if newdata not supplied ─────────────────
  model_frame <- stats::model.frame(mod)

  if (is.null(newdata)) {
    cond_vars <- all.vars(cond_frm)[-1]     # drop LHS (response)
    zi_vars   <- if (!is.null(zi$formula)) all.vars(zi$formula) else character(0)
    fixed_vars <- unique(c(cond_vars, zi_vars))

    # Normalize numeric_default: allow either a scalar spec applied
    # to all numeric covariates, or a per-variable named list.
    num_default_for <- function(v) {
      if (is.list(numeric_default)) {
        if (v %in% names(numeric_default)) numeric_default[[v]] else "mean"
      } else {
        numeric_default
      }
    }

    vals <- lapply(fixed_vars, function(v) {
      if (!is.null(at) && v %in% names(at)) {
        # User specified values explicitly — use them as-is
        at[[v]]
      } else {
        x <- model_frame[[v]]
        if (is.factor(x)) {
          levels(x)
        } else if (is.numeric(x)) {
          expand_numeric(x, num_default_for(v), numeric_length)
        } else if (is.character(x)) {
          unique(x)
        } else {
          stop("Unsupported column type for '", v, "'", call. = FALSE)
        }
      }
    })
    names(vals) <- fixed_vars
    newdata <- do.call(expand.grid,
                       c(vals, list(stringsAsFactors = FALSE)))

    if (nrow(newdata) > 10000) {
      stop("Auto-generated grid has ", nrow(newdata), " rows. ",
           "Supply `newdata` explicitly, narrow via `at`, or reduce ",
           "numeric_length.", call. = FALSE)
    }
  }

  # Coerce factor columns to have model-original levels — avoids
  # "contrasts can be applied only to factors with 2 or more levels".
  for (col in intersect(names(newdata), names(model_frame))) {
    if (is.factor(model_frame[[col]])) {
      newdata[[col]] <- factor(newdata[[col]],
                               levels = levels(model_frame[[col]]))
    }
  }

  # ── Design matrices ────────────────────────────────────────
  # Rebuild an intercept-preserving RHS for the conditional model.
  cond_rhs  <- stats::reformulate(
    attr(stats::terms(cond_frm), "term.labels"),
    intercept = attr(stats::terms(cond_frm), "intercept")
  )
  X_cond <- stats::model.matrix(cond_rhs, data = newdata)

  X_zi <- NULL
  if (zi$mode == "covariate") {
    zi_rhs <- stats::reformulate(
      attr(stats::terms(zi$formula), "term.labels"),
      intercept = attr(stats::terms(zi$formula), "intercept")
    )
    X_zi <- stats::model.matrix(zi_rhs, data = newdata)
  }

  # ── Parameter draws ────────────────────────────────────────
  # Joint asymptotic normal at the fitted point estimate. This
  # captures fixed-effect uncertainty (including any correlation
  # between fixed effects) but treats variance components as
  # plug-in — see file header for rationale.
  cond_draws <- MASS::mvrnorm(n_sim, mu = bhat, Sigma = Vhat)
  if (is.null(dim(cond_draws))) cond_draws <- matrix(cond_draws, ncol = 1)

  zi_draws <- NULL
  if (zi$mode != "none") {
    zi_draws <- MASS::mvrnorm(n_sim, mu = zi$fixef, Sigma = zi$vcov)
    if (is.null(dim(zi_draws))) zi_draws <- matrix(zi_draws, ncol = 1)
  }

  # ── Response-scale transform ───────────────────────────────
  link_fun <- switch(link,
                     log      = exp,
                     logit    = stats::plogis,
                     identity = identity)

  # ── Simulation loop ────────────────────────────────────────
  resp <- matrix(NA_real_, nrow = nrow(newdata), ncol = n_sim)
  for (j in seq_len(n_sim)) {
    eta_cond <- as.vector(X_cond %*% cond_draws[j, ])
    mu_cond  <- link_fun(eta_cond)

    p_zi <- switch(zi$mode,
                   none      = 0,
                   intercept = stats::plogis(zi_draws[j, 1]),
                   covariate = stats::plogis(as.vector(X_zi %*% zi_draws[j, ])))

    resp[, j] <- mu_cond * (1 - p_zi)
  }

  # ── Jensen's inequality correction ─────────────────────────
  # Applies for log-link and log-transformed-response models with
  # variance components on the log scale. jensen_correct() returns
  # 1 for models where no correction is warranted (fixed-effects
  # log-link GLMs, logit models, identity-response models with no
  # log transform), so the multiplication is safe as a no-op in
  # those cases. Warning is only raised for the affirmatively-
  # inapplicable case: non-log links where the user explicitly
  # requested correction.
  jf <- 1
  if (bias_adjust == "manual") {
    if (link == "logit") {
      warning("bias_adjust = 'manual' has no effect for logit-link models. ",
              "Jensen's inequality correction only applies where predictions ",
              "are back-transformed via the exponential.", call. = FALSE)
    } else {
      jf   <- jensen_correct(mod)
      resp <- resp * jf
    }
  }

  # ── Summarize and return ───────────────────────────────────
  out <- cbind(newdata,
               boot_mean = rowMeans(resp),
               boot_se   = apply(resp, 1, stats::sd),
               boot_lwr  = apply(resp, 1, stats::quantile, probs = alpha / 2),
               boot_upr  = apply(resp, 1, stats::quantile, probs = 1 - alpha / 2))

  attr(out, "draws")         <- resp
  attr(out, "jensen_factor") <- jf
  attr(out, "bias_adjust")   <- bias_adjust
  attr(out, "link")          <- link
  attr(out, "n_sim")         <- n_sim
  attr(out, "alpha")         <- alpha

  out
}
