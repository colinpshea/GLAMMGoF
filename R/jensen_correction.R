#' Jensen's inequality correction factor(s) for log-scale (G)LM(M) predictions
#'
#' Computes the multiplicative retransformation-bias correction factor
#' \eqn{\exp(V/2)} for back-transformed marginal predictions from log-link and
#' natural-log-response models. Complements [jensen_correct()] by handling the
#' cases the scalar path cannot: glmmTMB models with a non-trivial
#' \code{dispformula} (where \code{sigma()} returns \code{NA} and the correction
#' is per-observation rather than scalar). The returned numeric is suitable for
#' passing to [boot_predict()] via the \code{correction_factor} argument, or
#' for multiplying externally against response-scale predictions.
#'
#' @section What determines the correction:
#' The variance term \eqn{V} depends on where the additive error enters
#' relative to the exponential back-transformation.
#'
#' \describe{
#'   \item{Natural-log-response models (\code{log(y) ~ .})}{Gaussian on the
#'     log scale. Both the residual and random-effect variances are additive
#'     on the log scale, inside \eqn{\exp()}, so
#'     \eqn{V = \sigma^2_{resid} + \sum_k \sigma^2_k}. When the model has a
#'     non-trivial dispformula, \eqn{\sigma^2_{resid}(x)} depends on the
#'     covariates and \code{jensen_correction()} returns a length-\code{nrow(newdata)}
#'     vector, extracting the per-row dispersion from
#'     \code{predict(model, type = "disp", newdata = newdata)}.}
#'   \item{Log-link GLMMs (Poisson, negative binomial, Tweedie, Gamma, or
#'     \code{gaussian(link = "log")}, with random effects)}{Residual/dispersion
#'     is additive on the response scale, \emph{outside} \eqn{\exp()}, so
#'     \eqn{V = \sum_k \sigma^2_k}. A non-trivial dispformula affects
#'     \eqn{Var(Y)} but not \eqn{E[Y]}: the correction stays a scalar
#'     RE-only factor. An informational message is emitted so callers know
#'     the dispformula was seen and deliberately not used.}
#'   \item{glmmTMB \code{family = lognormal(link = "log")}}{Parameterizes the
#'     linear predictor to target \eqn{E[Y]} on the response scale directly
#'     (residual variance is absorbed into the mean structure by the MLE),
#'     so \eqn{V = \sum_k \sigma^2_k}. As with log-link GLMMs, a non-trivial
#'     dispformula is data-scale SD and does not enter the correction.}
#' }
#'
#' @section Scalar vs vector return:
#' The return is a single positive number when the correction is constant
#' across rows: any model without a non-trivial dispformula, and any log-link
#' or lognormal-family model regardless of dispformula. The return is a
#' numeric vector of length \code{nrow(newdata)} only for the one case where
#' the correction genuinely varies across rows: a natural-log-response model
#' (\code{log(y) ~ .}, family Gaussian) with a non-trivial dispformula. In
#' that case \code{newdata} is required.
#'
#' @param model A fitted \code{glmmTMB}, \code{lme4} (\code{merMod}),
#'   \code{glm}, or \code{lm} model. For log-link corrections the model must
#'   use a log link; for lognormal corrections the response must be
#'   natural-log transformed on the formula LHS.
#' @param newdata Optional data frame. Required when \code{model} is a
#'   glmmTMB fit to a natural-log-transformed response with a non-trivial
#'   dispformula; ignored (with a message) otherwise. Must contain all
#'   covariates referenced in the dispformula.
#' @param include_re Logical. If \code{TRUE} (default) the correction
#'   includes \eqn{\sum_k \sigma^2_k} from all random-intercept blocks. If
#'   \code{FALSE}, only the residual/dispersion contribution is included --
#'   useful when RE variance is handled elsewhere, or when only the
#'   dispformula contribution is wanted as a diagnostic.
#'
#' @return A positive numeric of length 1 (scalar case) or length
#'   \code{nrow(newdata)} (natural-log-response with non-trivial dispformula).
#'
#' @seealso [jensen_correct()] for the scalar-only correction and prediction
#'   application in one call; [boot_predict()] for parametric-bootstrap
#'   response-scale predictions that accept this correction via the
#'   \code{correction_factor} argument.
#'
#' @examples
#' \dontrun{
#' # Scalar case: log(y) ~ Gaussian with a single-covariate dispformula
#' # collapses to the trivial dispformula case and returns a scalar
#' m0 <- glmmTMB(log(y) ~ x + (1 | g), family = gaussian, data = dat)
#' jensen_correction(m0)                     # exp((sigma^2 + tau^2) / 2)
#'
#' # Vector case: dispformula varies with covariates
#' m1 <- glmmTMB(log(y) ~ x + (1 | g), dispformula = ~ x,
#'               family = gaussian, data = dat)
#' cf <- jensen_correction(m1, newdata = pred_grid)   # length nrow(pred_grid)
#' bp <- boot_predict(m1, newdata = pred_grid, correction_factor = cf)
#'
#' # log-link count model: RE-only, dispformula ignored with a message
#' m2 <- glmmTMB(y ~ x + (1 | g), dispformula = ~ x,
#'               family = nbinom2, data = dat)
#' jensen_correction(m2)                     # scalar exp(tau^2 / 2)
#' }
#'
#' @importFrom stats formula model.matrix predict sigma family
#' @importFrom reformulas nobars
#' @importFrom glmmTMB fixef
#' @export
jensen_correction <- function(model, newdata = NULL, include_re = TRUE) {

  # -- Validate model class ------------------------------------------------
  if (!inherits(model, c("glmmTMB", "merMod", "glm", "lm")))
    stop("jensen_correction() supports glmmTMB, lme4 (merMod), glm, and lm ",
         "model objects. For other model types, extract the relevant ",
         "variances manually and compute exp(sum(variances) / 2).",
         call. = FALSE)

  # -- Classify the model into one of four paths ---------------------------
  # path_A     : log(y) ~ ..., residual variance INSIDE exp() -- enters correction
  # path_B     : family = lognormal, mean parameterization absorbs residual
  # log_link   : link = "log", residual is OUTSIDE exp() and does not enter
  # unsupported: identity link, no log-transform, or non-log link (logit, etc.)
  fam    <- tryCatch(stats::family(model), error = function(e) NULL)
  link   <- if (!is.null(fam)) fam$link   else "identity"
  family <- if (!is.null(fam)) fam$family else "gaussian"
  lr     <- .detect_log_response(model)

  if (lr$logged && !lr$natural)
    stop("Only natural-log response transforms are supported (found '",
         lr$base, "'). The back-transformation and its variance scaling ",
         "differ for other bases; refit on the natural-log scale or ",
         "compute the correction manually.", call. = FALSE)

  path <- if (lr$natural) {
    "path_A"
  } else if (identical(family, "lognormal") && inherits(model, "glmmTMB")) {
    "path_B"
  } else if (identical(link, "log")) {
    "log_link"
  } else if (identical(link, "identity") && identical(family, "gaussian")) {
    stop("Identity link with an untransformed Gaussian response: no ",
         "back-transformation bias exists and no correction is warranted. ",
         "If the response is a pre-computed log column (e.g. 'logy'), ",
         "refit as log(y) ~ . so the transformation is visible.",
         call. = FALSE)
  } else {
    stop("No log back-transformation is present (link = '", link, "', ",
         "family = '", family, "'). jensen_correction() applies to log-link ",
         "and natural-log-response models only.", call. = FALSE)
  }

  # -- Random-effect variance sum (with slope guard) -----------------------
  rev    <- .re_variance_sum(model)
  re_sum <- if (include_re) rev$re_var else 0
  if (include_re && rev$slopes)
    warning("Random slopes detected. jensen_correction() uses random-",
            "intercept variances only and ignores slope variance/",
            "covariance; the returned factor is an approximation. For ",
            "glmmTMB random-slope models a more exact correction is ",
            "available via predict(., do.bias.correct = TRUE).",
            call. = FALSE)

  # -- Dispformula detection (glmmTMB only) --------------------------------
  has_dispformula <- FALSE
  if (inherits(model, "glmmTMB")) {
    disp_frm <- tryCatch(stats::formula(model, component = "disp"),
                         error = function(e) NULL)
    if (!is.null(disp_frm))
      has_dispformula <- !identical(deparse(disp_frm), "~1")
  }

  # -- Paths where dispformula is IRRELEVANT to E[Y] -----------------------
  # log_link : residual is outside exp(); dispformula affects Var(Y) not E[Y]
  # path_B   : glmmTMB lognormal absorbs residual into the mean parameterization
  # For both, the correction is RE-only regardless of dispformula.
  if (path %in% c("log_link", "path_B")) {
    if (has_dispformula) {
      msg <- if (path == "log_link") {
        paste0("Non-trivial dispformula detected. For log-link models, ",
               "residual/dispersion variance is additive on the response ",
               "scale (outside exp()) and does not enter the Jensen ",
               "correction. Returning RE-only scalar correction.")
      } else {
        paste0("Non-trivial dispformula detected. For glmmTMB ",
               "family = lognormal, the linear predictor targets E[Y] on ",
               "the response scale directly and residual variance is ",
               "absorbed into the mean parameterization. Returning RE-only ",
               "scalar correction.")
      }
      message(msg)
    }
    if (path == "log_link" && rev$re_var == 0)
      message("No random effects found with a log link: the correction ",
              "factor is 1. Log-link GLMs (no random effects) require no ",
              "Jensen correction.")
    return(exp(re_sum / 2))
  }

  # -- Path A: natural-log response, residual variance enters the correction -
  # Two subcases:
  #   (a) trivial dispformula (or non-glmmTMB backend): sigma() returns a
  #       scalar; return scalar exp((sigma^2 + RE_sum) / 2).
  #   (b) non-trivial dispformula (glmmTMB only): sigma() returns NA; extract
  #       per-row sigma via predict(type = "disp"); return length-nrow(newdata)
  #       vector exp((sigma^2(x) + RE_sum) / 2).

  if (!has_dispformula) {
    resid_var <- tryCatch(stats::sigma(model)^2, error = function(e) NA_real_)
    if (!is.finite(resid_var))
      stop("Could not extract a scalar residual variance via sigma() for ",
           "this natural-log-response model. If this is a glmmTMB fit with ",
           "a non-trivial dispformula, supply newdata so per-row dispersion ",
           "can be extracted via predict(type = 'disp').", call. = FALSE)
    return(exp((resid_var + re_sum) / 2))
  }

  # Dispformula subcase -- newdata is required
  if (is.null(newdata))
    stop("Model has a non-trivial dispformula; the residual variance ",
         "depends on covariates and the Jensen correction is a per-row ",
         "vector, not a scalar. Supply `newdata` at the covariates where ",
         "the correction should be evaluated.", call. = FALSE)

  # Extract per-row sigma^2 from the disp linear predictor.
  #
  # For glmmTMB gaussian family (verified against glmmTMB.cpp gaussian_family
  # case, line ~1014): eta_disp = X_disp %*% beta_disp; phi = exp(eta_disp);
  # dnorm(y, mu, phi) uses phi as SD, so sigma_i = exp(X_disp %*% beta_disp)
  # and sigma_i^2 = exp(2 * X_disp %*% beta_disp).
  #
  # We reconstruct this directly rather than calling predict(mod, type = "disp",
  # newdata = ...) because predict.glmmTMB validates newdata against the whole
  # model (including cond-model RE grouping variables), and callers of this
  # helper -- including boot_predict()'s auto-grid -- reasonably supply newdata
  # with only the dispformula covariates.
  disp_frm <- tryCatch(stats::formula(model, component = "disp"),
                       error = function(e) NULL)
  if (is.null(disp_frm))
    stop("Failed to extract dispformula from the model.", call. = FALSE)

  # Strip any random-effect bars (rare in dispformulas but possible); the
  # helper handles fixed-effect dispersion only. If bars are stripped, warn.
  disp_rhs <- reformulas::nobars(disp_frm)
  if (!identical(deparse(disp_frm), deparse(disp_rhs)))
    warning("Random effects in dispformula are ignored by jensen_correction(); ",
            "only the fixed-effect part of the dispersion model contributes to ",
            "the per-row correction.", call. = FALSE)

  beta_disp <- glmmTMB::fixef(model)$disp
  if (is.null(beta_disp) || length(beta_disp) == 0L)
    stop("Could not extract dispersion fixed effects via fixef(model)$disp. ",
         "This should not happen for a glmmTMB model with a non-trivial ",
         "dispformula -- please report with a reproducible example.",
         call. = FALSE)

  X_disp <- tryCatch(
    stats::model.matrix(disp_rhs, data = newdata),
    error = function(e) {
      stop("Failed to build dispersion design matrix from `newdata`: ",
           conditionMessage(e), "\nCheck that `newdata` contains all ",
           "covariates referenced in the dispformula (",
           deparse(disp_rhs), ") and that factor levels are consistent ",
           "with those used at fit time.", call. = FALSE)
    }
  )
  if (!identical(colnames(X_disp), names(beta_disp)))
    stop("Dispersion design matrix columns from `newdata` do not match the ",
         "dispersion coefficients from the fitted model. This usually means ",
         "`newdata` is missing factor levels that were present at fit time, ",
         "or that factor level ordering differs.\n",
         "  Design matrix columns: ", paste(colnames(X_disp), collapse = ", "),
         "\n  Model coefficient names: ", paste(names(beta_disp), collapse = ", "),
         call. = FALSE)

  log_sigma <- as.vector(X_disp %*% beta_disp)
  if (any(!is.finite(log_sigma)))
    stop("Non-finite log-sigma values from the dispersion linear predictor. ",
         "Check the fitted dispersion model for convergence problems or ",
         "extreme coefficient values.", call. = FALSE)

  resid_var_vec <- exp(2 * log_sigma)
  exp((resid_var_vec + re_sum) / 2)
}
