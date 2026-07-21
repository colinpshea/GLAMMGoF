# -----------------------------------------------------------------------------
# Internal helpers (not exported) shared by jensen_correct() and bias_precision()
# -----------------------------------------------------------------------------

#' Internal: was the response natural-log transformed in the model formula?
#'
#' Inspects the formula LHS only. Returns \code{natural = TRUE} only for a bare
#' \code{log(var)} (natural log of a single variable). \code{log2()},
#' \code{log10()}, \code{log1p()}, \code{log(var, base = )}, and compound
#' expressions such as \code{log(var + 1)} are reported as logged but
#' non-natural, so callers can reject them rather than mis-scale the
#' back-transformation. Detection deliberately keys on the response transform,
#' never on the family/link, which is what distinguishes a lognormal model
#' (\code{log(y) ~ .}, residual variance lives inside \code{exp()}) from a
#' Gaussian model with a log \emph{link} (\code{y ~ ., family = gaussian("log")},
#' residual variance is additive on the response scale and does not enter).
#'
#' @keywords internal
#' @noRd
.detect_log_response <- function(model) {
  lhs <- tryCatch(stats::formula(model)[[2]], error = function(e) NULL)
  out <- list(logged = FALSE, natural = FALSE, base = NA_character_)
  if (is.call(lhs)) {
    fn <- as.character(lhs[[1]])
    if (fn == "log" && length(lhs) == 2L && is.name(lhs[[2]])) {
      out <- list(logged = TRUE, natural = TRUE,  base = "e")
    } else if (fn %in% c("log", "log2", "log10", "log1p")) {
      out <- list(logged = TRUE, natural = FALSE, base = fn)
    } else {
      out <- list(logged = TRUE, natural = FALSE, base = fn)
    }
  }
  out
}

#' Internal: sum of random-effect intercept variances, with a random-slope flag
#'
#' Returns the sum of the intercept variances across all random-effect blocks
#' (the quantity that enters the scalar Jensen factor), and a flag indicating
#' whether any block carries random slopes. Both glmmTMB and lme4 backends are
#' reduced to intercept variances only, so the two agree exactly for
#' random-intercept models; slope variance/covariance is intentionally ignored
#' because the scalar factor \code{exp(sum(sigma^2)/2)} is not valid when the
#' linear-predictor variance depends on covariates. Models with no random
#' effects (\code{lm}, \code{glm}) return \code{re_var = 0}.
#'
#' @keywords internal
#' @noRd
.re_variance_sum <- function(model) {
  re_var <- 0
  slopes <- FALSE
  if (inherits(model, "glmmTMB")) {
    vc <- VarCorr(model)$cond
    if (length(vc) > 0) {
      re_var <- sum(vapply(vc, function(x) x[1, 1], numeric(1)))
      slopes <- any(vapply(vc, function(x) nrow(x) > 1L, logical(1)))
    }
  } else if (inherits(model, "merMod")) {
    vc <- VarCorr(model)
    if (length(vc) > 0) {
      re_var <- sum(vapply(vc, function(x) attr(x, "stddev")[1]^2, numeric(1)))
      slopes <- any(vapply(vc, function(x) length(attr(x, "stddev")) > 1L, logical(1)))
    }
  }
  list(re_var = re_var, slopes = slopes)
}

#' Internal: resolve model type and compute the retransformation-bias factor
#'
#' Central source of truth for the correction. \code{V} is the sum of the
#' variance components that sit \emph{inside} the exponential back-transform:
#' random-effect (intercept) variances always, plus the residual variance for a
#' lognormal (natural-log-transformed response) model. The returned
#' \code{factor} is \code{exp(V / 2)}.
#'
#' @keywords internal
#' @noRd
.jensen_factor <- function(model, type = c("auto", "log_link", "lognormal")) {
  type <- match.arg(type)

  fam    <- tryCatch(stats::family(model), error = function(e) NULL)
  link   <- if (!is.null(fam)) fam$link   else "identity"
  family <- if (!is.null(fam)) fam$family else "gaussian"
  lr     <- .detect_log_response(model)

  # --- resolve which correction applies ---
  if (type == "auto") {
    if (lr$natural) {
      resolved <- "lognormal"
    } else if (lr$logged) {
      stop("Only natural-log response transforms are supported (found '",
           lr$base, "'). The back-transformation and its variance scaling ",
           "differ for other bases; refit on the natural-log scale or correct ",
           "manually.")
    } else if (link == "log") {
      resolved <- "log_link"                    # incl. gaussian(link = 'log')
    } else if (link == "identity" && family == "gaussian") {
      stop("Cannot determine the correction automatically: identity link and no ",
           "visible log-transform of the response.\n",
           "  * Ordinary linear/Gaussian model -> no Jensen correction is ",
           "needed; do not use jensen_correct().\n",
           "  * Response is a pre-computed log column (e.g. 'logy') -> pass ",
           "type = 'lognormal' so the residual variance is included.")
    } else {
      stop("No log back-transformation is present (link = '", link, "'). ",
           "jensen_correct() applies to log-link and lognormal models only.")
    }
  } else {
    resolved <- type
    if (resolved == "lognormal" && lr$logged && !lr$natural)
      stop("type = 'lognormal' requires a natural-log response scale, but the ",
           "formula uses '", lr$base, "'. Refit on the natural-log scale.")
  }

  # --- random-effect variance (+ slope guard) ---
  rev <- .re_variance_sum(model)
  if (rev$slopes)
    warning("Random slopes detected. The scalar factor exp(sum(sigma^2)/2) uses ",
            "intercept variances only and ignores slope variance/covariance, so ",
            "it is approximate. For random-slope models use bias_adjust = 'tmb' ",
            "(glmmTMB) or predict(., do.bias.correct = TRUE).")

  # --- residual variance (inside exp() only for lognormal) ---
  resid_var <- 0
  if (resolved == "lognormal") {
    resid_var <- tryCatch(stats::sigma(model)^2, error = function(e) NA_real_)
    if (is.na(resid_var))
      stop("Could not extract a residual standard deviation via sigma() for the ",
           "lognormal correction. Supply the correction manually.")
  }

  if (resolved == "log_link" && rev$re_var == 0)
    warning("No random effects found with a log link: the correction factor is ",
            "1. Log-link GLMs (no random effects) require no Jensen correction.")

  V <- rev$re_var + resid_var
  list(factor       = exp(V / 2),
       V            = V,
       re_var       = rev$re_var,
       resid_var    = resid_var,
       type         = resolved,
       is_lognormal = identical(resolved, "lognormal"),
       has_re       = rev$re_var > 0)
}


#' Retransformation bias correction for log-link and lognormal (G)LM(M) predictions
#'
#' Computes the multiplicative bias correction factor \eqn{\exp(V/2)} for
#' marginal predictions that are back-transformed to the response scale via
#' \code{exp()}, and optionally applies it to a vector of predictions, their
#' standard errors, or confidence interval bounds. The correction addresses the
#' systematic underestimation of the arithmetic mean that arises from Jensen's
#' inequality when back-transforming population-level (marginal) predictions.
#'
#' The variance term \eqn{V} depends on where the additive error enters relative
#' to the exponential back-transformation:
#'
#' \describe{
#'   \item{Log-link GLMM}{(Poisson, negative binomial, Tweedie, gamma, or
#'     \code{gaussian(link = "log")} with random effects): the random effects are
#'     additive on the log scale, inside \eqn{\exp()}, but the conditional
#'     dispersion is not. \eqn{V = \sum_k \sigma^2_k}, the sum of random-effect
#'     variances. \eqn{\exp(\mathbf{X}\boldsymbol{\beta})} is the geometric mean;
#'     the arithmetic mean requires the correction.}
#'   \item{Lognormal LM / LMM}{(a natural-log-transformed response,
#'     \code{log(y) ~ .}, Gaussian on the log scale): the residual error is also
#'     additive on the log scale, inside \eqn{\exp()}, so
#'     \eqn{V = \sigma^2_{resid} + \sum_k \sigma^2_k}. A lognormal model with
#'     random effects is therefore \emph{doubly} biased.}
#' }
#'
#' A Gaussian model with a log \emph{link} on an untransformed response
#' (\code{y ~ ., family = gaussian("log")}) is a deliberate exception: its
#' residual error is additive on the response scale, \emph{outside} \eqn{\exp()},
#' and so it receives \eqn{V = \sum_k \sigma^2_k} (random effects only, no
#' residual term) exactly like the other log-link models. Detection keys on the
#' response transform in the formula, never on the family, so this case is
#' handled correctly and automatically.
#'
#' When \code{scale = "response"} (the default), predictions and associated
#' quantities are assumed to already be on the response scale (i.e.,
#' back-transformed via \code{exp()}). The correction factor is applied as a
#' simple multiplication: \eqn{\hat{\mu}_{adj} = \hat{\mu} \times c}, where
#' \eqn{c = \exp(V/2)}. Standard errors and confidence interval bounds scale
#' linearly by the same factor.
#'
#' When \code{scale = "link"}, predictions are assumed to be on the log scale
#' (the linear predictor for a log-link model, or the natural-log response scale
#' for a lognormal model). The function back-transforms them via \code{exp()} and
#' applies the correction in one step: \eqn{\hat{\mu}_{adj} = \exp(\hat{\eta})
#' \times c}. Standard errors on the link scale are propagated to the response
#' scale via the delta method before applying the correction:
#' \eqn{SE_{adj} = \exp(\hat{\eta}) \times SE_{\eta} \times c}. Confidence
#' intervals on the link scale are back-transformed symmetrically and then
#' multiplied by \eqn{c}.
#'
#' @param model A fitted \code{glmmTMB}, \code{lme4} (\code{merMod}), \code{glm},
#'   or \code{lm} model object. For log-link corrections the model must use a log
#'   link and have at least one random effect; for lognormal corrections the
#'   response must be natural-log transformed (\code{log(y) ~ .}), with or
#'   without random effects.
#' @param predictions An optional numeric vector of predictions to correct.
#'   Interpreted as response-scale predictions when \code{scale = "response"}
#'   (default), or as log-scale linear predictor values when
#'   \code{scale = "link"}. If \code{NULL} (default), only the scalar correction
#'   factor is returned.
#' @param se An optional numeric vector of standard errors, on the same scale as
#'   \code{predictions}. When \code{scale = "link"}, delta-method propagation is
#'   used to convert to the response scale before applying the correction.
#' @param lwr An optional numeric vector of lower confidence interval bounds, on
#'   the same scale as \code{predictions}.
#' @param upr An optional numeric vector of upper confidence interval bounds, on
#'   the same scale as \code{predictions}.
#' @param scale Character string, either \code{"response"} (default) or
#'   \code{"link"}. Specifies the scale of the supplied \code{predictions},
#'   \code{se}, \code{lwr}, and \code{upr}. All returned adjusted quantities are
#'   on the response scale.
#' @param type Character string, one of \code{"auto"} (default),
#'   \code{"log_link"}, or \code{"lognormal"}. With \code{"auto"} the correction
#'   is resolved from the model: a natural-log-transformed response gives the
#'   lognormal correction (residual + random-effect variance); a log link on an
#'   untransformed response gives the log-link correction (random-effect variance
#'   only). \code{"auto"} stops with an informative error when the response has an
#'   identity link and is not visibly log-transformed, since that case cannot be
#'   distinguished from an ordinary linear model. Set \code{type = "lognormal"}
#'   explicitly when the response is a pre-computed log column (e.g. \code{logy})
#'   whose transform is invisible to the formula parser. Set
#'   \code{type = "log_link"} to force the random-effect-only correction.
#'
#' @return If \code{predictions = NULL}, returns the scalar correction factor
#'   \eqn{\exp(V/2)} as a single numeric value. If \code{predictions} is
#'   supplied, returns a named list with elements \code{correction},
#'   \code{predictions_adjusted}, and (when supplied) \code{se_adjusted},
#'   \code{lwr_adjusted}, and \code{upr_adjusted}, all on the response scale.
#'
#' @note For \code{glmmTMB} models with random slopes or complex covariance
#'   structures (e.g., spatial random effects), the analytical scalar correction
#'   based on diagonal intercept variances is an approximation and a warning is
#'   issued; \code{predict(model, do.bias.correct = TRUE)} provides a more exact
#'   correction via TMB's automatic differentiation. The correction is not
#'   appropriate for models with an identity link and an untransformed response
#'   (ordinary Gaussian models), for which no back-transformation bias exists,
#'   nor for logit-link models where the direction and magnitude of Jensen's
#'   inequality bias varies with the linear predictor.
#'
#'   Whether to apply the bias correction depends on the nature of the random
#'   effects and the purpose of the predictions. Random effects that represent
#'   real, persistent features of the system (site effects, spatial fields,
#'   individual heterogeneity) warrant correction when making response-scale
#'   predictions. Random effects that represent study-specific nuisance variation
#'   may not require correction when the goal is population-level generalization.
#'   See the package vignette for a decision framework.
#'
#' @references Thorson, J.T. and Kristensen, K. (2016) Implementing a generic
#'   method for bias correction in statistical models using random effects, with
#'   spatial and population dynamics examples. \emph{Fisheries Research}, 175,
#'   66--74.
#'
#' @examples
#' \dontrun{
#' # Log-link GLMM: random-effect variance only
#' m <- glmmTMB(y ~ x + (1 | site), data = dat, family = nbinom2)
#' jensen_correct(m)
#'
#' # Lognormal LMM: residual + random-effect variance
#' m_ln <- lme4::lmer(log(y) ~ x + (1 | site), data = dat)
#' jensen_correct(m_ln)                       # exp((sigma^2_resid + sigma^2_RE) / 2)
#'
#' # Lognormal model fit on a pre-computed log column
#' dat$logy <- log(dat$y)
#' m_col <- lm(logy ~ x, data = dat)
#' jensen_correct(m_col, type = "lognormal")  # formula can't see the transform
#'
#' # Response-scale workflow
#' preds <- predict(m, newdata = dat, type = "response", re.form = ~0)
#' se    <- predict(m, newdata = dat, type = "response",
#'                  re.form = ~0, se.fit = TRUE)$se.fit
#' result <- jensen_correct(m, predictions = preds, se = se, scale = "response")
#'
#' # Link-scale workflow (back-transform + correct in one step)
#' preds_link <- predict(m, newdata = dat, type = "link", re.form = ~0)
#' result_link <- jensen_correct(m, predictions = preds_link, scale = "link")
#' }
#'
#' @importFrom nlme VarCorr
#' @importFrom stats family formula sigma
#' @export
jensen_correct <- function(model, predictions = NULL,
                           se = NULL, lwr = NULL, upr = NULL,
                           scale = c("response", "link"),
                           type  = c("auto", "log_link", "lognormal")) {

  scale <- match.arg(scale)
  type  <- match.arg(type)

  # --- Validate model class ---
  if (!inherits(model, c("glmmTMB", "merMod", "glm", "lm")))
    stop("jensen_correct() supports glmmTMB, lme4 (merMod), glm, and lm model ",
         "objects. For other model types, extract the relevant variances ",
         "manually and compute exp(sum(variances) / 2).")

  # --- Resolve type and compute the correction factor (shared logic) ---
  jf         <- .jensen_factor(model, type = type)
  correction <- jf$factor

  # --- Return scalar if no predictions supplied ---
  if (is.null(predictions)) return(correction)

  # --- Apply correction based on scale ---
  if (scale == "response") {

    # All inputs already on response scale -- multiply by correction
    out <- list(correction           = correction,
                predictions_adjusted = predictions * correction)
    if (!is.null(se))  out$se_adjusted  <- se  * correction
    if (!is.null(lwr)) out$lwr_adjusted <- lwr * correction
    if (!is.null(upr)) out$upr_adjusted <- upr * correction

  } else {

    # Link (log) scale -- backtransform via exp() then apply correction
    preds_resp <- exp(predictions) * correction

    out <- list(correction           = correction,
                predictions_adjusted = preds_resp)

    if (!is.null(se)) {
      # Delta method: SE_response = exp(eta) * SE_link, then multiply by c
      out$se_adjusted <- exp(predictions) * se * correction
    }

    if (!is.null(lwr)) out$lwr_adjusted <- exp(lwr) * correction
    if (!is.null(upr)) out$upr_adjusted <- exp(upr) * correction
  }

  out
}
