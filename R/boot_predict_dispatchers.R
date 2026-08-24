###############################################################
# boot_predict() backend dispatchers
#
# S3 methods that extract a common interface from each supported
# backend. Backend classes: glmmTMB, merMod (lme4), negbin (MASS),
# glm, lm. Each dispatcher returns a value with a shape that the
# main boot_predict() body can consume without knowing the backend.
#
# Coverage matches jensen_correct(): glmmTMB, lme4, lm/glm, and
# MASS::glm.nb (via the negbin class, which also inherits glm).
# GAMs are intentionally excluded — mgcv's variance-component
# handling via gam.vcomp() differs enough that a coherent Jensen
# correction requires its own design.
###############################################################

# ── Fixed-effect coefficients ──────────────────────────────────
# Returns a named numeric vector of conditional-model fixed effects.
# glmmTMB has a $cond component because of ZI; other backends are
# single-component.

.bp_fixef <- function(mod) UseMethod(".bp_fixef")

.bp_fixef.glmmTMB <- function(mod) glmmTMB::fixef(mod)$cond
.bp_fixef.merMod  <- function(mod) lme4::fixef(mod)
.bp_fixef.default <- function(mod) stats::coef(mod)   # lm, glm, negbin

# ── Fixed-effect covariance matrix ─────────────────────────────
# Returns a square matrix, same order as .bp_fixef(mod).
# This is the joint asymptotic covariance of the conditional-model
# fixed effects only — variance components are NOT included, since
# the Jensen correction handles sigma-hat as a plug-in scalar.

.bp_vcov_fixef <- function(mod) UseMethod(".bp_vcov_fixef")

.bp_vcov_fixef.glmmTMB <- function(mod) stats::vcov(mod)$cond
.bp_vcov_fixef.merMod  <- function(mod) as.matrix(stats::vcov(mod))
.bp_vcov_fixef.default <- function(mod) stats::vcov(mod)

# ── Link function name ─────────────────────────────────────────
# Returns a character scalar. boot_predict() then dispatches
# on it for the response-scale transform.

.bp_link <- function(mod) UseMethod(".bp_link")

.bp_link.glmmTMB <- function(mod) mod$modelInfo$family$link
.bp_link.merMod  <- function(mod) {
  # lmer models are Gaussian identity; glmer/glmer.nb carry family()
  if (inherits(mod, "lmerMod")) "identity" else stats::family(mod)$link
}
.bp_link.negbin  <- function(mod) mod$family$link   # inherits glm but
                                                    # family() sometimes
                                                    # returns the wrong
                                                    # structure on negbin
.bp_link.glm     <- function(mod) stats::family(mod)$link
.bp_link.lm      <- function(mod) "identity"

# ── Conditional-model formula (RHS-usable, no bars) ────────────
# Returns a formula suitable for building a design matrix on newdata.
# RE terms stripped; ZI terms are separate (see .bp_zi_info).

.bp_formula <- function(mod) UseMethod(".bp_formula")

.bp_formula.glmmTMB <- function(mod) {
  reformulas::nobars(stats::formula(mod, component = "cond"))
}
.bp_formula.merMod  <- function(mod) reformulas::nobars(stats::formula(mod))
.bp_formula.default <- function(mod) stats::formula(mod)  # lm/glm/negbin

# ── Zero-inflation structure ───────────────────────────────────
# Returns a list with three elements:
#   mode      : "none" | "intercept" | "covariate"
#   fixef     : ZI fixed-effect coefficients (or NULL if mode = "none")
#   vcov      : ZI fixed-effect covariance   (or NULL if mode = "none")
#   formula   : ZI formula, no bars          (or NULL if mode != "covariate")
#
# Only glmmTMB supports ZI in the current backend set. All other
# backends return mode = "none" and boot_predict() skips ZI machinery.

.bp_zi_info <- function(mod) UseMethod(".bp_zi_info")

.bp_zi_info.glmmTMB <- function(mod) {
  zi_fixef <- glmmTMB::fixef(mod)$zi
  n_zi     <- length(zi_fixef)

  if (n_zi == 0) {
    return(list(mode = "none", fixef = NULL, vcov = NULL, formula = NULL))
  }

  mode <- if (n_zi == 1) "intercept" else "covariate"
  zi_vcov <- stats::vcov(mod)$zi

  zi_formula <- NULL
  if (mode == "covariate") {
    zi_formula <- reformulas::nobars(mod$modelInfo$allForm$ziformula)
  }

  list(mode = mode, fixef = zi_fixef, vcov = zi_vcov, formula = zi_formula)
}

.bp_zi_info.default <- function(mod) {
  list(mode = "none", fixef = NULL, vcov = NULL, formula = NULL)
}

# ── Explicit rejection of unsupported backends ─────────────────
# pscl::zeroinfl(), brms models, etc. Fail early with a message
# that points users to the supported set.

.bp_reject_unsupported <- function(mod) {
  if (inherits(mod, "zeroinfl")) {
    stop("boot_predict() does not support pscl::zeroinfl() models. ",
         "Refit using glmmTMB with a ziformula for zero-inflated support.",
         call. = FALSE)
  }
  if (inherits(mod, c("brmsfit", "stanfit", "stanreg"))) {
    stop("boot_predict() does not support Bayesian model objects. ",
         "If posterior predictions are drawn directly on the response scale (e.g. brms::posterior_epred()), ",
         "no bias correction is needed. ",
         "If link-scale posterior draws are used (e.g. brms::posterior_linpred()), ",
         "apply GLAMMGoF::jensen_correct() to the back-transformed response-scale summaries.",
         call. = FALSE)
  }
  if (inherits(mod, c("gam", "bam", "gamm"))) {
    stop("boot_predict() does not currently support mgcv models. ",
         "GAM variance components require separate handling via gam.vcomp() ",
         "that is not yet implemented in jensen_correct().",
         call. = FALSE)
  }
  invisible(NULL)
}

# ── Sanity-check convenience for interactive use ───────────────
# Not exported — just prints what the dispatchers see for a given
# model. Useful during development to confirm the seams are right
# before wiring boot_predict()'s body to them.

.bp_dispatch_summary <- function(mod) {
  .bp_reject_unsupported(mod)
  list(
    class     = class(mod)[1],
    fixef     = .bp_fixef(mod),
    vcov_dim  = dim(.bp_vcov_fixef(mod)),
    link      = .bp_link(mod),
    formula   = .bp_formula(mod),
    zi        = .bp_zi_info(mod)
  )
}
