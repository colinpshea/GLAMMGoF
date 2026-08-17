#' Parametric-bootstrap predictions from GLMs and (G)LMMs
#'
#' Generate response-scale predictions with uncertainty intervals via
#' parametric bootstrap of the fixed-effect coefficients. Handles model
#' offsets (formula-side and argument-side) and supports an optional
#' analytical Jensen's inequality bias correction via `bias_adjust = "manual"`.
#'
#' @param mod A fitted model of a supported class (glmmTMB, lme4::glmerMod,
#'   lme4::lmerMod, lme4::glmer.nb, MASS::glm.nb, stats::glm, stats::lm).
#'   GAMs are not currently supported.
#' @param newdata Optional data frame at which to predict. If `NULL`, an
#'   auto-grid is constructed from the fixed-effect variables (offset
#'   variables are excluded and handled via the `offset` argument or
#'   defaults; see Details).
#' @param at Optional named list of values at which to evaluate specific
#'   variables in the auto-grid (e.g., `at = list(x = c(-1, 0, 1))`).
#' @param offset Optional numeric scalar or vector on the link scale
#'   (e.g., `log(effort)` for log-link models). When supplied, this value
#'   is added to the linear predictor for all prediction rows and takes
#'   precedence over any offset variables present in `newdata` or in the
#'   fitted model. Scalars are recycled to `nrow(newdata)`. See Details.
#' @param numeric_default,numeric_length Controls for auto-grid construction
#'   over numeric predictors.
#' @param n_sim Number of parametric bootstrap draws (default 5000).
#' @param alpha CI level (default 0.05 → 95% intervals).
#' @param bias_adjust One of `"none"` (default) or `"manual"`; `"manual"`
#'   multiplies each prediction by `jensen_correct(mod)` to correct for
#'   Jensen's inequality bias from the random-effect variance components.
#'   Has no effect on logit-link models (warning issued).
#' @param seed Optional integer seed for reproducibility.
#'
#' @details
#' **Offset handling.** If the fitted model contains an offset (either as
#' `offset()` in the formula or as an `offset =` argument), `boot_predict()`
#' will resolve it in the following order of precedence:
#'
#' 1. If the user supplies the `offset` argument, that value is used.
#' 2. Else if `newdata` contains the offset variable(s), the offset
#'    expression is evaluated on `newdata`.
#' 3. Else offset defaults to 0 (i.e., predictions at the per-unit-exposure
#'    reference), and a warning is issued.
#'
#' Auto-generated grids exclude offset variables; use the `offset` argument
#' to specify the reference exposure at which grid predictions should be made.
#'
#' **Output columns.** `boot_mean` is the mean of `exp(Xβ_j)` across draws
#' and includes a small implicit Jensen adjustment from parameter
#' uncertainty. `boot_median` approximates the plug-in point prediction
#' `exp(Xβ̂)`. For well-identified models these agree closely; substantial
#' divergence flags high parameter uncertainty (e.g., small n,
#' near-separation, sparse designs). Both are useful — the mean is the
#' correct marginal quantity when integrating over parameter uncertainty;
#' the median is closer to what most users mean by "the model's prediction."
#'
#' @return A data frame with `newdata` columns plus `boot_mean`,
#'   `boot_median`, `boot_se`, `boot_lwr`, `boot_upr`. The raw draws matrix
#'   and metadata are returned as attributes.
#'
#' @export
boot_predict <- function(mod,
                         newdata = NULL,
                         at = NULL,
                         offset = NULL,
                         numeric_default = "mean",
                         numeric_length = 20,
                         n_sim = 5000,
                         alpha = 0.05,
                         bias_adjust = c("none", "manual"),
                         seed = NULL) {

    bias_adjust <- match.arg(bias_adjust)
    .bp_reject_unsupported(mod)
    if (!is.null(seed)) set.seed(seed)

    bhat     <- .bp_fixef(mod)
    Vhat     <- .bp_vcov_fixef(mod)
    link     <- .bp_link(mod)
    cond_frm <- .bp_formula(mod)
    zi       <- .bp_zi_info(mod)

    if (!link %in% c("log", "logit", "identity")) {
        stop("boot_predict() currently supports log, logit, and identity ",
             "links only. Got link = '", link, "'.", call. = FALSE)
    }

    ## ---- Identify offset structure (used by auto-grid + offset resolver) ----
    tt_cond     <- stats::terms(cond_frm)
    off_idx     <- attr(tt_cond, "offset")
    off_arg     <- if (!is.null(mod$call) && !is.null(mod$call$offset)) mod$call$offset else NULL
    off_vars    <- character(0)
    if (!is.null(off_idx)) {
        off_expr <- attr(tt_cond, "variables")[[off_idx + 1]]
        off_vars <- c(off_vars, all.vars(off_expr))
    }
    if (!is.null(off_arg)) {
        off_vars <- c(off_vars, all.vars(off_arg))
    }
    off_vars         <- unique(off_vars)
    has_model_offset <- !is.null(off_idx) || !is.null(off_arg)

    expand_numeric <- function(x, spec, len) {
        if (is.function(spec)) return(spec(x))
        switch(spec,
               mean      = mean(x, na.rm = TRUE),
               median    = stats::median(x, na.rm = TRUE),
               range     = seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = len),
               quartiles = as.numeric(stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)),
               stop("Unknown numeric_default spec: '", spec, "'", call. = FALSE))
    }

    model_frame <- stats::model.frame(mod)

    if (is.null(newdata)) {
        cond_vars  <- all.vars(cond_frm)[-1]
        zi_vars    <- if (!is.null(zi$formula)) all.vars(zi$formula) else character(0)
        fixed_vars <- unique(c(cond_vars, zi_vars))
        ## Exclude offset variables from the auto-grid; they're handled separately
        fixed_vars <- setdiff(fixed_vars, off_vars)

        num_default_for <- function(v) {
            if (is.list(numeric_default)) {
                if (v %in% names(numeric_default)) numeric_default[[v]] else "mean"
            } else {
                numeric_default
            }
        }
        vals <- lapply(fixed_vars, function(v) {
            if (!is.null(at) && v %in% names(at)) {
                at[[v]]
            } else {
                x <- model_frame[[v]]
                if      (is.factor(x))    levels(x)
                else if (is.numeric(x))   expand_numeric(x, num_default_for(v), numeric_length)
                else if (is.character(x)) unique(x)
                else stop("Unsupported column type for '", v, "'", call. = FALSE)
            }
        })
        names(vals) <- fixed_vars
        newdata <- do.call(expand.grid, c(vals, list(stringsAsFactors = FALSE)))
        if (nrow(newdata) > 10000) {
            stop("Auto-generated grid has ", nrow(newdata), " rows. ",
                 "Supply `newdata` explicitly, narrow via `at`, or reduce ",
                 "numeric_length.", call. = FALSE)
        }
    }

    ## Coerce factor columns to match model frame levels
    for (col in intersect(names(newdata), names(model_frame))) {
        if (is.factor(model_frame[[col]])) {
            newdata[[col]] <- factor(newdata[[col]], levels = levels(model_frame[[col]]))
        }
    }

    ## ---- Resolve offset vector for prediction ----
    ## Precedence: user-supplied `offset` > offset variables in newdata > 0
    off_vec <- rep(0, nrow(newdata))
    if (!is.null(offset)) {
        if (!is.numeric(offset)) {
            stop("`offset` must be numeric (on the link scale).", call. = FALSE)
        }
        off_vec <- rep_len(as.numeric(offset), nrow(newdata))
        if (!has_model_offset) {
            warning("`offset` argument supplied but the fitted model has no ",
                    "offset term; the value will still be added to the linear ",
                    "predictor. Confirm this is intentional.", call. = FALSE)
        }
    } else if (has_model_offset) {
        ## Only attempt eval if ALL offset variables are present in newdata.
        ## Restrict enclos to baseenv() so eval() cannot fall through to the
        ## user's global environment and silently pick up a stale variable.
        if (length(off_vars) > 0 && all(off_vars %in% names(newdata))) {
            off_expr <- if (!is.null(off_idx)) {
                attr(tt_cond, "variables")[[off_idx + 1]]
            } else {
                off_arg
            }
            ## Strip the offset() wrapper if present — we want the inner expression
            if (is.call(off_expr) && identical(off_expr[[1]], as.name("offset"))) {
                off_expr <- off_expr[[2]]
            }
            off_resolved <- tryCatch(
                eval(off_expr, envir = newdata, enclos = baseenv()),
                error = function(e) NULL
            )
            if (is.null(off_resolved) || length(off_resolved) != nrow(newdata)) {
                warning("Failed to resolve model offset from `newdata` ",
                        "(evaluation error or length mismatch). Using offset = 0.",
                        call. = FALSE)
            } else {
                off_vec <- as.numeric(off_resolved)
            }
        } else {
            warning("Model has an offset but offset variable(s) [",
                    paste(off_vars, collapse = ", "),
                    "] are missing from `newdata` and no `offset` argument was ",
                    "supplied. Predictions will use offset = 0 (per-unit ",
                    "reference). Supply `offset` or include offset variables in ",
                    "`newdata` to change this.", call. = FALSE)
        }
    }

    ## ---- Design matrices ----
    cond_rhs <- stats::reformulate(attr(stats::terms(cond_frm), "term.labels"),
                                   intercept = attr(stats::terms(cond_frm), "intercept"))
    X_cond <- stats::model.matrix(cond_rhs, data = newdata)

    X_zi <- NULL
    if (zi$mode == "covariate") {
        zi_rhs <- stats::reformulate(attr(stats::terms(zi$formula), "term.labels"),
                                     intercept = attr(stats::terms(zi$formula), "intercept"))
        X_zi <- stats::model.matrix(zi_rhs, data = newdata)
    }

    ## ---- Parametric bootstrap draws ----
    cond_draws <- MASS::mvrnorm(n_sim, mu = bhat, Sigma = Vhat)
    if (is.null(dim(cond_draws))) cond_draws <- matrix(cond_draws, ncol = 1)

    zi_draws <- NULL
    if (zi$mode != "none") {
        zi_draws <- MASS::mvrnorm(n_sim, mu = zi$fixef, Sigma = zi$vcov)
        if (is.null(dim(zi_draws))) zi_draws <- matrix(zi_draws, ncol = 1)
    }

    link_fun <- switch(link, log = exp, logit = stats::plogis, identity = identity)

    resp <- matrix(NA_real_, nrow = nrow(newdata), ncol = n_sim)
    for (j in seq_len(n_sim)) {
        eta_cond <- as.vector(X_cond %*% cond_draws[j, ]) + off_vec
        mu_cond  <- link_fun(eta_cond)
        p_zi <- switch(zi$mode,
                       none      = 0,
                       intercept = stats::plogis(zi_draws[j, 1]),
                       covariate = stats::plogis(as.vector(X_zi %*% zi_draws[j, ])))
        resp[, j] <- mu_cond * (1 - p_zi)
    }

    ## ---- Jensen correction (RE-variance based) ----
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

    ## ---- Assemble output ----
    out <- cbind(newdata,
                 boot_mean   = rowMeans(resp),
                 boot_median = apply(resp, 1, stats::median),
                 boot_se     = apply(resp, 1, stats::sd),
                 boot_lwr    = apply(resp, 1, stats::quantile, probs = alpha/2),
                 boot_upr    = apply(resp, 1, stats::quantile, probs = 1 - alpha/2))

    attr(out, "draws")         <- resp
    attr(out, "jensen_factor") <- jf
    attr(out, "bias_adjust")   <- bias_adjust
    attr(out, "link")          <- link
    attr(out, "n_sim")         <- n_sim
    attr(out, "alpha")         <- alpha
    attr(out, "offset")        <- off_vec
    out
}
