# GLAMMGoF 1.2.0

## New features

* Bumped to 1.2.0

* `jensen_correct()` now supports lognormal models (a natural-log-transformed
  response, `log(y) ~ .`), applying the combined correction
  `exp((sigma^2_residual + sum(sigma^2_RE)) / 2)`. Because both the residual
  error and the random effects are additive on the log scale, lognormal mixed
  models are subject to retransformation bias from both sources simultaneously;
  a correction built from `VarCorr()` alone recovers only part of the bias.

* `jensen_correct()` gains a `type` argument (`"auto"`, `"log_link"`,
  `"lognormal"`). The default, `"auto"`, resolves the appropriate correction by
  inspecting the response transformation in the model formula. Set
  `type = "lognormal"` explicitly when the response was log-transformed before
  fitting and stored as its own column (e.g. `logy ~ x`), where the
  transformation is invisible to the formula parser.

* `jensen_correct()` now accepts `lm` and `glm` objects in addition to `glmmTMB`
  and `lme4` fits, supporting lognormal models without random effects
  (correction `exp(sigma^2_residual / 2)`).

* `bias_precision()` now supports lognormal models. Predictions are
  back-transformed to the original response scale before metrics are computed,
  and `bias_adjust = "manual"` applies the combined residual + random-effect
  correction.

* Detection of log-transformed responses keys on the left-hand side of the model
  formula rather than the family or link, so `gaussian(link = "log")` (log link,
  untransformed response, residual additive on the response scale) is correctly
  distinguished from `log(y) ~ .` (log-transformed response, residual additive
  on the log scale) and receives no residual term.

## Behavior changes

* **Results for lognormal models will change substantially.** Previously,
  `bias_precision()` compared predictions from a `log(y) ~ .` model against the
  untransformed response, so predictions on the log scale were evaluated against
  observations on the original scale. The resulting RRMSE, RMAE, RMedAE, and
  RBIAS values were not meaningful. Predictions are now correctly
  back-transformed. Any previously reported metrics for log-transformed-response
  models should be regarded as invalid and recomputed.

* `jensen_correct()` now throws an informative error, rather than returning a
  correction factor, when passed a model with an identity link and an
  untransformed Gaussian response (e.g. `lmer(y ~ x + (1 | site))`). No
  retransformation bias exists for such models and the previously returned
  factor was not meaningful. This case cannot be distinguished automatically
  from a pre-logged response column; pass `type = "lognormal"` if the response
  is in fact on the log scale.

* For `lme4` models with random slopes, `jensen_correct()` now uses the random
  intercept variance only, and issues a warning. Previously all diagonal
  standard deviations, including slope terms, were summed. Results are unchanged
  for random-intercept models. The scalar correction is not valid for random
  slopes in either case, since the variance of the linear predictor then depends
  on the covariate; use `bias_adjust = "tmb"` or
  `predict(., do.bias.correct = TRUE)`.

* `bias_adjust = "manual"` in `bias_precision()` now warns when random slopes
  are detected, consistent with `jensen_correct()`.

* Responses transformed with a non-natural logarithm (`log10()`, `log2()`,
  `log1p()`, or `log(x, base = )`) are now rejected with an informative error in
  both functions, since each implies a different back-transformation.

* `bias_adjust = "tmb"` now throws an error for lognormal models. TMB's bias
  correction integrates over the random effects but does not include the
  residual retransformation term, so it would under-correct. Use
  `bias_adjust = "manual"`.

## Documentation

* The vignette gains a decision table distinguishing log-*link* models from
  log-transformed *responses*, and a worked example demonstrating that a
  correction built from random effect variance alone leaves substantial residual
  bias in a lognormal mixed model.

## Known limitations

* `conditional_predictions = TRUE` remains incompatible with
  `bias_adjust = "manual"`. For lognormal models the reasoning differs from the
  log-link case: conditional predictions absorb the random effect contribution
  but not the residual contribution, so a residual-only correction of
  `exp(sigma^2_residual / 2)` would still be required. This is not yet
  implemented.


# GLAMMGoF 1.1.4

* Bumped to 1.1.4

# GLAMMGoF 1.1.3

* Bumped to 1.1.3
* Added `correction_factor` argument to `bias_precision()` allowing users to supply a known correction factor directly when `bias_adjust = "manual"`, bypassing internal VarCorr() computation. Primarily useful in simulation contexts where the true random effect variance is known and RE variance estimates from training subsets may be unstable at high sigma values. 
* Added argument validation checks for swapped testModel/testData arguments in bias_precision() and brier_auc().

# GLAMMGoF 1.1.2

* Added `bias_adjust = "manual"` support for `lme4` model objects (`glmerMod`, `lmerMod`) using `attr(VarCorr(model), "stddev")^2` for RE variance extraction
* Added `conditional_predictions = TRUE` support for `lme4` model objects
* Added informative `stop()` for `bias_adjust = "tmb"` with `lme4` models, pointing users to `"manual"` or refitting in `glmmTMB`
* Added `jensen_correct()` as an exported utility function for standalone lognormal bias correction of predictions, standard errors, and confidence intervals on either the response or link scale; supports both `glmmTMB` and `lme4` model objects
* Updated `@param`, `@note`, and vignette documentation to reflect `lme4` support for `bias_adjust = "manual"` and `conditional_predictions`

# GLAMMGoF 1.1.0

* Fixed: `glmmTMB` backend now correctly uses `re.form = ~0` for marginal predictions, consistent with `lme4` and `mgcv` backends and documented behavior
* Fixed: `nlme` added to `Imports` in DESCRIPTION
* Added: `bias_adjust` argument (`"none"`, `"manual"`, `"tmb"`) to `bias_precision()` and `brier_auc()` for Jensen's inequality bias diagnosis and correction in `glmmTMB` models; `"manual"` applies analytical lognormal correction `exp(sigma^2/2)`; `"tmb"` uses TMB's automatic differentiation
* Added: `conditional_predictions` argument to `bias_precision()` and `brier_auc()` for within-group predictive accuracy assessment using estimated random effects
* Added: `verbose` argument to `bias_precision()` to suppress Jensen diagnostic message in simulation or sweep contexts
* Added: runtime diagnostic message when both in-sample and out-of-sample RBIAS are consistently negative (< -10%) in a `glmmTMB` log-link model with `bias_adjust = "none"`, suggesting Jensen's inequality as a potential source
* Added: informative `stop()` for `bias_adjust = "manual"` with `conditional_predictions = TRUE` to prevent double-correction
* Added: informative `stop()` for `bias_adjust = "manual"` in `brier_auc()` since the analytical lognormal correction is not valid for logit-link models
* Added: Thorson & Kristensen (2016) reference to both functions
* Added: comprehensive `@note` sections documenting Jensen's inequality, substantive vs nuisance random effects, random slopes caveat, spatial RE considerations, correction factor leakage, and recommended workflow
* Added: stepwise diagnostic workflow in `@note` and vignette
* Added: combination matrix table in vignette summarising all valid `bias_adjust` x `conditional_predictions` combinations
* Added: lognormal vs log-link GLMM parallel table in vignette
* Added: substantive vs nuisance random effects decision framework in vignette
* Added: recommended workflow section in vignette (marginal for inference; bias-corrected for predictions and figures)
* Updated: `bias_adjust = "tmb"` documentation clarified as a computational necessity for TMB's AD rather than a scientific choice for group-specific conditional predictions

# GLAMMGoF 1.0.8

* Fixed: `glmmTMB` predict call updated to use `re.form = ~0` explicitly for marginal predictions
* Added: `bias_adjust` argument (initial `TRUE`/`FALSE` version, later replaced in 1.1.0)
* Added: correction factor pre-computed from full-data `testModel` via `VarCorr()` for stability across resampling replicates

# GLAMMGoF 1.0.7

* Added: `group` argument for group-level resampling via shared internal `resample_split()` helper
* Added: DHARMa zero-inflation test via `testZI` argument (default `TRUE`)
* Added: full NA handling via complete-cases check on all model frame variables
* Added: `jensen_correct()` internal function (later exported in 1.1.2)

# GLAMMGoF 1.0.0

* Initial release
* `bias_precision()` for continuous and integer response models; returns RRMSE, RMAE, RMedAE, and RBIAS
* `brier_auc()` for binary response models; returns AUC, Brier score, and log loss
* Supports `glmmTMB`, `lme4`, `mgcv`, `MASS`, and `stats` model objects
* Holdout and bootstrap resampling methods
* Optional DHARMa residual diagnostics
* Example datasets `countData` and `logitData` with six fitted example models
