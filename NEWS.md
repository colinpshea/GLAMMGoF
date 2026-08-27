GLAMMGoF 1.4.4

## Documentation

* Clarified that `bias_precision()`, `boot_predict()`, and `jensen_correct()`
  all correctly handle `glmmTMB` models fit with
  `family = lognormal(link = "log")`. This family targets `E[Y]` on the
  response scale via internal MLE handling of residual variance, so only
  random-effect variance requires external correction. An informational
  message is emitted when the correction is applied, clarifying what is
  and is not being corrected.

* Vignette expanded with a note distinguishing `glmmTMB`'s lognormal family
  parameterization from `log(y) ~ ., family = gaussian`, and a new section
  on potential sources of residual bias after Jensen correction.

## Internal

* Improved consistency of informational messages across functions when the
  Jensen correction is applied to lognormal-family models.

GLAMMGoF 1.4.3

## New features

* `boot_predict()` gains a `correction_factor` argument for supplying a
  user-computed Jensen correction on the response scale, overriding the
  internal `jensen_correct(mod)` computation. Useful when only a subset of
  random-effect variance components should contribute to the correction
  (e.g., substantive REs but not nuisance REs), or when the correction has
  been computed externally.

GLAMMGoF 1.4.2

## New features

* boot_predict() now uses joint MVN from the conditional and zero-inflation models when bootstrapping predictions from a zero-inflated model. Previous versions used independent multivariate normal draws for each component. 

GLAMMGoF 1.4.1

## New features

* boot_predict() now returns a boot_median column alongside boot_mean. The mean of exponentiated draws includes a small implicit Jensen adjustment from parameter uncertainty, while boot_median approximates the plug-in point prediction exp(Xβ̂) returned by predict(). Divergence between the two flags high parameter uncertainty. See ?boot_predict for details.

* boot_predict() now recognizes model offsets (formula-side and argument-side) and includes them in the linear predictor. Previously offsets were silently ignored, producing predictions off by a multiplicative factor of exp(offset). A new offset argument lets users specify a reference offset directly. See ?boot_predict for details.


# GLAMMGoF 1.4.0

## New features

* Added a boot_predict() function that allows for bootstrapping of response
  scale predictions (means, standard errors, and confidence intervals at 
  specified alpha labels), with additional options for correction of bias
  from Jensen's inequality due to log-transformed responses or the presence
  of random effects in models with a log-link function. See ?boot_predict for
  details.

# GLAMMGoF 1.3.1

## New features

* Added a simulation figure to the vignette illustrating the effect
  of Jensen's inequality bias on out-of-sample RBIAS across a range
  of random-effect standard deviations, and the effect of applying
  the analytical correction via `bias_adjust = "manual"`. Panels
  cover a GLM without random effects, a GLMM with one random effect,
  and a GLMM with two crossed random effects.

* Added a second simulation figure comparing holdout and case bootstrap
  resampling under a 2x2 design (model specification x observations
  per random-effect group). Demonstrates that case bootstrap
  systematically inflates per-refit random-effect variance when
  groups are sparse, causing `bias_adjust = "manual"` to under-correct;
  the effect is largely independent of model specification.

## Documentation

* Vignette now recommends `method = "holdout"` for random-effect models
  with sparse groups, with the underlying mechanism explained.
  `@param method` in `bias_precision()` and `brier_auc()` updated to
  match.

* Random-slopes note in `bias_precision()` rewritten to reflect that
  `bias_adjust = "manual"` uses a scalar approximation for random-slope
  models; the exact marginal correction under random slopes is
  per-observation.

## Tests

* Regression guards added: `bias_adjust = "tmb"` now errors via
  `match.arg()` in `bias_precision()` and `brier_auc()`.

# GLAMMGoF 1.3.0

## Breaking changes

* `bias_adjust = "tmb"` has been removed from `bias_precision()` and
  `brier_auc()`. Diagnostics on the single-fit path revealed that
  `predict.glmmTMB(..., re.form = ~0, do.bias.correct = TRUE)` does not
  apply a marginal-mean Jensen correction: glmmTMB itself emits the
  warning `'bias.correct' does nothing without random effects` on that
  call, and the returned "bias-corrected" estimate is identical to the
  raw marginal prediction. The option therefore ran expensive TMB
  machinery for no effect. Users should use `bias_adjust = "manual"`
  for log-link GLMMs and lognormal LMMs. For `brier_auc()`, no bias
  correction is currently supported for binomial GLMMs: the logit link
  does not admit a closed-form scalar marginal-mean correction.

## Documentation

* Random slopes note in `bias_precision()` rewritten to reflect that
  `bias_adjust = "manual"` uses a scalar approximation for random-slope
  models. The exact marginal correction under random slopes is
  per-observation (`exp(x_i' Sigma x_i / 2)`) and is planned for a
  future release. The scalar approximation is reasonable when slope
  variance is small or covariates are centered near zero.

* Recommended workflow, examples, and cross-references updated
  throughout to remove the TMB path.

## Tests

* Regression guards added: `bias_adjust = "tmb"` now errors via
  `match.arg()` in both `bias_precision()` and `brier_auc()`.

# GLAMMGoF 1.2.1

## New features

* Bumped to 1.2.1

# GLAMMGoF 1.2.0

## New features

* Bumped to 1.2.0

* `jensen_correct()` now supports log-normal models (a natural-log-transformed
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
