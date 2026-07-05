# GLAMMGoF 1.1.2

* Updated `@param`, `@note`, and vignette documentation to reflect `lme4` support for `bias_adjust = "manual"` and `conditional_predictions`

# GLAMMGoF 1.1.1

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
