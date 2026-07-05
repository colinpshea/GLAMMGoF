## R CMD check results
0 errors | 1 warning | 0 notes

### Warning
* "package 'ggplot2' was built under R version 4.5.3" — this is a local 
  machine issue and not related to the package itself.

## Test environments
* Local Windows 11, R 4.5.x
* win-builder (devel)
* win-builder (release)

## Downstream dependencies
This is a new submission. There are no downstream dependencies.

## Notes on package data
Example datasets (`countData`, `logitData`) and fitted model objects 
(`countModel_GLM`, `countModel_GLMM`, etc.) are included in the `data/` 
directory to allow users to immediately explore package functionality 
without needing to fit their own models first. Since the core functions 
require an existing fitted model object as input, these examples provide 
a convenient starting point for new users.

## Resampling methods
Both functions support two resampling methods via the `method` argument: 
repeated random holdout (Monte Carlo cross-validation, `method = "holdout"`, 
the default) and bootstrap resampling with out-of-bag evaluation 
(`method = "bootstrap"`).

## Possibly misspelled words
The following flagged words are intentional technical terminology and not 
misspellings:
* RRMSE (Relative Root Mean Square Error)
* RMAE (Relative Mean Absolute Error)
* RMedAE (Relative Median Absolute Error)
* RBIAS (Relative Bias)
* AUC (Area Under Curve statistic, a standard predictive performance metric)
* Brier (Brier score, a standard predictive performance metric)
* log loss (a standard predictive performance metric)
* DHARMa (R package name)
* glmmTMB (R package name)
* lme (as in lme4, an R package name)
* mgcv (R package name)

## Notes to CRAN reviewers

This is the first submission of GLAMMGoF to CRAN. The package is currently
distributed via R-universe (https://cshea15.r-universe.dev/GLAMMGoF) and has
been in active use by ecological researchers prior to this submission.

The package provides resampling-based predictive validation for generalized
linear and additive models, with particular support for detecting and correcting
Jensen's inequality bias in log-link GLMMs with random effects. Supported model
backends include glmmTMB, lme4, mgcv, MASS, and stats.

All examples in the documentation are wrapped in \dontrun{} to avoid long
runtimes during checking, as the core functions involve repeated model refitting
via Monte Carlo cross-validation or bootstrap resampling. A vignette
demonstrating full package functionality is included.

The package imports DHARMa for optional residual diagnostics. DHARMa is listed
in Imports rather than Suggests because its simulateResiduals() function is
called within the main exported functions when DHARMaPlot = TRUE (the default).
