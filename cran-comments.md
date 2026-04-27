## R CMD check results
0 errors | 1 warning | 0 notes

### Warning
* "package 'ggplot2' was built under R version 4.5.3" — this is a local 
  machine issue and not related to the package itself.

### Note
* "unable to verify current time" — this is a network issue on the check 
  machine and not related to the package itself.

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
