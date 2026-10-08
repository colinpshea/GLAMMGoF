# AI-Assisted Development Disclosure

## Overview

The `GLAMMGoF` R package was developed iteratively with assistance from a generative AI tool (Claude, Anthropic; most recently Claude Opus 4.7, with earlier Claude versions used during earlier stages of development). This document is intended to provide transparency about the nature and extent of AI involvement in the package's development in a manner consistent with the AI disclosure statement in the accompanying manuscript (Shea, XXXX). 

## Nature of AI assistance

AI assistance was used during package development for the following activities: 

- **Consultation on R programming approaches** - discussion of implementation options and code design for package functions
- **Code review and suggestions for reviewing code efficiency** - identifying redundant code, suggesting vectorized alternatives, and recommending refactoring to improve clarity, efficiency, and generalization of the package's functions
- **Debugging assistance** - diagnosing error messages, identifying root causes of unexpected or incorrect behavior, and suggesting fixes
- **Development of tests for the package's functions** - suggesting the design of test structures, including edge cases, and reviewing the comprehensiveness of overall test coverage for the package's functions
- **Documentation drafting and review** - refining roxygen documentation to improve clarity of function descriptions and examples

## Nature of author-driven work

The following aspects of this R package and the accompanying paper were author-driven: 

- **Study design and methodological choices** - the decision to pursue this project; the conception and initial implementation of `GLAMMGoF` and its core functions (`bias_precision()`, `brier_auc()`, and `boot_predict()`); the scope of how Jensen bias applies to statistical ecology practice; and the choice of model classes and R packages to support
- **Mathematical derivations** - the Jensen correction formula; the Gaussian moment-generating function application; the initial derivation of the row-wise correction for `dispformula` covariates; and the initial design of the package's bootstrapping function used for prediction
- **Selection and validation of the package's core algorithms** - the choice of making `holdout` and `bootstrap` resampling schemes available in `bias_precision()` and `brier_auc()`; the choice and implementation of performance measures for the package's core functions; the parametric bootstrap approach in `boot_predict()`; and the handling of variance components and conditional-at-zero, conditional, and population-averaged predictions across model classes
- **Interpretation of results and ecological context** - the framing of the practice gap; the arguments about when the correction does and does not apply; the limitations and additional considerations when applying the Jensen correction; and the recommended workflows detailed in the accompanying paper and package vignette
- **Scope decisions** - the choice to focus on log-link GLMMs and log-transformed LMMs with independent Gaussian random intercepts; and the decision not to extend to models with more complex random effects structures such as those with spatially or temporally correlated random effects or random slopes

## Verification

All AI-assisted code was reviewed, tested, and verified by the author. The package includes more than 200 automated tests that are designed to verify function behavior across supported model classes and use cases. The package's behavior on standard datasets has been compared with similar implementations in the R packages `emmeans`, `ggeffects`, `modelbased`, and `marginaleffects` (see companion manuscript for details).

## Function-level annotation

Due to the iterative nature of AI assistance during the development of `GLAMMGoF`, AI-assisted contributions cannot be cleanly separated from author-written code within specific functions. AI assistance was distributed across the package's development rather than being confined to identifiable code blocks: a given function may have had its logic authored by the user, its structure refined through AI-assisted refactoring, its documentation refined with AI assistance, and its tests designed with AI help. Function-level annotation would therefore have required flagging essentially the entire package (thereby losing the informative value of the annotation) or arbitrarily deciding which functions to flag as AI-assisted. This package-level disclosure accurately reflects the development of `GLAMMGoF`. 

## Responsibility

The author takes full responsibility for the package's code and functionality, including any errors. Users encountering issues or suggesting improvements should file issues at the package's repository. 

## Package information

- **Package:** GLAMMGoF
- **Concept DOI:** https://doi.org/10.5281/zenodo.22666285 (resolves to the latest version)
- **Author:** Colin P. Shea, Florida Fish and Wildlife Conservation Commission / Fish and Wildlife Research Institute
- **Contact:** colin.shea@myfwc.com

## Citation

If you use GLAMMGoF in published work, please cite:

**Shea, C.P. (XXXX). [Manuscript title]. [Journal], [volume/pages/DOI].**

The package can be cited using its Zenodo concept DOI:

**Shea, C.P. (XXXX). GLAMMGoF: Resampling-based predictive validation for generalized linear and generalized additive models. R package version X.X.X. Zenodo. https://doi.org/10.5281/zenodo.22666285**

## Document version

- **Version:** 2.0
- **Last updated:** October 8, 2026
