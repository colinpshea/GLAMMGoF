# Pre-computes bias_precision() summaries for the countModel_GLMM
# worked example in the tutorial vignette. Cached to inst/extdata/ so
# the vignette renders quickly on R-universe and CRAN.

# Rerun this script and rebuild the vignette after any change to
# countData or countModel_GLMM in R/data-*.R.

library(GLAMMGoF)

set.seed(123)

bp_countModelGLMM_uncorrected <- bias_precision(
  nReps       = 1000,
  testModel   = countModel_GLMM,
  testData    = countData,
  propTrain   = 0.8,
  DHARMaPlot  = FALSE,
  seed        = 123,
  method      = "holdout",
  bias_adjust = "none"
)

bp_countModelGLMM_corrected <- bias_precision(
  nReps       = 1000,
  testModel   = countModel_GLMM,
  testData    = countData,
  propTrain   = 0.8,
  DHARMaPlot  = FALSE,
  seed        = 123,
  method      = "holdout",
  bias_adjust = "manual"
)

saveRDS(bp_countModelGLMM_uncorrected,
        file = "C:/Users/Colin.Shea/OneDrive - Florida Fish and Wildlife Conservation/R Packages/GLAMMGoF/inst/extdata/bp_countModelGLMM_uncorrected.rds")
saveRDS(bp_countModelGLMM_corrected,
        file = "C:/Users/Colin.Shea/OneDrive - Florida Fish and Wildlife Conservation/R Packages/GLAMMGoF/inst/extdata/bp_countModelGLMM_corrected.rds")


