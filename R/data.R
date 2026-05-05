#' Simulated count data
#'
#' A simulated dataset with an integer count response and predictor variables
#' for demonstrating the bias_precision function.
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{y}{Integer count response variable}
#'   \item{Season}{Factor with 4 levels: Autumn, Spring, Summer, Winter}
#'   \item{Temp}{Continuous temperature predictor}
#'   \item{Site}{Factor with 10 levels representing sampling sites}
#'   \item{Year}{Factor with 8 levels representing sampling years (2015-2022)}
#' }
"countData"

#' Simulated binary data
#'
#' A simulated dataset with a binary response and predictor variables
#' for demonstrating the brier_auc function.
#'
#' @format A data frame with 1000 rows and 5 variables:
#' \describe{
#'   \item{y}{Binary response variable (0 or 1)}
#'   \item{Season}{Factor with 4 levels: Autumn, Spring, Summer, Winter}
#'   \item{Temp}{Continuous temperature predictor}
#'   \item{Site}{Factor with 10 levels representing sampling sites}
#'   \item{Year}{Factor with 8 levels representing sampling years (2015-2022)}
#' }
"logitData"

#' Simulated count GLM example model
#'
#' A glmmTMB model fit to countData with fixed effects only.
#' Formula: y ~ Season + Temp
#' @format A glmmTMB model object
"countModel_GLM"

#' Simulated count GLMM example model
#'
#' A glmmTMB model fit to countData with one random effect.
#' Formula: y ~ Season + Temp + (1 | Site)
#' @format A glmmTMB model object
"countModel_GLMM"

#' Simulated count GLMM2 example model
#'
#' A glmmTMB model fit to countData with two random effects.
#' Formula: y ~ Season + Temp + (1 | Site) + (1 | Year)
#' @format A glmmTMB model object
"countModel_GLMM2"

#' Simulated count GAM example model
#'
#' A mgcv GAM fit to countData with fixed effects only.
#' Formula: y ~ Season + s(Temp)
#' @format A gam model object
"countModel_GAM"

#' Simulated count GAMM example model
#'
#' A mgcv GAM fit to countData with one random effect.
#' Formula: y ~ Season + s(Temp) + s(Site, bs = "re")
#' @format A gam model object
"countModel_GAMM"

#' Simulated count GAMM2 example model
#'
#' A mgcv GAM fit to countData with two random effects.
#' Formula: y ~ Season + s(Temp) + s(Site, bs = "re") + s(Year, bs = "re")
#' @format A gam model object
"countModel_GAMM2"

#' Simulated binary GLM example model
#'
#' A glmmTMB model fit to logitData with fixed effects only.
#' Formula: y ~ Season + Temp
#' @format A glmmTMB model object
"logitModel_GLM"

#' Simulated binary GLMM example model
#'
#' A glmmTMB model fit to logitData with one random effect.
#' Formula: y ~ Season + Temp + (1 | Site)
#' @format A glmmTMB model object
"logitModel_GLMM"

#' Simulated binary GLMM2 example model
#'
#' A glmmTMB model fit to logitData with two random effects.
#' Formula: y ~ Season + Temp + (1 | Site) + (1 | Year)
#' @format A glmmTMB model object
"logitModel_GLMM2"

#' Simulated binary GAM example model
#'
#' A mgcv GAM fit to logitData with fixed effects only.
#' Formula: y ~ Season + s(Temp)
#' @format A gam model object
"logitModel_GAM"

#' Simulated binary GAMM example model
#'
#' A mgcv GAM fit to logitData with one random effect.
#' Formula: y ~ Season + s(Temp) + s(Site, bs = "re")
#' @format A gam model object
"logitModel_GAMM"

#' Simulated binary GAMM2 example model
#'
#' A mgcv GAM fit to logitData with two random effects.
#' Formula: y ~ Season + s(Temp) + s(Site, bs = "re") + s(Year, bs = "re")
#' @format A gam model object
"logitModel_GAMM2"
