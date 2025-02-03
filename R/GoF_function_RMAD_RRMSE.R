#' Bootstrap RRMSE and RMAD fit statistics
#'
#' @description Bootstrap RRMSE and RMAD fit statistics for generalized linear models with continuous or integer response variables and with or without random effects. The fit statistics are relative root mean squared error (RRMSE), calculated as sqrt(mean((observed - predicted)^2))/mean(observed)*100, and relative median absolute deviation (RMAD), calculated as median(abs((observed - predicted)))/mean(observed)*100. The RMAD is generally less sensitive than RRMSE to extreme values.
#' @param nReps Desired number of bootstrap replicates. The default value is 100.
#' @param testModel A regression model fit to testData in `glmmTMB` (with or without random effects), `glmer` (with random effects), or `glm`/`lm` (without random effects). The response variable can be continuous or an integer, and possible error distributions include Poisson, negative binomial, gamma, tweedie, and gaussian.
#' @param testData A data frame with a continuous or integer response variable and continuous and/or categorical predictors.
#' @param propTrain Proportion of testData that is used for model-fitting and in-sample predictive performance (the remaining % is used to assess out-of-sample predictive performance). The default value is 0.8.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the `simulateResiduals()` function of the `DHARMa` package? The default is "Yes".
#' @return This function returns four objects: a data frame with all of the bootstrapping results (i.e., all nReps bootstrapped values for each performance statistic), a data frame with a summary (mean and 95% CLs) of all bootstrap replicates for each performance statistic, a histogram of values for each performance statistic, and a goodness-of-fit plot based on scaled residuals from the `simulateResiduals()` function of the `DHARMa` package.
#'
#' This package contains an example data set for a negative binomial or Poisson regression called countData (but data with a continuous response variable could also be used). Two example negative binomial regression model objects are also included called countModel1, which includes a random effect, and countModel2, which does not; both models were fitted using glmmTMB, but countModel1 could also be a `glmer` model object (fitted using `glmer` or `glmer.nb` from `lme4`) and countModel2 could also be a `glm.nb` (from the `MASS` package) model object:
#'
#' countModel1 <- glmmTMB(y ~ Season + River + Temp + Snags + Year + AvgDepth + (1|RiverSeasonYear), data = countDat, family = nbinom2)
#'
#' countModel2 <- glmmTMB(y ~ Season + River + Temp + Snags + Year + AvgDepth, data = countDat, family = nbinom2)
#'
#' Bootstrapping the fit statistics requires specifying the data and model being tested, the desired number of bootstrap replicates, and the proportion of data used in the training (in-sample performance) data set:
#'
#' RRMSE_RMAD(nReps = 100, testModel = NULL, testData = countData, propTrain = 0.8)
#' @importFrom magrittr %>%
#' @importFrom dplyr select group_by summarize summarise mutate bind_rows n
#' @importFrom tidyr pivot_longer pivot_wider separate
#' @importFrom DHARMa simulateResiduals
#' @importFrom glmmTMB ranef glmmTMB
#' @importFrom lme4 glmer lmer glmer.nb
#' @importFrom MASS glm.nb
#' @export
RRMSE_RMAD <- function(nReps = 100, testModel = NULL, testData = NULL, propTrain = 0.8, DHARMaPlot = "Yes"){
  fit_cost_rrmse <- function(y, yhat){sqrt(mean((y - yhat)^2))/mean(y)*100}
  fit_cost_rmad <- function(y, yhat){median(abs((y - yhat)))/mean(y)*100}
  cost_test_fin_RRMSE = NULL
  cost_train_fin_RRMSE = NULL
  cost_test_fin_RMAD = NULL
  cost_train_fin_RMAD = NULL
  testResp <- function(data){length(unique(data))==2 && all(data %in% c(0, 1))}
  stopifnot("Response variable is binary! Use BRIER_AUC() instead" = testResp(unname(unlist(eval(as.symbol(paste0("testData")))[,all.vars(formula(testModel))[1]])))=="FALSE")
  for (j in 1:nReps){
    smp_size <- floor(propTrain*nrow(testData))
    train_ind <- sample(seq_len(nrow(testData)), size = smp_size)
    train <- testData[train_ind, ]
    test <-  testData[-train_ind, ]
    if ("glmmTMB" %in% class(testModel)) {m_train <- glmmTMB(formula(testModel), family = family(testModel), data = train)}
    if ("glmerMod" %in% class(testModel)) {
      if ((grepl("Negative Binomial", family(testModel)$family))) {
        try(m_train <- glmer.nb(formula(testModel), data = train))
      }
    }
    if ("glmerMod" %in% class(testModel)) {
      if (!(grepl("Negative Binomial", family(testModel)$family))) {
        try(m_train <- glmer(formula(testModel), family = family(testModel), data = train))
      }
    }
    if ("lmerMod" %in% class(testModel)) {try(m_train <- lmer(formula(testModel), data = train))}
    if ("negbin" %in% class(testModel)) {m_train <- glm.nb(formula(testModel), data = train)}
    if ("glm" %in% class(testModel)) {m_train <- glm(formula(testModel), family = family(testModel), data = train)}
    if ("lm" %in% class(testModel)) {m_train <- lm(formula(testModel), data = train)}
    if ("glmmTMB" %in% class(testModel)) {
      if (sum(ranef(testModel)=="list()")<length(ranef(testModel))){
      train_pred <- train
      train_pred[,which(names(train_pred) %in% names(ranef(testModel)$cond))] <- NA
      test_pred <- test
      test_pred[,which(names(test_pred) %in% names(ranef(testModel)$cond))] <- NA
      cost_train_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("train_pred")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = train_pred))
      cost_test_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("test_pred")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = test_pred))
      cost_train_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("train_pred")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = train_pred))
      cost_test_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("test_pred")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = test_pred))
      }
       if (sum(ranef(testModel)=="list()")==length(ranef(testModel))){
          cost_train_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = train))
          cost_test_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = test))
          cost_train_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = train))
          cost_test_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = test))
        }
      }
      if (any(c("negbin", "lm", "glm") %in% class(testModel))){
        cost_train_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = train))
        cost_test_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = test))
        cost_train_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = train))
        cost_test_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", newdata = test))
      }
    if (any(c("glmerMod", "lmerMod") %in% class(testModel))){
      cost_train_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", re.form = ~0, newdata = train))
      cost_test_fin_RRMSE[j] <- fit_cost_rrmse(y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", re.form = ~0, newdata = test))
      cost_train_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", re.form = ~0, newdata = train))
      cost_test_fin_RMAD[j] <- fit_cost_rmad(y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), yhat = predict(m_train, type = "response", re.form = ~0, newdata = test))
    }
  }
  results_list <- list(train_RRMSE = cost_train_fin_RRMSE, test_RRMSE = cost_test_fin_RRMSE, train_RMAD = cost_train_fin_RMAD, test_RMAD = cost_test_fin_RMAD)

  results_df <- bind_rows(results_list, .id = "column_label") %>% mutate(simRep = 1:n()) %>% pivot_longer(cols = -simRep, values_to = "value", names_to = "metric") %>% separate(metric, into = c("Group", "Metric")) %>% mutate(Group = factor(Group, levels = c("train", "test"), labels = c("In-sample performance", "Out-of-sample performance")), Metric = factor(Metric, levels = c("RRMSE", "RMAD"), labels = c("RRMSE", "RMAD")))

  results_summary <- results_df %>% group_by(Group, Metric) %>% summarise(mn = mean(value), lwr95 = quantile(value, 0.025), upr95 = quantile(value, 0.975))

  if (DHARMaPlot=="Yes"){
    dharmaPlot <- simulateResiduals(n = 1000, testModel, plot = T)
  results_plot <- ggplot(results_df, aes(x = value)) + geom_histogram(color = "black", fill = "grey") + facet_grid(Group~Metric, scales = "free") + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey90", linetype = "solid"), panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), axis.text = element_text(colour = "black")) + labs(x = "% relative to true mean", y = "Frequency") + theme(panel.spacing = unit(1.5, "lines"))
  return(list(rrmse_rmad_results = results_df, rrmse_rmad_hist = results_plot, rrmse_rmad_summary = results_summary, dharmaPlot = dharmaPlot))
  }
  if (DHARMaPlot=="No"){
  return(list(rrmse_rmad_results = results_df, rrmse_rmad_hist = results_plot, rrmse_rmad_summary = results_summary))
  }
}
