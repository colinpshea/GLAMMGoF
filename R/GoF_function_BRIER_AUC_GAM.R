#' Bootstrap Brier score and AUC fit statistics
#'
#' @description Bootstrap Brier score and AUC fit statistics (see `rms` package documentation for details) for generalized linear and generalized additive models with binary response variables and with or without random effects. Brier scores range from 0 to 1, with values closer to 0 indicating a better-predicting model, and where sqrt(Brier score) is the average difference, across all observations, between the predicted probability and the observed value (0 or 1). Conversely, AUC statistics range from 0 to 1, where values closer to 1 indicate a better-predicting model.
#' @param nReps Desired number of bootstrap replicates. The default value is 100, but this number should be at least 1000 in practice.
#' @param testModel A logistic regression model fitted to testData using `glmmTMB` (with or without random effects), `glmer` (with random effects), or `glm` (without random effects).
#' @param testData A data frame with a binary response variable and continuous and/or categorical predictors variables.
#' @param propTrain Proportion of `testData` that is used for model-fitting and in-sample predictive performance (the default value is 0.8). The remaining % is used to assess out-of-sample predictive performance.
#' @param DHARMaPlot Do you want to return a goodness-of-fit plot from the `simulateResiduals()` function of the `DHARMa` package? The default is `TRUE`. You can also specify `DHARMaReps` if you want something other than the default of 1000 simulation replicates.
#' @return This function returns four objects: a data frame with all of the bootstrapping results (i.e., all `nReps` bootstrapped values for each performance statistic), a data frame with a summary (mean and 95% CLs) of all bootstrap replicates for each performance statistic, a histogram of values for each performance statistic, and a goodness-of-fit plot based on scaled residuals from the `simulateResiduals()` function of the `DHARMa` package. If DHARMaPlot = `FALSE`, then `simulateResiduals()` isn't used to assess the model's residuals and only three of the four objects are returned.
#'
#' This package contains an example data set to fit a logistic regression called logitData. Two example logistic regression model objects are also included: logitModel1 includes a random effect, and logitModel2 does not; both models were fitted in `glmmTMB`, but logitModel1 could also be a `glmer` (from `lme4`) or gam (from `mgcv`) model object and logitModel2 could also be `glm` or `gam` (from `mgcv`) model object:
#'
#' logitModel1 <- glmmTMB(y ~ totalLengthcm + Zone + (1|Year), family = binomial, data = logitData)
#'
#' logitModel2 <- glmmTMB(y ~ totalLengthcm + Zone, family = binomial, data = logitData)
#'
#' Bootstrapping the performance statistics requires specifying the data and model being tested, the desired number of bootstrap replicates (the default is 100 but it should be higher in practice), the proportion of data used in the training (in-sample performance) data set, whether you want use DHARMa to assess the residuals (the default is TRUE), and how many simulation replicates you want to use in DHARMa's simulateResiduals() function (the default is 1000):
#'
#' BRIER_AUC(nReps = 100, testModel = logitModel, testData = logitData, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000)
#' @importFrom magrittr %>%
#' @importFrom dplyr select group_by summarize summarise mutate bind_rows n
#' @importFrom tidyr pivot_longer pivot_wider separate
#' @importFrom DHARMa simulateResiduals
#' @importFrom rms val.prob
#' @importFrom glmmTMB ranef glmmTMB
#' @importFrom lme4 glmer lmer glmer.nb
#' @importFrom mgcv gam predict.gam
#' @export
BRIER_AUC_GAM <- function(nReps = 100, testModel = NULL, testData = NULL, propTrain = 0.8, DHARMaPlot = TRUE, DHARMaReps = 1000){
  auc_train = NULL
  brier_train = NULL
  auc_test = NULL
  brier_test = NULL
  testResp <- function(data) {
    length(unique(data)) == 2 && all(data %in% c(0, 1))
  }
  stopifnot(`Response variable is not binary! Use RRMSE_RMAD() instead` = testResp(unname(unlist(eval(as.symbol(paste0("testData")))[,
                                                                                                                                     all.vars(formula(testModel))[1]]))) == "TRUE")
  for (j in 1:nReps) {
    smp_size <- floor(propTrain * nrow(testData))
    train_ind <- sample(seq_len(nrow(testData)), size = smp_size)
    train <- testData[train_ind, ]
    test <- testData[-train_ind, ]
    if ("glmmTMB" %in% class(testModel)) {
      m_train <- glmmTMB(formula(testModel), family = family(testModel),
                         data = train)
    }
    if ("gam" %in% class(testModel)){
      m_train <- gam(formula(testModel), family = family(testModel), data = train)
    }
    if ("glmerMod" %in% class(testModel)) {
      try(m_train <- glmer(formula(testModel), family = family(testModel),
                           data = train))
    }
    if ("glm" %in% class(testModel) & !("gam" %in% class(testModel))) {
      m_train <- glm(formula(testModel), family = family(testModel),
                     data = train)
    }
    if ("glmmTMB" %in% class(testModel)) {
      if (sum(ranef(testModel) == "list()") < length(ranef(testModel))) {
        train_pred <- train
        train_pred[, which(names(train_pred) %in% names(ranef(testModel)$cond))] <- NA
        test_pred <- test
        test_pred[, which(names(test_pred) %in% names(ranef(testModel)$cond))] <- NA
        auc_train[j] <- val.prob(p = predict(m_train,
                                             type = "response", newdata = train_pred), y = unname(unlist(eval(as.symbol(paste0("train_pred")))[,
                                                                                                                                               all.vars(formula(testModel))[1]])), smooth = FALSE,
                                 pl = FALSE)["C (ROC)"]
        brier_train[j] <- val.prob(p = predict(m_train,
                                               type = "response", newdata = train_pred), y = unname(unlist(eval(as.symbol(paste0("train_pred")))[,
                                                                                                                                                 all.vars(formula(testModel))[1]])), smooth = FALSE,
                                   pl = FALSE)["Brier"]
        auc_test[j] <- val.prob(p = predict(m_train,
                                            type = "response", newdata = test_pred), y = unname(unlist(eval(as.symbol(paste0("test_pred")))[,
                                                                                                                                            all.vars(formula(testModel))[1]])), smooth = FALSE,
                                pl = FALSE)["C (ROC)"]
        brier_test[j] <- val.prob(p = predict(m_train,
                                              type = "response", newdata = test_pred), y = unname(unlist(eval(as.symbol(paste0("test_pred")))[,
                                                                                                                                              all.vars(formula(testModel))[1]])), smooth = FALSE,
                                  pl = FALSE)["Brier"]
      }
      if (sum(ranef(testModel) == "list()") == length(ranef(testModel))) {
        auc_train[j] <- val.prob(p = predict(m_train,
                                             type = "response", newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,
                                                                                                                                     all.vars(formula(testModel))[1]])), smooth = FALSE,
                                 pl = FALSE)["C (ROC)"]
        brier_train[j] <- val.prob(p = predict(m_train,
                                               type = "response", newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,
                                                                                                                                       all.vars(formula(testModel))[1]])), smooth = FALSE,
                                   pl = FALSE)["Brier"]
        auc_test[j] <- val.prob(p = predict(m_train,
                                            type = "response", newdata = test), y = unname(unlist(eval(as.symbol(paste0("test")))[,
                                                                                                                                  all.vars(formula(testModel))[1]])), smooth = FALSE,
                                pl = FALSE)["C (ROC)"]
        brier_test[j] <- val.prob(p = predict(m_train,
                                              type = "response", newdata = test), y = unname(unlist(eval(as.symbol(paste0("test")))[,
                                                                                                                                    all.vars(formula(testModel))[1]])), smooth = FALSE,
                                  pl = FALSE)["Brier"]
      }
    }
    if ("glmerMod" %in% class(testModel)) {
      auc_train[j] <- val.prob(p = predict(m_train, type = "response",
                                           newdata = train, re.form = NA), y = unname(unlist(eval(as.symbol(paste0("train")))[,
                                                                                                                              all.vars(formula(testModel))[1]])), smooth = FALSE,
                               pl = FALSE)["C (ROC)"]
      brier_train[j] <- val.prob(p = predict(m_train, type = "response",
                                             newdata = train, re.form = NA), y = unname(unlist(eval(as.symbol(paste0("train")))[,
                                                                                                                                all.vars(formula(testModel))[1]])), smooth = FALSE,
                                 pl = FALSE)["Brier"]
      auc_test[j] <- val.prob(p = predict(m_train, type = "response",
                                          newdata = test, re.form = NA), y = unname(unlist(eval(as.symbol(paste0("test")))[,
                                                                                                                           all.vars(formula(testModel))[1]])), smooth = FALSE,
                              pl = FALSE)["C (ROC)"]
      brier_test[j] <- val.prob(p = predict(m_train, type = "response",
                                            newdata = test, re.form = NA), y = unname(unlist(eval(as.symbol(paste0("test")))[,
                                                                                                                             all.vars(formula(testModel))[1]])), smooth = FALSE,
                                pl = FALSE)["Brier"]
    }

    if (!("gam" %in% class(testModel)) & any(c("glm", "lm") %in% class(testModel))){
      auc_train[j] <-  val.prob(p = predict(m_train, type="response", newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["C (ROC)"]
      brier_train[j] <- val.prob(p = predict(m_train, type="response", newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["Brier"]
      auc_test[j] <- val.prob(p = predict(m_train, type="response", newdata = test), y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["C (ROC)"]
      brier_test[j] <- val.prob(p = predict(m_train, type="response", newdata = test), y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["Brier"]
    }
    if ("gam" %in% class(testModel)){
      if (length(testModel$smooth[lengths(lapply(testModel$smooth, function(x) x$random==TRUE))>0]) > 0){
        re_name <- testModel$smooth[lengths(lapply(testModel$smooth, function(x) x$random==TRUE))>0][[1]]$label
        train[which(colnames(train) %in% str_match(re_name, "\\((.*)\\)"))] <- 0
        test[which(colnames(test) %in% str_match(re_name, "\\((.*)\\)"))] <- 0
        auc_train[j] <-  val.prob(p = predict.gam(m_train, type="response", newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["C (ROC)"]
        brier_train[j] <- val.prob(p = predict.gam(m_train, type="response",newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["Brier"]
        auc_test[j] <- val.prob(p = predict.gam(m_train, type="response", newdata = test, newdata.guaranteed = TRUE), y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["C (ROC)"]
        brier_test[j] <- val.prob(p = predict.gam(m_train, type="response", newdata = test, newdata.guaranteed = TRUE), y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["Brier"]
      }
    }
    if ("gam" %in% class(testModel)){
      if (length(testModel$smooth[lengths(lapply(testModel$smooth, function(x) x$random==TRUE))>0]) == 0){
        auc_train[j] <- val.prob(p = predict(m_train, type="response", newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["C (ROC)"]
        brier_train[j] <- val.prob(p = predict(m_train, type="response", newdata = train), y = unname(unlist(eval(as.symbol(paste0("train")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["Brier"]
        auc_test[j] <- val.prob(p = predict(m_train, type="response", newdata = test), y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["C (ROC)"]
        brier_test[j] <- val.prob(p = predict(m_train, type="response", newdata = test), y = unname(unlist(eval(as.symbol(paste0("test")))[,all.vars(formula(testModel))[1]])), smooth = FALSE, pl = FALSE)["Brier"]
      }
    }
  }
  results_list <- list(auc_train = unname(auc_train), brier_train = unname(brier_train),
                       auc_test = unname(auc_test), brier_test = unname(brier_test))
  results_df <- bind_rows(results_list, .id = "column_label") %>%
    mutate(simRep = 1:n()) %>% pivot_longer(cols = -simRep,
                                            values_to = "value", names_to = "metric") %>% separate(metric,
                                                                                                   into = c("Metric", "Group")) %>% mutate(Group = factor(Group,
                                                                                                                                                          levels = c("train", "test"), labels = c("In-sample performance",
                                                                                                                                                                                                  "Out-of-sample performance")), Metric = factor(Metric,
                                                                                                                                                                                                                                                 levels = c("auc", "brier"), labels = c("AUC statistic",
                                                                                                                                                                                                                                                                                        "Brier score")))
  results_summary <- results_df %>% group_by(Group, Metric) %>%
    summarise(mn = mean(value), lwr95 = quantile(value, 0.025),
              upr95 = quantile(value, 0.975))
  results_plot <- ggplot(results_df, aes(x = value)) + geom_histogram(color = "black",
                                                                      fill = "grey") + facet_grid(Group ~ Metric) + theme_bw() +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(colour = "grey90",
                                                                                  linetype = "solid"), panel.grid.minor.y = element_line(colour = "grey90",
                                                                                                                                         linetype = "dashed"), axis.text = element_text(colour = "black")) +
    labs(x = "% relative to true mean", y = "Frequency") +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0,
                                                      1, 0.2), expand = expansion(add = c(0.05, 0.05))) +
    theme(panel.spacing = unit(1.5, "lines"))
  library(DHARMa)
  if (DHARMaPlot == TRUE) {
    dharmaPlot <- simulateResiduals(n = DHARMaReps, testModel,
                                    plot = T)
        return(list(brier_auc_results = results_df, brier_auc_hist = results_plot,
                    brier_auc_summary = results_summary, dharmaPlot = dharmaPlot))
  }
  if (DHARMaPlot == FALSE) {
       return(list(brier_auc_results = results_df, brier_auc_hist = results_plot,
                    brier_auc_summary = results_summary))
  }
  }
