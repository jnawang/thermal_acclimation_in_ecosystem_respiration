# This script identifies drivers of apparent and actual thermal response strength
# Authors: Junna Wang, October, 2025

# It takes 10min to run this script. 

library(librarian)
shelf(dplyr, ggplot2, corrplot, MuMIn, fBasics, car, lme4, vip, randomForest, caret, VSURF, rsample)
rm(list=ls())

#
acclimation <- read.csv('data/acclimation_data_TA.csv')

acclimation <- acclimation %>% filter(!site_ID %in% c('IT-Lav', "IT-TrF")) # The two sites night NEE has problems. 

#------------------------do some explorations first-----------------------------
summary(lm(data=acclimation, TAS ~ IGBP))
summary(lm(data=acclimation, TAS ~ Climate_class))
summary(aov(TAS ~ IGBP, data=acclimation))
summary(aov(TAS ~ Climate_class, data=acclimation))

# total TRS--strongly affected by both climate (p = 0.002592) and vegetation (p = 5.77e-05)
summary(lm(data=acclimation, TAS_tot ~ IGBP))
summary(lm(data=acclimation, TAS_tot ~ Climate_class))

# Do I want to combine climate and vegetation groups?
# combine climate class
table(acclimation$Climate_class)
# tropical (Am, n=1)
# semi-arid (Bs, n = 8)
# arid (Bw, n = 1)
# hot summer Mediterranean (Csa, n = 10)
# warm summer Mediterranean (Csb, n = 4)
# humid subtropical (Cfa, Cfb, Cfc, n = 27)
# humid continental (Dfa, Dfb, n = 37)
# subarctic (Dfc, Dfd, Dwc, n = 19)
# tundra (ET, n = 3)

# combine vegetation class
table(acclimation$IGBP)
# combine DBF and DNF: these are decidulous forests
# combine SAV and WSA: these are savannah

acclimation <- acclimation %>% mutate(Climate_class_new = case_when(Climate_class %in% c("Bsh", "Bsk") ~ 'Bs', 
                                                                    Climate_class %in% c("Bwk") ~ 'Bw',
                                                                    Climate_class %in% c("Cfa", "Cfb", "Cfc") ~ 'Cf',
                                                                    Climate_class %in% c("Csa") ~ 'Csa',
                                                                    Climate_class %in% c("Csb") ~ 'Csb',
                                                                    Climate_class %in% c("Dfa", "Dfb") ~ 'Df',
                                                                    Climate_class %in% c("Dfc, Dfd", "Dwc") ~ 'Subartic',
                                                                    Climate_class %in% c("ET") ~ 'ET')) %>%
  mutate(IGBP_new = case_when(IGBP %in% c("CSH") ~ 'CSH', IGBP %in% c("DBF") ~ 'DBF', IGBP %in% c("EBF") ~ 'EBF',
                              IGBP %in% c("ENF") ~ 'ENF', IGBP %in% c("GRA") ~ 'GRA', IGBP %in% c("MF") ~ 'MF',
                              IGBP %in% c("OSH") ~ 'OSH', IGBP %in% c("SAV", "WSA") ~ 'SAV', IGBP %in% c("WET") ~ 'WET'))

summary(lm(data=acclimation, TAS_tot ~ Climate_class_new)) # p-value: 1.795e-06
summary(lm(data=acclimation, TAS_tot ~ IGBP_new))          # p-value: 0.001007

summary(lm(data=acclimation, TAS ~ Climate_class_new)) # p-value: 0.05199
summary(lm(data=acclimation, TAS ~ IGBP_new))          # p-value: 0.0425

# Yes, the significance levels are almost the same with ungrouped cases. 


#------------------------correlation among predictor variables------------------
data.cor <- acclimation[, c(4, 7:20, 22)]
cor(data.cor)

#---------------contribution of direct and apparent thermal responses to TAS_tot----------
acclimation$TAS_app <- acclimation$TAS_tot - acclimation$TAS
(var(acclimation$TAS_app) + cov(acclimation$TAS_app, acclimation$TAS)) / var(acclimation$TAS_tot)  
(var(acclimation$TAS) + cov(acclimation$TAS_app, acclimation$TAS)) / var(acclimation$TAS_tot)

#----------------------------try simple linear regression-----------------------
summary(lm(data=acclimation, TAS ~ MATA + ELEV + LAI + SOC))          # Only productivity is most important. 
summary(lm(data=acclimation, TAS_tot ~ MATA + ELEV + LAI + SOC))      # The first 3 are all important. 

#-------------------------------------------------------------------------------
#----------------------try random forests, used in manuscript-------------------
set.seed(222)
data_TAS <- acclimation[, c("TAS", "ELEV", "MATA", "LAI", "SOC")]

# calculate vif of the four predictors
lm_model <- lm(data=data_TAS, TAS ~ .)
car::vif(lm_model)
# vif: variance inflation factors
#    ELEV     MATA     LAI     SOC 
# 1.283012 1.163074 1.241224 1.149224
#
# stratify the data, and create folds
acclimation$Climate_reclass <- acclimation$Climate_class
acclimation$Climate_reclass[acclimation$Climate_class %in% c("Bwk", "Csa")] <- 'arid_hot_med'
acclimation$Climate_reclass[acclimation$Climate_class %in% c("Bsk", "Bsh", "Csb", "Bwk")] <- 'semi_arid_warm_med'
acclimation$Climate_reclass[acclimation$Climate_class %in% c("Am", "Cfa", "Cfb", "Cfc", "Dfa", "Dfb", "Dfc", "Dfd", "Dwc")] <- 'humid'
# create folds
k_folds <- 5
folds <- createFolds(acclimation$Climate_reclass, k = k_folds, list = TRUE, returnTrain = FALSE)
#
#
data  <- data_TAS
# 
perf = data.frame(RMSE_train=double(), R2_train=double(), MAE_train=double(), RMSE_test=double(), R2_test=double(), MAE_test=double())
for (i in 1:40) {
  print(paste0('nodesize', i))
  y0.train <- c()
  y0.test  <- c()
  yobs.train <- c()
  yobs.test  <- c()
  for (j in 1:k_folds) {
    rf0 <- randomForest(formula = TAS ~ ., data=data[-folds[[j]], ], do.trace=FALSE, mtry=1, nodesize=i, ntree=500)
    y0  <- predict(rf0, data)    
    y0.train <- c(y0.train, y0[-folds[[j]]])
    y0.test  <- c(y0.test,  y0[folds[[j]]])
    #
    yobs.train <- c(yobs.train, data$TAS[-folds[[j]]])
    yobs.test  <- c(yobs.test, data$TAS[folds[[j]]])
  }
  #
  perf[i,1:3] <- postResample(pred=y0.train, obs=yobs.train)
  perf[i,4:6] <- postResample(pred=y0.test,  obs=yobs.test)
}
plot(perf[,1] - perf[,4])
plot(perf[,2] - perf[,5])
plot(perf[,3] - perf[,6])
plot(perf[,4])
plot(perf[,5])
plot(perf[,6])
#
##############################
# This is another way to set up cross-validation; but it can only tell the proper mtry. 
train_control <- trainControl(method = "cv",  # Cross-validation
                              number = 5)    # Number of folds
#
set.seed(985)
rf0 <- caret::train(form = TAS ~ ., data=data, method='rf', trControl = train_control, tuneLength = 5)
print(rf0)
rf0$bestTune    # mtry = 2
rf0$results


#####the model to use in manuscript, with bootstrapping
varImp_output <- data.frame(TRS_type = character(), var = character(), varImp = double())
#                          # low025 = double(), low050 = double(), med500 = double(), high950 = double(), high975 = double())
###########################################################For direct TAS###########################################################
library(rsample)
acclimation$climate_vegetation <- paste0(acclimation$Climate_class, '_', acclimation$IGBP)
data_strat <- acclimation[, c("TAS", "ELEV", "MATA", "LAI", "SOC", "climate_vegetation")]
strat_bootstrap <- bootstraps(data_strat, times = 200, strata = climate_vegetation)
# strat_bootstrap <- bootstraps(data_strat, times = 200, apparent = FALSE)
# This also gives the out of bag sampling using: first_oob <- assessment(strat_bootstrap$splits[[1]])  
# Access the first bootstrap sample
partial <- data.frame(var=character(), x=double(), y025=double(), y050=double(), y=double(), y950=double(), y975=double())
partial[1:51, 1]    <- "elev"
partial[52:102, 1]  <- "mat"
partial[103:153, 1] <- "lai"
partial[154:204, 1] <- "soc"
matrix.boot <- matrix(0, nrow=204, ncol=200)

for (i in 1:200) {
  print(i)
  data_sample <- analysis(strat_bootstrap$splits[[i]])
  data_sample <- data_sample[, 1:5]
  rf <- randomForest(formula = TAS ~ ., data=data_sample, do.trace=FALSE, mtry=1, nodesize=30, ntree=500, importance=TRUE)   # I have to add importance=TRUE here, to get type 1 RI values!
  # partial plot values
  tmp1 <- partialPlot(rf, pred.data=data, x.var="ELEV", plot=FALSE)
  tmp2 <- partialPlot(rf, pred.data=data, x.var="MATA", plot=FALSE)
  tmp3 <- partialPlot(rf, pred.data=data, x.var="LAI", plot=FALSE)
  tmp4 <- partialPlot(rf, pred.data=data, x.var="SOC", plot=FALSE)
  if (i == 1) {
    partial[1:51, 2]    <- tmp1$x
    partial[52:102, 2]  <- tmp2$x
    partial[103:153, 2] <- tmp3$x
    partial[154:204, 2] <- tmp4$x
  }
  matrix.boot[1:51, i]    <- tmp1$y
  matrix.boot[52:102, i]  <- tmp2$y
  matrix.boot[103:153, i] <- tmp3$y
  matrix.boot[154:204, i] <- tmp4$y
}

# confidence interval of bootstrap
for (i in 1:204) {
  partial[i, 3:7] <- quantile(matrix.boot[i,], c(0.025, 0.05, 0.5, 0.95, 0.975))
}
# standarize my x variable
partial$x_norm   <- partial$x
partial$x_norm[1:51] <- scale(partial$x[1:51])
partial$x_norm[52:102] <- scale(partial$x[52:102])
partial$x_norm[103:153] <- scale(partial$x[103:153])
partial$x_norm[154:204] <- scale(partial$x[154:204])
#
plot(partial$x[1:51], partial$y[1:51])
plot(partial$x[52:102], partial$y[52:102])
plot(partial$x[103:153], partial$y[103:153])
plot(partial$x[154:204], partial$y[154:204])
#

partial_output <- data.frame(TRS_type = 'TAS_direct', partial)
# write.csv(partial, 'data/partial_dependence_direct_TAS.csv', row.names = F)
#
#########################################################
####the last one best model for estimating direct TAS
data <- acclimation[, c("TAS", "ELEV", "MATA", "LAI", "SOC")]   # 
ri_var <- data.frame(v1 = double(), v2 = double(), v3 = double(), v4 = double())
for (i in 1:50) {
  rf0 <- randomForest(formula = TAS ~ ., data=data, do.trace=FALSE, mtry=1, nodesize=30, ntree=500, importance=TRUE)
  ri_var[i, 1:4] <- randomForest::importance(rf0, type=1, scale=TRUE) / sum(randomForest::importance(rf0, type=1, scale=TRUE)) * 100
}

y0  <- predict(rf0, data) 
print(postResample(pred=y0,  obs=data$TAS))
#
partialPlot(rf0, pred.data=data, x.var="ELEV")
partialPlot(rf0, pred.data=data, x.var="LAI")
#
varImpPlot(rf0, sort=TRUE, n.var=min(30, nrow(rf0$importance)),
           type=1, class=NULL, scale=TRUE,
           main='Relative importance')

varImp_output[1:4, 3] <- colMeans(ri_var)
varImp_output[1:4, 2] <- rownames(randomForest::importance(rf0, type=1, scale=TRUE))
varImp_output[1:4, 1] <- 'TAS_direct'

###########################################################For total TAS###########################################################
####the last one best model for estimating total TAS
data <- acclimation[, c("TAS_tot", "LAI", "ELEV", "MATA", "SOC")]   # "LAI" is slightly better than "GPP"
ri_var <- data.frame(v1 = double(), v2 = double(), v3 = double(), v4 = double())
for (i in 1:50) {
  rf0 <- randomForest(formula = TAS_tot ~ ., data=data, do.trace=FALSE, mtry=1, nodesize=30, ntree=500, importance=TRUE)
  ri_var[i, 1:4] <- randomForest::importance(rf0, type=1, scale=TRUE) / sum(randomForest::importance(rf0, type=1, scale=TRUE)) * 100
}

y0  <- predict(rf0, data) 
print(postResample(pred=y0,  obs=data$TAS_tot))

data <- data %>% mutate(y0 = y0, Climate = acclimation$Climate_class, IGBP = acclimation$IGBP, 
                        residual = acclimation$TAS_tot - y0, site_ID = acclimation$site_ID) 

data %>% group_by(IGBP) %>% summarise(residual_mean = mean(abs(residual)), n=n())
data %>% group_by(Climate) %>% summarise(residual_mean = mean(abs(residual)), n=n())
#
partialPlot(rf0, pred.data=data, x.var="LAI")
partialPlot(rf0, pred.data=data, x.var="ELEV")
partialPlot(rf0, pred.data=data, x.var="SOC")
partialPlot(rf0, pred.data=data, x.var="MATA")
#
varImpPlot(rf0, sort=TRUE, n.var=min(30, nrow(rf0$importance)),
           type=1, class=NULL, scale=TRUE,
           main='Relative importance')

varImp_output[5:8, 3] <- colMeans(ri_var)
varImp_output[5:8, 2] <- rownames(randomForest::importance(rf0, type=1, scale=TRUE))
varImp_output[5:8, 1] <- 'TAS_tot'

#################################################################################################
####get confidence intervals of the partial plot for total TAS
data_strat <- acclimation[, c("TAS_tot", "ELEV", "MATA", "LAI", "SOC", "climate_vegetation")] 
strat_bootstrap <- bootstraps(data_strat, times = 200, strata = climate_vegetation)
# This also gives the out of bag sampling using: first_oob <- assessment(strat_bootstrap$splits[[1]])
# Access the first bootstrap sample
partial <- data.frame(var=character(), x=double(), y025=double(), y050=double(), y=double(), y950=double(), y975=double())
partial[1:51, 1]    <- "elev"
partial[52:102, 1]  <- "mat"
partial[103:153, 1] <- "lai"
partial[154:204, 1] <- "soc"
matrix.boot <- matrix(0, nrow=204, ncol=200)
for (i in 1:200) {
  print(i)
  data_sample <- analysis(strat_bootstrap$splits[[i]])
  data_sample <- data_sample[, 1:5]
  rf <- randomForest(formula = TAS_tot ~ ., data=data_sample, do.trace=FALSE, mtry=1, nodesize=30, ntree=500, importance=TRUE)   # I have to add importance=TRUE here, to get type 1 RI values!
  # partial plot values
  tmp1 <- partialPlot(rf, pred.data=data, x.var="ELEV", plot=FALSE)
  tmp2 <- partialPlot(rf, pred.data=data, x.var="MATA", plot=FALSE)
  tmp3 <- partialPlot(rf, pred.data=data, x.var="LAI", plot=FALSE)
  tmp4 <- partialPlot(rf, pred.data=data, x.var="SOC", plot=FALSE)
  if (i == 1) {
    partial[1:51, 2]    <- tmp1$x
    partial[52:102, 2]  <- tmp2$x
    partial[103:153, 2] <- tmp3$x
    partial[154:204, 2] <- tmp4$x
  }
  matrix.boot[1:51, i]    <- tmp1$y
  matrix.boot[52:102, i]  <- tmp2$y
  matrix.boot[103:153, i] <- tmp3$y
  matrix.boot[154:204, i] <- tmp4$y
}

# confidence interval of bootstrap
for (i in 1:204) {
  partial[i, 3:7] <- quantile(matrix.boot[i,], c(0.025, 0.05, 0.5, 0.95, 0.975))
}
# standarize my x variable
partial$x_norm   <- partial$x
partial$x_norm[1:51] <- scale(partial$x[1:51])
partial$x_norm[52:102] <- scale(partial$x[52:102])
partial$x_norm[103:153] <- scale(partial$x[103:153])
partial$x_norm[154:204] <- scale(partial$x[154:204])
#
plot(partial$x[1:51], partial$y[1:51])
plot(partial$x[52:102], partial$y[52:102])
plot(partial$x[103:153], partial$y[103:153])
plot(partial$x[154:204], partial$y[154:204])
#
partial_output <- rbind(partial_output, data.frame(TRS_type = 'TAS_total', partial))

###------------------------------factors affecting apparent TAS-----------------
data <- acclimation[, c("TAS_app", "LAI", "ELEV", "MATA", "SOC")]   # "LAI" is slightly better than "GPP"
ri_var <- data.frame(v1 = double(), v2 = double(), v3 = double(), v4 = double())
for (i in 1:50) {
  rf0 <- randomForest(formula = TAS_app ~ ., data=data, do.trace=FALSE, mtry=1, nodesize=30, ntree=500, importance=TRUE)
  ri_var[i, 1:4] <- randomForest::importance(rf0, type=1, scale=TRUE) / sum(randomForest::importance(rf0, type=1, scale=TRUE)) * 100
}
  
y0  <- predict(rf0, data) 
print(postResample(pred=y0,  obs=data$TAS_app))
#       RMSE   Rsquared        MAE 
# 0.02352216 0.57249025 0.01600765
#
partialPlot(rf0, pred.data=data, x.var="LAI")
partialPlot(rf0, pred.data=data, x.var="ELEV")
partialPlot(rf0, pred.data=data, x.var="SOC")
partialPlot(rf0, pred.data=data, x.var="MATA")
varImpPlot(rf0, sort=TRUE, n.var=min(30, nrow(rf0$importance)),
           type=1, class=NULL, scale=TRUE,
           main='Relative importance')

varImp_output[9:12, 3] <- colMeans(ri_var)
varImp_output[9:12, 2] <- rownames(randomForest::importance(rf0, type=1, scale=TRUE))
varImp_output[9:12, 1] <- 'TAS_app'
# this is opposite: MATA, SOC, and LAI are most important. 
#################################################################################################
####get confidence intervals of the partial plot for total TAS
data_strat <- acclimation[, c("TAS_app", "ELEV", "MATA", "LAI", "SOC", "climate_vegetation")]
strat_bootstrap <- bootstraps(data_strat, times = 200, strata = climate_vegetation)
# This also gives the out of bag sampling using: first_oob <- assessment(strat_bootstrap$splits[[1]])
# Access the first bootstrap sample
partial <- data.frame(var=character(), x=double(), y025=double(), y050=double(), y=double(), y950=double(), y975=double())
partial[1:51, 1]    <- "elev"
partial[52:102, 1]  <- "mat"
partial[103:153, 1] <- "lai"
partial[154:204, 1] <- "soc"
matrix.boot <- matrix(0, nrow=204, ncol=200)
for (i in 1:200) {
  print(i)
  data_sample <- analysis(strat_bootstrap$splits[[i]])
  data_sample <- data_sample[, 1:5]
  rf <- randomForest(formula = TAS_app ~ ., data=data_sample, do.trace=FALSE, mtry=1, nodesize=30, ntree=500, importance=TRUE)   # I have to add importance=TRUE here, to get type 1 RI values!
  # partial plot values
  tmp1 <- partialPlot(rf, pred.data=data, x.var="ELEV", plot=FALSE)
  tmp2 <- partialPlot(rf, pred.data=data, x.var="MATA", plot=FALSE)
  tmp3 <- partialPlot(rf, pred.data=data, x.var="LAI", plot=FALSE)
  tmp4 <- partialPlot(rf, pred.data=data, x.var="SOC", plot=FALSE)
  if (i == 1) {
    partial[1:51, 2]    <- tmp1$x
    partial[52:102, 2]  <- tmp2$x
    partial[103:153, 2] <- tmp3$x
    partial[154:204, 2] <- tmp4$x
  }
  matrix.boot[1:51, i]    <- tmp1$y
  matrix.boot[52:102, i]  <- tmp2$y
  matrix.boot[103:153, i] <- tmp3$y
  matrix.boot[154:204, i] <- tmp4$y
}

# confidence interval of bootstrap
for (i in 1:204) {
  partial[i, 3:7] <- quantile(matrix.boot[i,], c(0.025, 0.05, 0.5, 0.95, 0.975))
}
# standarize my x variable
partial$x_norm   <- partial$x
partial$x_norm[1:51] <- scale(partial$x[1:51])
partial$x_norm[52:102] <- scale(partial$x[52:102])
partial$x_norm[103:153] <- scale(partial$x[103:153])
partial$x_norm[154:204] <- scale(partial$x[154:204])
#
plot(partial$x[1:51], partial$y[1:51])
plot(partial$x[52:102], partial$y[52:102])
plot(partial$x[103:153], partial$y[103:153])
plot(partial$x[154:204], partial$y[154:204])
#

partial_output <- rbind(partial_output, data.frame(TRS_type = 'TAS_app', partial))

write.csv(partial_output, 'plot_paper/partial_plot_TA.csv', row.names = F)
write.csv(varImp_output, 'plot_paper/varImp_plot_TA.csv', row.names = F)
