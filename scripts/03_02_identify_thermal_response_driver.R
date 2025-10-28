# This script identifies drivers of apparent and actual thermal response strength
# Authors: Junna Wang, October, 2025

library(librarian)
shelf(dplyr, ggplot2, corrplot, MuMIn, fBasics, car, lme4, vip, randomForest, caret, VSURF)
rm(list=ls())

#
acclimation <- read.csv('data/acclimation_data.csv')

#------------------------do some explorations first-----------------------------
# intrinsic TRS--more related to vegetation type
summary(lm(data=acclimation, TAS ~ IGBP))
summary(lm(data=acclimation, TAS ~ Climate_class))
summary(aov(TAS ~ IGBP, data=acclimation))
summary(aov(TAS ~ Climate_class, data=acclimation))

# apparent TRS--strongly affected by climate
summary(lm(data=acclimation, TAS_tot ~ IGBP))
summary(lm(data=acclimation, TAS_tot ~ Climate_class))

#------------------------correlation among predictor variables------------------
data.cor <- acclimation[, c(4, 7:20, 22)]
cor(data.cor)

#----------------------------try simple linear regression-----------------------
summary(lm(data=acclimation, TAS ~ MATA + ELEV + NEE_day + SOC))          # Only productivity is most important. 
summary(lm(data=acclimation, TAS_tot ~ MATA + ELEV + NEE_day + SOC))      # The first 3 are all important. 

options(na.action = na.omit)
TAS_tot_data <- acclimation[, c(4, 7:20)]
model <- lm(TAS_tot ~ ., data = TAS_tot_data)   # climate class is strongly correlated to MATA [,-c()]
summary(model)
dredge(model)
# elevation, MAT, NEE, NEE_day, NEE_night


TAS_data <- acclimation[, c(4, 7:19, 22)]
model <- lm(TAS ~ ., data = TAS_data)   # climate class is strongly correlated to MATA [,-c()]
summary(model)
dredge(model)
# NEE, IAT 

#-------------------------------------------------------------------------------
#----------------------------------try random forests---------------------------
set.seed(222)
acclimation.use <- acclimation[, c("TAS", "ELEV", "MATA", "LAI", "SOC")]
# acclimation.use$TAS <- acclimation$TAS_tot
# calculate vif of the four predictors
lm_model <- lm(data=acclimation.use, TAS ~ .)
car::vif(lm_model)
# vif: variance inflation factors
#     ELEV     MATA      Lai  SOC 
# 1.114907 1.306338 1.178243 1.172405
#
# stratify the data, and create folds
acclimation$Climate_reclass <- acclimation$Climate_class
acclimation$Climate_reclass[acclimation$Climate_class %in% c("Bwk", "Csa")] <- 'arid_hot_med'
acclimation$Climate_reclass[acclimation$Climate_class %in% c("Bsk", "Bsh", "Csb")] <- 'semi_arid_warm_med'
acclimation$Climate_reclass[acclimation$Climate_class %in% c("Am", "Cfa", "Cfb", "Dfa", "Dfb", "Dfc", "Dfd", "Dwc")] <- 'humid'
# create folds
k_folds <- 5
folds <- createFolds(acclimation$Climate_reclass, k = k_folds, list = TRUE, returnTrain = FALSE)
#
#
# var.use[1:22] , "Lai" , "MAP", "MATA",  # do not incorporate nitrogen because it is hard to justify the highly non-linear results 
data  <- acclimation.use # %>% filter(!is.na(age))  # adding ages also does not make sense to the results
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
################################################### 
# This is another way to set up cross-validation; but it can only tell the proper mtry. 
train_control <- trainControl(method = "cv",  # Cross-validation
                              number = 5)    # Number of folds
#
set.seed(985)
rf0 <- caret::train(form = TAS ~ ., data=data, method='rf', trControl = train_control, tuneLength = 5)
print(rf0)
rf0$bestTune
rf0$results

#####the model to use, with bootstrapping
library(rsample)
data_strat <- acclimation[, c("TAS", "ELEV", "MATA", "LAI", "SOC", "Climate_reclass")]
strat_bootstrap <- bootstraps(data_strat, times = 200, strata = "Climate_reclass")
# This also gives the out of bag sampling using: first_oob <- assessment(strat_bootstrap$splits[[1]])
# Access the first bootstrap sample
varImp  <- data.frame(elev=double(), mat=double(), lai=double(), ph=double())
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
  rf <- randomForest(formula = TAS ~ ., data=data_sample, do.trace=FALSE, mtry=1, nodesize=24, ntree=500, importance=TRUE)   # I have to add importance=TRUE here, to get type 1 RI values!
  Imp_tmp <- importance(rf, type=1, scale=TRUE)[1:4]
  varImp[i, 1:4] <- Imp_tmp / sum(Imp_tmp)
  # type=2: the decrease in node impurity; The node impurity is measured by the Gini index; May have some problems. so I use type I: permutation-based MSE reduction. 
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
# confidence interval of RI; the three are almost equally important
for (i in 1:4) {
  print(quantile(varImp[,i], c(0.025, 0.05, 0.5, 0.95, 0.975)))
}
#       2.5%        5%       50%       95%     97.5% 
# ele: 0.1886710 0.2129134 0.2940669 0.3742908 0.3956409  (1)
# mat: 0.1453366 0.1549683 0.2172829 0.2881334 0.2987378  (3)
# lai: 0.2187312 0.2311713 0.2988688 0.3824310 0.4117751  (2) 
# soc: 0.1053378 0.1169086 0.1846725 0.2545087 0.2620789  (4)

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
# write.csv(partial, '/Users/jw2946/Documents/stability_project/manuscripts/plot/partial_plot.csv', row.names = F)
#
#################################################################################################
####the last best model
rf0 <- randomForest(formula = TAS ~ ., data=data, do.trace=FALSE, mtry=1, nodesize=24, ntree=500, importance=TRUE)
print(rf0)
y0  <- predict(rf0, data) 
print(postResample(pred=y0,  obs=data$TAS))
#       RMSE   Rsquared        MAE 
# 0.05788558 0.63570658 0.04197571
# out of bag variation explained: 9.27% 
#
partialPlot(rf0, pred.data=data, x.var="ELEV")
partialPlot(rf0, pred.data=data, x.var="LAI")
# partialPlot(rf0, pred.data=data, x.var="MAP")
partialPlot(rf0, pred.data=data, x.var="MATA")
partialPlot(rf0, pred.data=data, x.var="SOC")
# partialPlot(rf0, pred.data=data, x.var="soc")
# partialPlot(rf0, pred.data=data, x.var="age")
#
varImpPlot(rf0, sort=TRUE, n.var=min(30, nrow(rf0$importance)),
           type=1, class=NULL, scale=TRUE,
           main='Relative importance')




