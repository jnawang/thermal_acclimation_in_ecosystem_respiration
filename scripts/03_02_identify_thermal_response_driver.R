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



