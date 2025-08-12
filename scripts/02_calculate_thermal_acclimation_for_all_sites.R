# This script calculates the strength of thermal acclimation for all study sites
# Junna Wang
# 7/24/2025
# pay attention to FR-Pue, sensitive to data in 2014; IT-Noe: they really likely have enhancing effects. 
# CZ-RAJ: change TAS sign, when considering 3 more years.
# it is important to see how adding 2023, 2024 data to impact NL-Loo results, particularly after considering GPP effects.  
# IL-Yat: last version TS data is wrong. The current TS data make sense. 
# FR-FBn: new data make this site show enhancing responses.
# 

library(librarian)
shelf(dplyr, lubridate, gslnls, caret, performance, ggpubr, ggplot2)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
####################End Attention

site_info <- read.csv('data/site_info.csv')
feature_gs <- read.csv('data/growing_season_feature.csv')
feature_gs_AmeriFlux <- read.csv('data/growing_season_feature_AmeriFlux.csv')
feature_gs <- rbind(feature_gs, feature_gs_AmeriFlux)

# a data frame to store model performance
ER_model_performance <- data.frame(RMSE=double(), Rsquared=double(), MAE=double(), C_pool=double(), alpha=double(), beta=double(), Hs=double())
# nrow(site_info)
for (id in 87:87) {
  # id = 110
  print(id)
  name_site <- site_info$site_ID[id]
  print(name_site)
  
  path <- file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.RDS'))
  if (file.exists(path)) {
    a_measure_night_complete <- readRDS(path)
    ac <- readRDS(file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.RDS')))
  } else {
    next
  }
  
  data <- a_measure_night_complete
# ---- get a TS ~ ER model for all the data -----
  # this result is used as a test of model performance
  # check data quality first
  plot(data$TS, data$NEE)
  plot(data$DOY, data$NEE)
  if (site_info$SWC_use[id] == 'YES') {
    plot(data$SWC, data$NEE)
  }
  
  # model: quadratic exponential (!the best model!)
  # Negative Hs, so I remove the soil moisture effects.
  if (site_info$SWC_use[id] == 'YES') {
    frmu   <- NEE ~ C_pool * exp(alpha * TS + beta*TS^2) * SWC / (Hs + SWC)     # in this formula: C_pool, alpha, beta, and Hs are unclear parameters. I know Cpool and Hs should be positive;
    stprm  <- c(C_pool=1, alpha=0.1, beta=-0.01, Hs = 20)
  } else {
    frmu   <- NEE ~ C_pool * exp(alpha * TS + beta*TS^2) 
    stprm  <- c(C_pool=1, alpha=0.1, beta=-0.01)
  }
  
  # some sites need unique initial parameter values to make models converge
  if (name_site == 'ZA-Kru') {
    if (site_info$SWC_use[id] == 'YES') {
      stprm  <- c(C_pool=10, alpha=0.1, beta=-0.01, Hs = 100)
    } else {
      stprm  <- c(C_pool=10, alpha=0.1, beta=-0.01)
    }
  }
  
  if (name_site == 'GF-Guy') {
    if (site_info$SWC_use[id] == 'YES') {
      stprm  <- c(C_pool=1, alpha=0.1, beta=-0.001, Hs = 20)
    } else {
      stprm  <- c(C_pool=1, alpha=0.1, beta=-0.001)
    }    
  }
  
  df_accurate <- data.frame(NEE=data$NEE, TS=data$TS, SWC=data$SWC)
  mod2 <- gsl_nls(fn=frmu, data=df_accurate, start=stprm)
  summary(mod2)
  NEE_pred <- predict(mod2, newdata=df_accurate)
  plot(df_accurate$TS, df_accurate$NEE)

  # check performance
  ER_model_performance[id, 1:3] <- postResample(pred=NEE_pred, obs=df_accurate$NEE)
  ER_model_performance[id, 4:6] <- coefficients(mod2)[1:3]
  
  # get Hs if soil water is modeled
  if (site_info$SWC_use[id] == 'YES') {
    Hs <- coefficients(mod2)[4]
    ER_model_performance[id, 7] <- Hs
  }


# ---- get Ts-ER curve for each year---- 
  # Add growing_year for sites in southern hemisphere
  a_measure_night_complete <- a_measure_night_complete  %>% 
    mutate(growing_year = case_when(DOY <= 366 ~ YEAR, 
                                    TRUE ~ YEAR - 1))

  iStart = a_measure_night_complete$YEAR[1]
  iEnd   = a_measure_night_complete$YEAR[nrow(a_measure_night_complete)]
  
  # sites in southern hemisphere loss one year, because growing season crosses two years
  if (name_site %in% c('AU-Tum', 'ZA-Kru', "BR-Sa1")) {
    iEnd <- iEnd - 1
  }
  
  coef <- data.frame(year=double(), C_pool=double(), alpha=double(), beta=double())
  perf <- data.frame(year=double(), RMSE=double(), Rsquared=double(), MAE=double())
  
  # decide which formula to use based on if water is used. 
  # model 2: quadratic exponential (!the best model!)
  if (site_info$SWC_use[id] == 'YES') {
    frmu   <- NEE ~ C_pool * exp(exp(alpha) * TS + beta*TS^2)  * SWC / (SWC+Hs)
  } else {
    frmu   <- NEE ~ C_pool * exp(exp(alpha) * TS + beta*TS^2)
  }  
  stprm  <- c(C_pool=20, alpha=-1, beta=-0.01)            
  
  # These initial values need to be adjusted for some sites. 
  if (name_site == 'GF-Guy') {
    stprm  <- c(C_pool=20, alpha=-1, beta=-0.001)
  }
  
  for (i in iStart:iEnd) {     
    data <- a_measure_night_complete %>% filter(growing_year==i)   #
    if (site_info$SWC_use[id] == 'YES') {
      data <- data %>% filter(!is.na(SWC))
    }

    if (nrow(data) <= 10) {
      next
    }
    
    print(i)
    df_accurate <- data.frame(NEE=data$NEE, TS=data$TS, SWC=data$SWC)      
    # mod2 <-try(nlxb(formula=frmu, start=stprm, data=df_accurate, lower=c(0, 0, -1000), upper=c(1000, 1000, 1000)))
    mod2 <- gsl_nls(fn=frmu, data=df_accurate, start=stprm)
    summary(mod2)
    NEE_pred <- predict(mod2, newdata=df_accurate)
    plot(df_accurate$TS, df_accurate$NEE, main=i)
    lines(df_accurate$TS, NEE_pred)
    # check performance
    performance_aicc(mod2)
    coef[i-iStart+1,] <- c(i, coefficients(mod2)[1], exp(coefficients(mod2)[2]), coefficients(mod2)[3])
    perf[i-iStart+1,] <- c(i, postResample(pred=NEE_pred, obs=df_accurate$NEE))  
  }
  
  # Remember to check performance before moving on!
  coef$Topt  <- - coef$alpha/2/coef$beta
  coef$Rpeak <- coef$C_pool * exp(-coef$alpha^2 / 4 / coef$beta)  


# ---- prepare for data for calculating strength of thermal acclimation ----
  gStart <- feature_gs$gStart[feature_gs$site_ID==name_site]
  gEnd <- feature_gs$gEnd[feature_gs$site_ID==name_site]
  tStart <- as.integer(feature_gs$tStart[feature_gs$site_ID==name_site] * 10) + 1
  tEnd <- as.integer(feature_gs$tEnd[feature_gs$site_ID==name_site] * 10) + 1
  
  ac <- ac %>% mutate(growing_year = case_when(DOY <= 366 ~ YEAR, 
                                               TRUE ~ YEAR - 1))
  NEE_annual <- ac %>% filter(between(growing_year, iStart, iEnd)) %>% group_by(growing_year) %>% summarise(TA=mean(TA, na.rm=T), TS=mean(TS, na.rm=T), SWC=mean(SWC, na.rm=T)) 
  NEE_grow   <- ac %>% filter(between(growing_year, iStart, iEnd)) %>% filter(DOY >= gStart & DOY <= gEnd) %>% group_by(growing_year) %>% summarise(TA=mean(TA, na.rm=T), TS=mean(TS, na.rm=T), SWC=mean(SWC, na.rm=T)) 
  
  result <- data.frame(TS=double(), NEE=double(), year=integer())
  TS    = seq(0, ceiling(tEnd/10)+1, 0.1)
  nTS   = length(TS)
  nyear = nrow(coef)
  for (i in 1:nyear) {
    #
    result[((i-1)*nTS+1):(i*nTS), 3] <- i + iStart - 1
    result[((i-1)*nTS+1):(i*nTS), 1] <- TS
    #
    # if (i %in% c()) {next}
    if (is.na(perf$Rsquared[i])) {
      next
    }
    if (perf$Rsquared[i] < 0.0) {             # I have to remove bad performance years. mean(perf$Rsquared, na.rm=T) * 0.5 
      next
    }
    alpha  = coef$alpha[i]
    beta   = coef$beta[i]
    C_pool = coef$C_pool[i]
    result[((i-1)*nTS+1):(i*nTS), 2] <- C_pool * exp(alpha * TS + beta * TS * TS)
  }
  
  
  #
  p <- ggplot(result, aes(x=TS, y=NEE, col=as.factor(year))) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = tStart/10) +
    geom_vline(xintercept = tEnd/10)
  print(p)

  
  # growing season 25-31 C; using a standard method, integration over the period. 
  # get the average NEE during the growing season.
  result_row <- tidyr::spread(result, key=year, value=NEE)
  # only choose the data with temperature between 11 and 18 C
  Ravg_GW    <- colSums(result_row[tStart:tEnd, 2:(nyear+1)]) * 0.1
  Ravg_GW    <- Ravg_GW / mean(Ravg_GW, na.rm=T) 
  
  #
  df <- merge(NEE_grow, data.frame(growing_year= as.integer(names(Ravg_GW)), alpha = Ravg_GW), by='growing_year')
  p <- ggplot(df, aes(x=TS, y=alpha)) + 
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(color='green') + 
    labs(x='Mean soil temperature (C)', y='normalized night respiration under growth temperature')
  # linear regression
  print(summary(lm(data=df, alpha~TS)))
  print(p)
  
  ggplot(df, aes(x=TA, y=alpha)) + 
    geom_point() +
    geom_smooth(method = lm) +
    stat_cor(color='green') + 
    labs(x='Mean air temperature (C)', y='normalized night respiration under growth temperature')
  # linear regression
  summary(lm(data=df, alpha~TA))
  
  #-----------------------------------------------------------
  # how about the relation during growing season?
  tmp <- ac %>% filter(TS > 2)     # only use the data with TS > 2C
  TS_TA <- summary(lm(data=tmp, TS~TA))$coefficients[2,1]
  
  # write important results into an output file
  outcome <- list(site_ID=name_site, gStart=gStart, gEnd=gEnd, tStart=tStart, tEnd=tEnd, TS_TA=TS_TA,
                  NEE_grow=NEE_grow, NEE_annual=NEE_annual, result=result_row, coef=coef, perf=perf)
  
  saveRDS(outcome, file=file.path(dir_rawdata, "RespirationData", paste0(name_site, '_outcome.RDS')))

}
