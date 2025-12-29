# This script estimates the effects of total thermal responses on future ecosystem respiration
# Authors: Junna Wang, October, 2025

# This script takes ~3 hours to run.

library(librarian)
shelf(dplyr, terra, ggplot2, caret, performance, zoo, bayesplot, brms, gslnls, lubridate)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
# dir_rawdata <- '/Users/junnawang/YaleLab/data_server/'
####################End Attention

options(na.action = "na.omit")
#
# read data used for this module
feature_gs <- read.csv(file.path('data', 'growing_season_feature_EuropFlux.csv'))
feature_gs_AmeriFlux <- read.csv(file.path('data', 'growing_season_feature_AmeriFlux.csv'))
feature_gs <- rbind(feature_gs, feature_gs_AmeriFlux)
#
acclimation <- read.csv(file.path("data", "acclimation_data.csv"))
#
Tmin_month <- read.csv(file.path('data', 'Tmin_month_ssp245_wc.csv'))
acclimation     <- acclimation %>% left_join(Tmin_month, by='site_ID')

# add 6 column about respiration
acclimation$TS_TA   <- 0            # conversion between TS and TA
acclimation$TSmin_c <- 0            # monthly min soil temperature change
acclimation$TSmin_c_gs <- 0         # gs: growing season
acclimation$NEE_night_mod_p <- 0
acclimation$NEE_night_mod_f <- 0
acclimation$NEE_night_mod_fa <- 0

#
priors_temp <- brms::prior("normal(2, 5)", nlpar = "C0", lb = 0, ub = 10) +
  brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 0.2) + 
  brms::prior("normal(-0.001, 0.1)", nlpar = "beta", lb = -0.01, ub = 0.0)

frmu <- NEE ~ exp(alpha * TS + beta*TS^2) * C0
param <- alpha+beta+C0 ~ 1
priors <- priors_temp

#
files  <- list.files(path = file.path(dir_rawdata, "RespirationData"), pattern = '_ac.csv$', full.names = FALSE)
# air and soil temperature patterns
# for (i in 1:length(files)) {
for (i in 56:length(files)) {
  ###
  # i = 74
  name_site <- substring(files[i], 1, 6)
  print(paste0(i, name_site))
  #
  iacclimation <- which(acclimation$site_ID==name_site)
  #
  ac <- read.csv(file.path(dir_rawdata, "RespirationData", files[i]))
  #
  # use only the years with qualified data
  a_measure_night_complete <- read.csv(file.path(dir_rawdata, "RespirationData", paste0(name_site, "_nightNEE.csv")))
  good_years <- unique(a_measure_night_complete$YEAR)
  # print(good_years)
  
  #
  gStart <- feature_gs$gStart[feature_gs$site_ID == name_site]
  gEnd <- feature_gs$gEnd[feature_gs$site_ID == name_site]
  
  ############I should use control year, and this control year should have enough data during every period#######
  yearly_TSgs <- a_measure_night_complete %>% filter(between(DOY, gStart, gEnd)) %>% group_by(YEAR) %>% summarise(n = n(), TSgs=mean(TS))
  yearly_TSgs <- yearly_TSgs %>% mutate(close_mean = abs(TSgs - median(yearly_TSgs$TSgs))) %>% arrange(close_mean)
  control_year <- yearly_TSgs$YEAR[which(yearly_TSgs$n[1:3] == sort(yearly_TSgs$n[1:3])[2])]
  # at some sites, the automatically selected years may not work the best, choose another year. 
  if (name_site == 'AU-Tum') {
    control_year = 2001
  } else if (name_site == 'CA-Cbo') {
    control_year = 2006
  } else if (name_site == 'CH-Fru') {
    control_year = 2008
  } else if (name_site == 'US-Wkg') {
    control_year = 2008
  } else if (name_site == 'ES-LJu') {
    control_year = 2010
  } 
  
  
  # control_year = outcome$control_year[outcome$site_ID == name_site]
  a_measure_night_complete_control <- a_measure_night_complete %>% filter(YEAR == control_year)
  # plot(a_measure_night_complete$TS[a_measure_night_complete$YEAR == 2009], a_measure_night_complete$NEE[a_measure_night_complete$YEAR == 2009], main = paste0(i, name_site))
  
  # recalculate TS_TA relationship using TA threshold
  if (name_site!='GF-Guy') {
    tmp <- ac %>% filter(TA > 0)     # only use the data with TA > 0C; use another range to calculate the coefficient
    acclimation$TS_TA[iacclimation] <- summary(lm(data=tmp, TS~TA))$coefficients[2,1]
  } else {
    acclimation$TS_TA[iacclimation] <- 0.311303822     # tropical site: its values are obtained using observed temperature only.
  }
  #
  night_pattern <- ac %>% filter(YEAR %in% good_years & !daytime) %>% group_by(DOY, HOUR, MINUTE, MONTH) %>% 
    summarise(TAp=mean(TA, na.rm=T), TSp=mean(TS, na.rm=T), SWC=mean(SWC, na.rm=T), .groups = 'drop')
  # remove duplicate rows
  night_pattern  <- night_pattern[!duplicated(night_pattern[,c('DOY', 'HOUR', 'MINUTE')]), ]
  # deal with occasional TA missing cases
  if (sum(is.na(night_pattern$TAp)) > 0) {
    night_pattern$TAp <- na.approx(night_pattern$TAp)
  }
  # deal with occasional TS missing cases
  if (sum(is.na(night_pattern$TSp)) > 0) {
    mod_lm <- lm(data=night_pattern, TSp ~ TAp)
    TSp_pred <- predict(mod_lm, night_pattern)
    night_pattern$TSp[is.na(night_pattern$TSp)] <- TSp_pred[is.na(night_pattern$TSp)]
  }
  #
  temp_change   <- data.frame(MONTH=1:12, TAmnc=as.numeric(acclimation[iacclimation, which(colnames(acclimation) %in% (paste0('Tmin', 1:12)))]))
  #
  night_pattern <- night_pattern %>% left_join(temp_change, by="MONTH")
  #
  night_pattern$TAf <- night_pattern$TAp + night_pattern$TAmnc
  night_pattern$TSf <- night_pattern$TSp + night_pattern$TAmnc * acclimation$TS_TA[iacclimation]
  #
  night_pattern$NEEp <- NA
  night_pattern$NEEf <- NA
  
  
  # use uniform window size: 2 weeks
  window_size <- 14
  
  # use non-overlapping windows and determine number of windows for growing season; decide to use overlapping windows
  nwindow <- max(round((gEnd - gStart + 1) / window_size), 1)  
  
  ER_obs_pred <- data.frame()
  for (iwindow in 1:nwindow) {
    window_start = gStart + window_size*(iwindow-1)
    window_end = min(gStart + window_size*iwindow, gEnd)
    if (iwindow == nwindow & window_end != gEnd) {
      window_end <- gEnd
    }
    
    # Skip a window if no need to estimate ER due to all daytime; 
    # This is the case for sites ('US-ICt', 'US-ICh', 'US-ICs', 'FI-Sod') with a period of the whole day is daytime. 
    # Data measurement sometime is not reliable, when TS is lower than 1C or 2C.
    id <- which(between(night_pattern$DOY, window_start, window_end) & night_pattern$TSp > 1.0)
    if ( length(id) == 0 ) { next }

    # extend the range to get TS~ER curve, because we need to predict for warming future
    data <- a_measure_night_complete_control %>% filter(between(DOY, window_start - 7, window_end + 7) & TS > 1.0)
    if (name_site %in% c('US-ICt', 'US-ICh', 'US-ICs')) {
      nobs_threshold = 90
    } else {
      nobs_threshold = 200
    }
    
    # if not enough observed data in control year, extend the window 
    extend_days <- 0
    while(nrow(data) < nobs_threshold) {
      extend_days <- extend_days + 7
      data <- a_measure_night_complete_control %>% filter(between(DOY, window_start - extend_days, window_end + extend_days) & TS > 1.0)
      # some conditions to break out to avoid dead loop
      if (extend_days > (gEnd - gStart) / 2.0) {break}
    }

    # call the brm model to estimate parameters; this step takes much longer time.
    mod <- brms::brm(brms::bf(frmu, param, nl = TRUE),
                     prior = priors, data = data, iter = 2000, cores = 4, chains = 4, backend = "cmdstanr",
                     control = list(adapt_delta = 0.95, max_treedepth = 15), refresh = 0) # , silent = 2
    # print(summary(mod), digits = 4)

    df_accurate <- data.frame(NEE=data$NEE, TS=data$TS)
    NEE_pred <- fitted(mod, newdata=df_accurate)[, "Estimate"]
    plot(df_accurate$TS, df_accurate$NEE, main = paste0(iwindow, name_site))
    ER_obs_pred <- rbind(ER_obs_pred, data.frame(pred=NEE_pred, obs=df_accurate$NEE))

    df_present <- data.frame(TS=night_pattern$TSp[id])
    night_pattern$NEEp[id] <- fitted(mod, newdata=df_present)[, "Estimate"]

    df_future <- data.frame(TS=night_pattern$TSf[id])
    night_pattern$NEEf[id] <- fitted(mod, newdata=df_future)[, "Estimate"]
    
    #----------------------------------------------------------------------------------------------------------
    # # Alternative model to estimate parameters, this method is fast, but need more data to get reliable results.
    # df_accurate <- data.frame(NEE=data$NEE, TS=data$TS)
    # mod2 <- gsl_nls(fn=frmu_nls, data=df_accurate, start=stprm)
    # print(summary(mod2))
    # 
    # NEE_pred <- predict(mod2, newdata=df_accurate)
    # plot(df_accurate$TS, df_accurate$NEE, main = paste0(name_site, iwindow))
    # ER_obs_pred <- rbind(ER_obs_pred, data.frame(pred=NEE_pred, obs=df_accurate$NEE))
    # #
    # # predict current night respiration
    # df_present <- data.frame(TS=night_pattern$TSp[id])
    # night_pattern$NEEp[id] <- predict(mod2, newdata=df_present)
    # #
    # # predict for future respiration
    # df_future <- data.frame(TS=night_pattern$TSf[id])
    # night_pattern$NEEf[id] <- predict(mod2, newdata=df_future)
  }
  # check performance
  # print(postResample(pred=ER_obs_pred$pred, obs=ER_obs_pred$obs))
  
  # get the required values
  acclimation$TSmin_c[iacclimation] <- mean(temp_change$TAmnc) * acclimation$TS_TA[iacclimation]
  #
  # NEE during growing season only
  acclimation$NEE_night_mod_p[iacclimation]  <- mean(night_pattern$NEEp[night_pattern$DOY >= gStart & night_pattern$DOY <= gEnd], na.rm=T)
  acclimation$NEE_night_mod_f[iacclimation]  <- mean(night_pattern$NEEf[night_pattern$DOY >= gStart & night_pattern$DOY <= gEnd], na.rm=T)
  NEE_gs <- acclimation$NEE_night_mod_f[iacclimation]
  acclimation$TSmin_c_gs[iacclimation] <- mean(night_pattern$TAmnc[night_pattern$DOY >= gStart & night_pattern$DOY <= gEnd], na.rm=T) * acclimation$TS_TA[iacclimation]
  acclimation$NEE_night_mod_fa[iacclimation] <- NEE_gs * exp(acclimation$TAS_tot[iacclimation] * acclimation$TSmin_c_gs[iacclimation])
}
#
write.csv(acclimation, file=file.path('data', 'acclimation_data_future_ssp245_56_110.csv'), row.names = F)

# 
(mean(acclimation$NEE_night_mod_f) - mean(acclimation$NEE_night_mod_fa)) / (mean(acclimation$NEE_night_mod_f) - mean(acclimation$NEE_night_mod_p))

# if we only consider the sites with significant effects, the ratio is about 0.1974396
acclimation <- acclimation %>% mutate(NEE_night_mod_fa_sig = ifelse(TAS_totp > 0.05, NEE_night_mod_f, NEE_night_mod_fa))
(mean(acclimation$NEE_night_mod_f) - mean(acclimation$NEE_night_mod_fa_sig)) / (mean(acclimation$NEE_night_mod_f) - mean(acclimation$NEE_night_mod_p))

