# This script estimates thermal acclimation strength in ecosystem respiration using moving window methods.
# Authors: Junna Wang, August and September, 2025
# Methods we tested:
# Step 1: determine the window size; 2 weeks ~ 2 month
# Step 2: for a specific year, do we need to  due to data limitation (such as no observations during a period). 
# Step 3: use all the data within that moving window to build an ER model. Use these parameter values as initial values of each year. 
# Step 4: find a method to estimate missed values. It is likely linear interpolation. 


# first try: I will use US-IB2 and US-Kon as an example
# to be determined, do I want to use nls to get initial values of brms
# how to save time? should we use overlapped windows or non-overlapping windows.
# We need to differentiate sites with water and sites without water; they use different equations. 

library(librarian)
shelf(dplyr, lubridate, gslnls, caret, performance, ggpubr, ggplot2, zoo, bayesplot)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
# dir_rawdata <- '/Users/junnawang/YaleLab/data_server/'
####################End Attention

site_info <- read.csv('data/site_info.csv')
feature_gs <- read.csv('data/growing_season_feature_EuropFlux.csv')
feature_gs_AmeriFlux <- read.csv('data/growing_season_feature_AmeriFlux.csv')
feature_gs <- rbind(feature_gs, feature_gs_AmeriFlux)

# outcome data frame
outcome <- data.frame(site_ID = character(), RMSE = double(), R2 = double(), control_year = double(), window_size = integer(), 
                      nwindow = integer(), TAS = double(), TASp = double())
outcome_siteyear <- data.frame()

priors_temp <- brms::prior("normal(2, 5)", nlpar = "C0", lb = 0, ub = 10) +
  brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 0.2) + 
  brms::prior("normal(-0.001, 0.1)", nlpar = "beta", lb = -0.01, ub = 0.0)

for (id in 1:nrow(site_info)) {
  # id = 25
  print(id)
  name_site <- site_info$site_ID[id]
  print(name_site)
  
  #-------------------------------------------DATA PREPARATION--------------------------
  # read data
  path <- file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.RDS'))
  if (file.exists(path)) {
    a_measure_night_complete <- readRDS(path)
    ac <- readRDS(file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.RDS')))
  }
  if (site_info$SWC_use[id] == 'YES') {
    a_measure_night_complete <- a_measure_night_complete %>% filter(!is.na(SWC))
  }
  
  # calculate daily daytime NEE and rolling average
  dt = 30  # minute
  if (ac$MINUTE[2] - ac$MINUTE[1] != 30) { dt = 60 }   # minutes 
  # unit is umol / m2 / s
  ac_day <- ac %>% filter(daytime) %>% group_by(YEAR, DOY) %>% 
    summarise(NEE_daytime1 = max(-sum(NEE_uStar_f, na.rm=T) * dt * 60 / 86400, 0.0), .groups = "drop")
  # get rolling average of the prior three days
  ac_day$NEE_daytime <- rollapply(ac_day$NEE_daytime1, width = 3, FUN = mean,
                                  align = "right", fill = ac_day$NEE_daytime1[1])
  
  # attach it to nighttime data
  if (!"NEE_daytime1" %in% names(a_measure_night_complete)) {
    a_measure_night_complete <- a_measure_night_complete %>% mutate(DOY_gpp = case_when(HOUR >= 12 ~ DOY,
                                                                                        HOUR < 12 & DOY == min(ac$DOY) ~ DOY, 
                                                                                        HOUR < 12 & DOY != min(ac$DOY) ~ DOY - 1)) %>%
      left_join(ac_day, by = c("YEAR", "DOY_gpp"="DOY")) %>% select(-"DOY_gpp") %>% filter(!is.na(NEE_daytime1))
  }
  
  # Deal with sites in southern hemisphere
  # add growing_year
  a_measure_night_complete <- a_measure_night_complete  %>% 
    mutate(growing_year = case_when(DOY <= 366 ~ YEAR, 
                                    TRUE ~ YEAR - 1))
  ac <- ac %>% 
    mutate(growing_year = case_when(DOY <= 366 ~ YEAR, 
                                    TRUE ~ YEAR - 1))
  
  # remove the growing_year with incomplete data
  # sites in southern hemisphere lose one year, because growing season crosses two years
  if (name_site %in% c('AU-Tum', 'ZA-Kru')) {
    a_measure_night_complete <- a_measure_night_complete %>% filter(between(growing_year, ac$YEAR[1], ac$YEAR[nrow(ac)] - 1))
    ac <- ac %>% filter(between(growing_year, ac$YEAR[1], ac$YEAR[nrow(ac)] - 1))
  }
  years <- sort(unique(a_measure_night_complete$growing_year))
  
  #---------------------------------DECIDE CONTROL YEAR, MOVING WINDOW SIZE, WINDOW NUMBERS------------------------------------
  # get start and end dates of growing season
  gStart <- feature_gs$gStart[feature_gs$site_ID == name_site]
  gEnd <- feature_gs$gEnd[feature_gs$site_ID == name_site]
  tStart <- feature_gs$tStart[feature_gs$site_ID == name_site]
  tEnd <- feature_gs$tEnd[feature_gs$site_ID == name_site]
  
  # decide control year: the year with growing-season TS closest to long-term mean. 
  ac_yearly_gs <- ac %>% filter(between(DOY, gStart, gEnd)) %>% filter(growing_year %in% years) %>% group_by(growing_year) %>% summarise(TS=mean(TS, na.rm=T))
  control_year <- ac_yearly_gs$growing_year[which.min(abs(ac_yearly_gs$TS - mean(ac_yearly_gs$TS)))]

  # determine moving window size and number of windows
  if (dt == 30) {
    nobs_threshold <- 100
  } else {
    nobs_threshold <- 60
  }

  # use uniform window size: 2 weeks
  window_size <- 14
  
  # use non-overlapping windows and determine number of windows for growing season; decide to use overlapping windows
  nwindow <- max(round((gEnd - gStart + 1) / window_size), 1)

  #------------------------------------PREPARE FOR formula, stprm, priors of ER models---------------------------------
  # for brm models
  frmu <- NEE ~ exp(alpha * TS + beta*TS^2) * C0
  param <- alpha+beta+C0 ~ 1
  priors <- priors_temp
  # for nls models
  # alpha and CO are both positive and beta is negative
  frmu_nls <- NEE ~ exp(exp(alpha_ln) * TS - exp(beta_ln)*TS^2) * (exp(C0_ln))
  stprm <- c(C0_ln = 0.7, alpha_ln = -2.99, beta_ln = -6.9)
  
  df_site_year_window <- data.frame(site_ID = character(), growing_year = integer(), window = character(), nobsv = integer(), extend_days = integer(),  
                                    alpha = double(), beta = double(), C0 = double(), Hs = double(), k2 = double(), TS = double(), ERref = double(), lnRatio = double())
  ER_obs_pred <- data.frame()

#------------------------------------FIT ER MODELS FOR EACH YEAR AND WINDOWS WITH BRM METHOD---------------------------------
  icount = 0
  for (iwindow in 1:nwindow) {
    window_start = gStart + window_size*(iwindow-1)
    window_end = min(gStart + window_size*iwindow, gEnd)
    
    ac_yearly_window <- ac %>% filter(between(DOY, window_start, window_end)) %>% group_by(growing_year) %>% summarise(TS=mean(TS, na.rm=T))
    if (!between(mean(ac_yearly_window$TS, na.rm=T), max(tStart, 2.0), tEnd) & nwindow > 1) {
      next
    }
    
    data <- a_measure_night_complete %>% filter(between(DOY, window_start, window_end))
    # skip a window if no enough data; this is for sites ('US-ICt', 'US-ICh', 'US-ICs', 'FI-Sod') with a period of the whole day is daytime. 
    if (nrow(data) < 100) { next }
    
    # try gsl_nls first because it is fast, and then update priors of brm models based on mod_nls
    # this can give abnormal initial values. 
    mod_nls <- try(gsl_nls(fn=as.formula(frmu_nls), data=data, start=stprm))
    if (!inherits(mod_nls, "try-error")) {
      priors$prior[priors$nlpar == 'alpha'] <- paste0("normal(", min(exp(coefficients(mod_nls)["alpha_ln"]), 0.2), ", 1.0)")
      priors$prior[priors$nlpar == 'beta'] <- paste0("normal(", -min(exp(coefficients(mod_nls)["beta_ln"]), 0.01), ", 0.1)")
      priors$prior[priors$nlpar == 'C0'] <- paste0("normal(", min(exp(coefficients(mod_nls)["C0_ln"]), 10), ", 5)")
    }

    # call the brm model to estimate parameters; this step takes much longer time.
    mod0 <- brms::brm(brms::bf(frmu, param, nl = TRUE),
                      prior = priors, data = data, iter = 2000, cores =4, chains = 4, backend = "cmdstanr",
                      control = list(adapt_delta = 0.95, max_treedepth = 15), refresh = 0) # , silent = 2
    # print(summary(mod0), digits = 3)
    
    # use this result as prior of each year
    priors$prior[priors$nlpar == 'alpha'] <- paste0("normal(", brms::fixef(mod0)["alpha_Intercept", "Estimate"], ", 1.0)")
    priors$prior[priors$nlpar == 'beta'] <- paste0("normal(", brms::fixef(mod0)["beta_Intercept", "Estimate"], ", 0.1)")
    priors$prior[priors$nlpar == 'C0'] <- paste0("normal(", brms::fixef(mod0)["C0_Intercept", "Estimate"], ", 5)")

    # get reference temperature, SWC, and NEEday of each window
    TSref <- mean(ac$TS[between(ac$DOY, window_start, window_end)], na.rm=T)
    NEEdayref <- mean(ac_day$NEE_daytime[between(ac_day$DOY, window_start, window_end)], na.rm=T)
    if (site_info$SWC_use[id] == 'YES') {
      SWCref <- mean(ac$SWC[between(ac$DOY, window_start, window_end)], na.rm=T)
      data_ref <- data.frame(TS=TSref, NEE_daytime=NEEdayref, SWC=SWCref)
    } else {
      data_ref <- data.frame(TS=TSref, NEE_daytime=NEEdayref)
    }
    
    #----loop through each year----    
    for (iyear in years) {
      print(paste(name_site, iwindow, iyear, sep='_'))
      icount = icount + 1
      df_site_year_window[icount, 1] <- name_site
      df_site_year_window$growing_year[icount] <- iyear
      df_site_year_window$window[icount] <- paste(window_start, window_end, sep='_')
      
      # determine if there are enough data for the regression of each year
      data_subset <- data[data$growing_year == iyear, ]
      # two rules are needed:
      # rule 1: total number of points > 100. 
      # rule 2: TSref is within the 0.025 and 0.975 quantiles. 
      # if the two rules are violated, extend window size.  # , 
      extend_days <- 0
      while(nrow(data_subset) < nobs_threshold | !between(TSref, quantile(data_subset$TS, 0.025, na.rm=T), quantile(data_subset$TS, 0.975, na.rm=T))) {
        extend_days <- extend_days + 3
        data_subset <- a_measure_night_complete %>% filter(growing_year == iyear) %>% 
          filter(between(DOY, window_start - extend_days, window_end + extend_days))
        # some conditions to break out to avoid dead loop
        if (nrow(data_subset) >= nobs_threshold & (window_end + extend_days) >= gEnd & TSref >= max(data_subset$TS, na.rm=T)) {break}
        if (extend_days > (gEnd - gStart) / 2.0) {break}   # max window size for some years = 91; 3 month
      }
      
      df_site_year_window[icount, 4] <- nrow(data_subset)
      df_site_year_window[icount, 5] <- extend_days
        
      mod <- brms::brm(brms::bf(frmu, param, nl = TRUE),
                       prior = priors, data = data_subset, iter = 1000, cores =4, chains = 4, backend = "cmdstanr", 
                       control = list(adapt_delta = 0.90, max_treedepth = 15), refresh = 0) # , silent = 2

      if (!inherits(mod, "try-error")) {
        np <- nuts_params(mod)
        n_divergent <- sum(subset(np, Parameter == "divergent__")$Value)
        failed_brm <- n_divergent > 0 
      } else {
        failed_brm <- TRUE
      }
      
      
      # if brm models fail or have divergent transitions, try another time
      if (failed_brm) {
        mod <- try(brms::brm(brms::bf(frmu, param, nl = TRUE),
                             prior = priors, data = data_subset, iter = 4000, cores =4, chains = 4, backend = "cmdstanr", 
                             control = list(adapt_delta = 0.98, max_treedepth = 15), refresh = 0)) # , silent = 2
      }
      
      # extract model results
      data_subset$NEE_pred <- fitted(mod)[, "Estimate"]
      ER_obs_pred <- rbind(ER_obs_pred, data_subset[between(data_subset$DOY, window_start, window_end), ])
      
      # print(plot(data_subset$TS, data_subset$NEE, main = paste(name_site, iwindow, iyear, sep = '_')))
      # lines(data_subset$TS, data_subset$NEE_pred)
      
      # model parameters
      df_site_year_window[icount, sub("_Intercept$", "", names(brms::fixef(mod)[, "Estimate"]))] <- brms::fixef(mod)[, "Estimate"]
      
      # average TS of moving window at each year
      df_site_year_window$TS[icount] <- ac_yearly_window$TS[ac_yearly_window$growing_year == iyear]
      
      # ER at reference temperature
      df_site_year_window$ERref[icount] <- fitted(mod, newdata=data_ref)[, "Estimate"]
      
      if (iyear == control_year) {
        ERref_control <- df_site_year_window$ERref[icount]
      }
    }
    # end of each year loop
    
    df_site_year_window$lnRatio[(icount-length(years) + 1):icount] <- log(df_site_year_window$ERref[(icount-length(years) + 1):icount] / ERref_control)
    
    # remove extreme values due to potentially large gaps
    x <- df_site_year_window$lnRatio[(icount-length(years) + 1):icount]
    outlier <- boxplot.stats(x, coef = 3)$out
    id.remove <- match(outlier[abs(outlier) > 1.5], x)
    if (length(id.remove) > 0) {
      df_site_year_window <- df_site_year_window[-(icount - length(years) + id.remove), ]
      icount <- icount - length(id.remove)
    }    

  }
  # end of each window
  outcome_siteyear <- rbind(outcome_siteyear, df_site_year_window)
  
  tmp <- df_site_year_window %>% group_by(growing_year) %>% summarise(TS = mean(TS), lnRatio = sum(lnRatio*ERref)/sum(ERref)) # %>% left_join(ac_yearly_gs, by="growing_year")
  
  plot <- df_site_year_window %>% ggplot(aes(x=TS, y=lnRatio, col=window)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    labs(title = name_site)
  ggsave(paste0('graphs/', name_site, '_TAS_temp.png'))
  print(plot)
  
  outcome[id, "site_ID"] <- name_site
  outcome[id, c("RMSE", "R2")] <- postResample(pred = ER_obs_pred$NEE_pred, obs = ER_obs_pred$NEE)[1:2]
  outcome[id, c("control_year", "window_size", "nwindow")] <- c(control_year, window_size, nwindow)
  # outcome[id, c("TAS", "TASp")] <- summary(lm(data=tmp, lnRatio ~ TS))$coefficients[2, c(1, 4)]
  if (length(unique(df_site_year_window$window)) > 1) {
    outcome[id, c("TAS", "TASp")] <- summary(lm(data=df_site_year_window, lnRatio ~ TS + window, na.action = na.exclude))$coefficients[2, c(1, 4)]
  } else {
    outcome[id, c("TAS", "TASp")] <- summary(lm(data=df_site_year_window, lnRatio ~ TS, na.action = na.exclude))$coefficients[2, c(1, 4)]
  }
}

# end of each site
write.csv(outcome, 'data/outcome_window_temp.csv', row.names = F)
write.csv(outcome_siteyear, 'data/outcome_siteyear_window_temp.csv', row.names = F)

