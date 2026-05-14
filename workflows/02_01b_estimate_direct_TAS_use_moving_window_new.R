# This script estimates thermal acclimation strength in ecosystem respiration using moving window methods.
# Authors: Junna Wang, August and September, 2025
# Methods we tested:
# Step 1: determine the window size; 2 weeks ~ 4 weeks
# Step 2: for a specific year, do we need to skip due to data limitation (such as no observations during a period). 
# Step 3: use all the data within that moving window to build an ER model. Use these parameter values as initial values of each year. 
# Step 4: find a method to estimate missed values. It is likely linear interpolation. 

# It takes about 2 days to finish the estimate of TAS; but on other computers it could take one week, depending on number of cores used.  
# please install 'cmdstanr' before you run this script. 

library(librarian)
shelf(dplyr, lubridate, gslnls, caret, performance, ggpubr, ggplot2, zoo, bayesplot, brms, nlme)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
dir_rawdata <- '/Users/junnawang/YaleLab/data_server/'
####################End Attention

site_info <- read.csv(file.path('data', 'site_info.csv'))
feature_gs <- read.csv(file.path('data', 'growing_season_feature_EuropFlux.csv'))
feature_gs_AmeriFlux <- read.csv(file.path('data', 'growing_season_feature_AmeriFlux.csv'))
feature_gs <- rbind(feature_gs, feature_gs_AmeriFlux)

swc_ERA5 <- read.csv(file.path(dir_rawdata, "ERA5_daily_swc_1990_2024_allsites.csv"))
swc_ERA5$date <- as.Date(swc_ERA5$date)
swc_ERA5$YEAR <- year(swc_ERA5$date)
swc_ERA5$MONTH <- month(swc_ERA5$date)
swc_ERA5$DAY <- day(swc_ERA5$date)

# outcome data frame
outcome <- data.frame(site_ID = character(), RMSE = double(), R2 = double(), control_year = double(), window_size = integer(), 
                      nwindow = integer(), TAS = double(), TASp = double())
outcome_siteyear <- data.frame()

# priors of ER models:
priors_temp <- brms::prior("normal(2, 5)", nlpar = "C0", lb = 0, ub = 10) +
  brms::prior("normal(0.1, 1)", nlpar = "alpha", lb = 0, ub = 0.2) + 
  brms::prior("normal(-0.001, 0.1)", nlpar = "beta", lb = -0.01, ub = 0.0)

priors_water <- brms::prior("normal(10, 10)", nlpar = "Hs", lb = 0, ub = 1000)

priors_gpp <- brms::prior("normal(0.5, 2)", nlpar = "k2", lb = 0, ub = 10)

# Below are sites whose inter-annual TS was not strongly correlated to inter-annual TA, so use TA to obtain TS by linear regression. 
site_TS_issue <- c("BE-Bra", "CA-Cbo", "CA-Gro", "CA-Mer", "CA-Obs", "CA-TP3", "CH-Lae", "DE-RuC", "DE-SfS", "FI-Sod",
                   "GF-Guy", "IT-Ren", "NL-Loo", "US-Bar", "US-BZB", "US-BZF", "US-BZS", "US-CMW", "US-GLE", "US-Ha2",
                   "US-IB2", "US-Jo2", "US-KL2", "US-Kon", "US-LL1", "US-MBP", "US-Myb", "US-NC4", "US-Tw1", "US-ICt",
                   "BE-Dor", "CA-TP4", "UK-AMo", "Ru-Fyo", "ZA-Kru", "IT-Tor")
  
for (id in 1:nrow(site_info)) {
  # id = 27  # 1, 6, 77, 89, 64, 1:nrow(site_info)
  print(id)
  name_site <- site_info$site_ID[id]
  print(name_site)
  
  #-------------------------------------------DATA PREPARATION--------------------------
  # read data
  path <- file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.csv'))
  if (file.exists(path)) {
    a_measure_night_complete <- read.csv(path)
    ac <- read.csv(file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.csv')))
  }
  if (site_info$SWC_use[id] == 'YES') {
    a_measure_night_complete <- a_measure_night_complete %>% filter(!is.na(SWC))
  }
  
  # if no measured SWC data use daily SWC from ERA5 land
  ####################################################
  if (site_info$SWC_use[id] == 'NO') {
    a_measure_night_complete$SWC <- NULL
    ac$SWC <- NULL
    # use SWC data from ERA5 land climate reanalysis
    swc_ERA5_site <- swc_ERA5[swc_ERA5$name == name_site, ]
    # attach to a_measure_night_complete and ac
    a_measure_night_complete <- a_measure_night_complete %>% left_join(swc_ERA5_site[, c('YEAR', 'MONTH', 'DAY', 'SWC')], by = c('YEAR', 'MONTH', 'DAY'))
    ac <- ac %>% left_join(swc_ERA5_site[, c('YEAR', 'MONTH', 'DAY', 'SWC')], by = c('YEAR', 'MONTH', 'DAY'))
    site_info$SWC_use[id] = 'YES'
  }
  ###########################################

  # For the 36 sites, obtain TS from simple linear regression of TA, in order to follow reviewer's suggestion
  if (name_site %in% site_TS_issue) {
    mod_lm <- lm(data = a_measure_night_complete, TS ~ TA, na.action = na.omit)
    TS_pred <- predict(mod_lm, newdata = data.frame(TA = a_measure_night_complete$TA), na.action = na.pass)
    a_measure_night_complete$TS[!is.na(TS_pred)] <- TS_pred[!is.na(TS_pred)]
    TS_pred <- predict(mod_lm, newdata = data.frame(TA = ac$TA), na.action = na.pass)
    ac$TS[!is.na(TS_pred)] <- TS_pred[!is.na(TS_pred)]
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
  
  # use non-overlapping windows and determine number of windows for growing season
  nwindow <- max(round((gEnd - gStart + 1) / window_size), 1)

  #------------------------------------PREPARE FOR formula, stprm, priors of ER models---------------------------------
  if (site_info$SWC_use[id] == 'YES') {
    frmu <- NEE ~ exp(alpha * TS + beta*TS^2) * SWC / (Hs + SWC) * (C0 + NEE_daytime * k2)
    param <- list(alpha ~ 0 + year_group,
                  beta ~ 0 + year_group,
                  C0 ~ 0 + year_group,
                  Hs+k2 ~ 1)
    priors <- priors_temp + priors_water + priors_gpp
  } else {
    frmu <- NEE ~ exp(alpha * TS + beta*TS^2) * (C0 + NEE_daytime * k2)
    param <- list(alpha ~ 0 + year_group,
                  beta ~ 0 + year_group,
                  C0 ~ 0 + year_group,
                  k2 ~ 1)
    priors <- priors_temp + priors_gpp
  }
  
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
    
    # get reference temperature, SWC, and NEEday of each window
    TSref <- mean(ac$TS[between(ac$DOY, window_start, window_end)], na.rm=T)
    NEEdayref <- mean(ac_day$NEE_daytime[between(ac_day$DOY, window_start, window_end)], na.rm=T)
    if (site_info$SWC_use[id] == 'YES') {
      SWCref <- mean(ac$SWC[between(ac$DOY, window_start, window_end)], na.rm=T)
      data_ref <- data.frame(TS=TSref, NEE_daytime=NEEdayref, SWC=SWCref)
    } else {
      data_ref <- data.frame(TS=TSref, NEE_daytime=NEEdayref)
    }
    
    data_window <- data.frame()
    # collect data for each year
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
      # if the two rules are violated, extend window size. 
      extend_days <- 0
      while(nrow(data_subset) < nobs_threshold | !between(TSref, quantile(data_subset$TS, 0.025, na.rm=T), quantile(data_subset$TS, 0.975, na.rm=T))) {
        extend_days <- extend_days + 3
        data_subset <- a_measure_night_complete %>% filter(growing_year == iyear) %>% 
          filter(between(DOY, window_start - extend_days, window_end + extend_days))
        # some conditions to break out to avoid dead loop
        if (nrow(data_subset) >= nobs_threshold & (window_end + extend_days) >= gEnd & TSref >= max(data_subset$TS, na.rm=T)) {break}
        if (extend_days >= 24) { break }   # max window size: two months
      }
      df_site_year_window[icount, 4] <- nrow(data_subset)
      df_site_year_window[icount, 5] <- extend_days
      
      # ensure enough observations
      if (nrow(data_subset) <= 25) { next }   # this is mainly for US-Uaf, which has few data
      
      # ensure enough temperature range
      if (!between(TSref, quantile(data_subset$TS, 0.025, na.rm=T), quantile(data_subset$TS, 0.975, na.rm=T))) { next }
      
      data_window <- bind_rows(data_window, data_subset)
    }
    data_window$year_group <- as.factor(data_window$growing_year)
    
    # run a simple model first; if failed, adjust parameter values
    mod <- try(brms::brm(brms::bf(frmu, param, nl = TRUE),
                         prior = priors, data = data_window, iter = 2000, cores =4, chains = 4, backend = "cmdstanr", 
                         control = list(adapt_delta = 0.90, max_treedepth = 15), refresh = 0)) # , silent = 2
    
    if (!inherits(mod, "try-error")) {
      # cmdfit <- mod$fit
      # diag <- cmdfit$sampler_diagnostics()
      # n_divergent <- sum(diag[, , "divergent__",])
      np <- nuts_params(mod)
      n_divergent <- sum(subset(np, Parameter == "divergent__")$Value)
      failed_brm <- n_divergent > 0 
    } else {
      failed_brm <- TRUE
    }
    
    # if brm models fail or have divergent transitions, try another time
    if (failed_brm) {
      mod <- try(brms::brm(brms::bf(frmu, param, nl = TRUE),
                           prior = priors, data = data_window, iter = 4000, cores =4, chains = 4, backend = "cmdstanr", 
                           control = list(adapt_delta = 0.98, max_treedepth = 15), refresh = 0)) # , silent = 2
    }    
    
    # extract model results
    data_window$NEE_pred <- fitted(mod)[, "Estimate"]
    ER_obs_pred <- rbind(ER_obs_pred, data_window[between(data_window$DOY, window_start, window_end), ])
    
    # print(plot(data_window$TS, data_window$NEE, main = paste(name_site, iwindow, iyear, sep = '_')))
    # lines(data_window$TS, fitted(mod)[, "Estimate"])  
    param_names_comb <- names(brms::fixef(mod)[, "Estimate"])
    index <- grepl("_year_group", param_names_comb)
    parameters <- data.frame(growing_year = as.integer(sub(".*_year_group", "", param_names_comb[index])), 
                     param_name = sub("_year_group.*", "", param_names_comb[index]),
                     param_value = brms::fixef(mod)[index, "Estimate"])
    #
    parameters <- pivot_wider(parameters, names_from = param_name, values_from = param_value)
     
    index <- which(grepl("_Intercept$", param_names_comb))
    for (icol in index) {
      parameters[[sub('_Intercept', '', param_names_comb[icol])]] <- brms::fixef(mod)[icol, "Estimate"]
    }
  
    # add TS
    parameters <- left_join(parameters, ac_yearly_window[, c("growing_year", "TS")], by = 'growing_year')
    
    # get ER at reference temperature, water, and NEEday conditions
    data_ref <- crossing(data_ref, year_group = unique(data_window$year_group))
    data_ref$ERref <- fitted(mod, newdata=data_ref)[, "Estimate"]
    #
    data_ref$year_group <- as.integer(as.character(data_ref$year_group))
    
    if (control_year %in% data_ref$year_group) {
      ERref_control <- data_ref$ERref[control_year == data_ref$year_group]
    } else {
      # if no data in a control year during this window, use average ER across years as reference ER
      ERref_control <- mean(data_ref$ERref, na.rm=T)
    }

    # add reference ER
    parameters <- left_join(parameters, data_ref[, c("year_group", "ERref")], by = c('growing_year' = "year_group"))
    
    parameters$lnRatio <- log(parameters$ERref / ERref_control)
    
    # remove extreme values due to potentially large gaps
    x <- parameters$lnRatio
    outlier <- boxplot.stats(x, coef = 3)$out
    id.remove <- match(outlier[abs(outlier) > 1.5], x)
    if (length(id.remove) > 0) {
      parameters <- parameters[-id.remove, ]
    }
    
    # update df_site_year_window
    df_site_year_window[(icount-length(years) + 1):icount, ] <- rows_update(df_site_year_window[(icount-length(years) + 1):icount, ], parameters, by = 'growing_year')   
  }    

  # end of each window
  outcome_siteyear <- rbind(outcome_siteyear, df_site_year_window)
  
  plot <- df_site_year_window %>% ggplot(aes(x=TS, y=lnRatio, col=window)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    labs(title = name_site)
  print(plot)
  
  outcome[id, "site_ID"] <- name_site
  outcome[id, c("RMSE", "R2")] <- postResample(pred = ER_obs_pred$NEE_pred, obs = ER_obs_pred$NEE)[1:2]
  outcome[id, c("control_year", "window_size", "nwindow")] <- c(control_year, window_size, nwindow)
  # if (length(unique(df_site_year_window$window)) > 1) {
  #   outcome[id, c("TAS", "TASp")] <- summary(lm(data=df_site_year_window, lnRatio ~ TS + window, na.action = na.exclude))$coefficients[2, c(1, 4)]
  # } else {
  #   outcome[id, c("TAS", "TASp")] <- summary(lm(data=df_site_year_window, lnRatio ~ TS, na.action = na.exclude))$coefficients[2, c(1, 4)]
  # }
  
  # take account of potential autocorrelation across years
  mod_ar1 <- gls(lnRatio ~ TS + window,
                 data = df_site_year_window,
                 correlation = corAR1(form = ~ growing_year | window), 
                 na.action = na.omit)
  outcome[id, c("TAS", "TASp")] <- summary(mod_ar1)$tTable["TS", c("Value", "p-value")]
}
# end of each site

write.csv(outcome, file = file.path('data', 'outcome_temp_water_gpp_new.csv'), row.names = F)
write.csv(outcome_siteyear, file = file.path('data', 'outcome_siteyear_temp_water_gpp_new.csv'), row.names = F)

