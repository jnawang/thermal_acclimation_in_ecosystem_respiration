# I create an artificial time series using pre-defined parameter values, and test if methods I used can recover these parameters.
# example site: BE-Bra: mixed forest; US-Ho1: ENF; DE-Tha ENF; DE-Hai: DBF; 

# create an artificial time series
library(librarian)
shelf(dplyr, lubridate, gslnls, caret, performance, ggpubr, ggplot2, zoo, brms)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
####################End Attention

name_site <- 'BE-Bra'
path <- file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.RDS'))
if (file.exists(path)) {
  a_measure_night_complete <- readRDS(path)
  ac <- readRDS(file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.RDS')))
} else {
  next
}

  
dt = 30  # minute
if (ac$MINUTE[2] - ac$MINUTE[1] != 30) { dt = 60 }   # minutes 
window_sizes <- c(2, 3, 5, 7, 10, 14, 20)

# add rolling average GPP
if ("GPP_DT" %in% names(ac)) {
  ac_day <- ac %>% filter(daytime) %>% group_by(YEAR, DOY) %>% 
    summarise(GPP_DT1 = sum(GPP_DT, na.rm=T) * dt * 60 * 12 / 1000000, .groups = "drop") 
  
  ac_day[paste0("GPP_DT", window_sizes)] <- lapply(
    window_sizes,
    function(w) rollapply(ac_day$GPP_DT1, width = w, FUN = mean,
                          align = "right", fill = ac_day$GPP_DT1[1])
  )
}

# attach it to nighttime data
if (!"GPP_DT1" %in% names(a_measure_night_complete)) {
  a_measure_night_complete <- a_measure_night_complete %>% mutate(DOY_gpp = case_when(HOUR >= 12 ~ DOY,
                                                                                      HOUR < 12 ~ DOY - 1)) %>%
    left_join(ac_day, by = c("YEAR", "DOY_gpp"="DOY")) %>% select(-"DOY_gpp")
}

window_sizes <- c(1, window_sizes)
if ("GPP_DT" %in% names(ac)) {
  idx <- which(!is.na(a_measure_night_complete$GPP_DT1))[1]
  a_measure_night_complete[1:idx, paste0("GPP_DT", window_sizes)] <- a_measure_night_complete$GPP_DT1[idx]
}


# create an artificial time series
a_measure_night_complete$NEE_art <- 0.1 * a_measure_night_complete$GPP_DT2 * exp(a_measure_night_complete$TS * 0.03) + 
  5 * exp(a_measure_night_complete$TS * 0.08 - a_measure_night_complete$TS^2 * 0.001) * a_measure_night_complete$SWC / (a_measure_night_complete$SWC + 10) + 
  rnorm(nrow(a_measure_night_complete), mean = 0, sd = 1)
# I also need to add noise into it. 

frmu <- NEE_art ~ exp(alpha * TS + beta*TS^2) * C0 * SWC / (Hs + SWC) + k2 * GPP_DT2 * exp(gamma * TS)
stprm <- c(C0=1, alpha=-1, beta=-0.01, Hs=3, k2=0.5, gamma = 0.1)

# Under which cases, i can recover parameter values?
# time window? number of data? 

iStart = a_measure_night_complete$YEAR[1]
iEnd   = a_measure_night_complete$YEAR[nrow(a_measure_night_complete)]
gStart = 150
gEnd   = 170
for (iyear in iStart:iEnd) {
  data <- a_measure_night_complete %>% filter(YEAR==iyear) %>% filter(between(DOY, gStart, gEnd))
  
  data <- data %>% filter(!is.na(SWC))
    
  if (nrow(data) <= 10) { next }
  print(paste0(iyear, "--", nrow(data)))

  # mod2 <- gsl_nls(fn=frmu, data=data, start=stprm)
  # print(summary(mod2))
  # data$NEE_mod2 <- predict(mod2, newdata=data)
  # print(postResample(pred=data$NEE_mod2, obs=data$NEE_art)[1:2])
  
  # use brm
  priors.er <- brms::prior("normal(0.05, 0.5)", nlpar = "alpha", lb = 0, ub = 0.2) +
    brms::prior("normal(-0.001, 0.1)", nlpar = "beta", lb = -0.01, ub = 0.01) +
    brms::prior("normal(2, 10)", nlpar = "C0", lb = 0, ub = 20) +
    brms::prior("normal(10, 10)", nlpar = "Hs", lb = 0, ub = 1000) +
    brms::prior("normal(0.2, 1)", nlpar = "k2", lb = 0, ub = 1) + 
    brms::prior("normal(0.05, 0.5)", nlpar = "gamma", lb = 0, ub = 0.2)
  model.brms <- brms::brm(brms::bf(frmu, alpha+beta+C0+Hs+k2+gamma ~ 1, nl = TRUE),
                          prior = priors.er, data = data, iter = 1000, cores = 4, chains = 4, backend = "cmdstanr", 
                          control = list(adapt_delta = 0.9))
  
  print(summary(model.brms), digits = 3)
  
}

# for the gsl_nls method: it is fast, however, 
# (1) We need at least 100 data to recover these parameter values under no noise cases. 30 day window worked most times. 
# Results from 15, 20 days may not. 

# (2) If there is noise, it is so hard to reproduce the result! actually the effects of SWC is the easiest to recover, then alpha, and C0. 
# If a parameter is not significant, its recovered value is not reliable. 

# (3) even if R2 is super good, it is still hard to recover the correct parameter values if there is large noise. 

# in contrast: 
# surprisingly, brms can recover parameter values way better than the gsl_nls method, particularly when there is noise. 
# so i chose the Bayesian method. 




