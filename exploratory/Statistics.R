library(librarian)
shelf(dplyr, tidyr, lubridate, ggpubr, ggplot2, ggrepel, ggsignif, png, grid, vip, trend)
#
rm(list=ls())

options(na.action = na.omit)

acclimation <- read.csv('data/acclimation_data.csv')
acclimation <- acclimation %>% filter(!site_ID %in% c('IT-Lav', "IT-TrF") )
acclimation$nyear <- NA

data.future <- read.csv('data/acclimation_data_future_ssp245.csv')
data.future$dTS_gs <- 0

feature_gs <- read.csv('data/growing_season_feature_EuropFlux.csv')
feature_gs_AmeriFlux <- read.csv('data/growing_season_feature_AmeriFlux.csv')
feature_gs <- rbind(feature_gs, feature_gs_AmeriFlux)

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Users/junnawang/YaleLab/data_server/'
####################End Attention

files <- list.files(file.path(dir_rawdata, 'RespirationData'), pattern = '_ac.csv$', full.names = TRUE)
##
## calculate warming trend
nyears = 0    # 1645 site years
df.warm <- data.frame(trend_ta = double(), p_ta = double(), trend_ts = double(), p_ts = double(), 
                      R2_all = double(), R2_annual = double(), p_annual = double(), 
                      iaTA = double(), iaTS = double(), slope = double(), site_ID = character())
for (i in 1:length(files)) {
  
  # if (i %in% c(7, 8,11,14,15,35,40,49,62,63,66,73,76,77,78,79,80,84,96,100,102,110)) {
  #   name_site <- substring(files[i], nchar(files[i])-12, nchar(files[i])-7)
  #   print(paste0(i, name_site))
  #   next
  # }
  # i = 76
  
  name_site <- substring(files[i], nchar(files[i])-12, nchar(files[i])-7)
  print(paste0(i, name_site))
  ac <- read.csv(files[i])
  a_measure_night_complete <- read.csv(file.path(dir_rawdata, 'RespirationData', paste0(name_site, '_nightNEE.csv')))
  good_years <- unique(a_measure_night_complete$YEAR)
  nyears = nyears + length(good_years)
  acclimation$nyear[acclimation$site_ID == name_site] <- length(good_years)
  
  # plot(ac$YEAR+ac$DOY/366, ac$TA)
  # plot(ac$YEAR+ac$DOY/366, ac$TS)
  
  gStart <- feature_gs$gStart[feature_gs$site_ID == name_site]
  gEnd <- feature_gs$gEnd[feature_gs$site_ID == name_site]
  
  Ta <- a_measure_night_complete %>% filter(YEAR %in% good_years) %>% filter(between(DOY, gStart, gEnd)) %>% 
    group_by(YEAR) %>% summarise(TA = mean(TA, na.rm=T)) %>% filter(!is.na(TA))
  res <- sens.slope(Ta$TA, conf.level = 0.95)
  df.warm[i, 1] <- res$estimates
  df.warm[i, 2] <- res$p.value
  
  Ts <- a_measure_night_complete %>% filter(YEAR %in% good_years) %>% filter(between(DOY, gStart, gEnd)) %>% 
    group_by(YEAR) %>% summarise(TS = mean(TS, na.rm=T)) %>% filter(!is.na(TS))
  res <- sens.slope(Ts$TS, conf.level = 0.95)
  df.warm[i, 3] <- res$estimates
  df.warm[i, 4] <- res$p.value
  
  df.warm[i, 5] <- summary(lm(data = a_measure_night_complete, TS ~ TA))$r.squared
  
  TT <- a_measure_night_complete %>% filter(YEAR %in% good_years) %>% filter(between(DOY, gStart, gEnd)) %>% 
    group_by(YEAR) %>% summarise(TS = mean(TS, na.rm=T), TA = mean(TA, na.rm=T)) %>% filter(!is.na(TS) & !is.na(TA))
  
  df.warm[i, 6] <- summary(lm(data = TT, TS ~ TA))$r.squared
  df.warm[i, 7] <- summary(lm(data = TT, TS ~ TA))$coefficients[2, 4]
  
  #
  df.warm[i, 8] <- sd(TT$TA)
  df.warm[i, 9] <- sd(TT$TS)
  
  df.warm[i, 10] <- summary(lm(data = a_measure_night_complete, TS ~ TA))$coefficients[2, 1]
  
  df.warm[i, 11] <- name_site
}
##
warming_sites <- df.warm$site_ID[df.warm$p_ta < 0.05 & df.warm$trend_ta > 0]
acclimation_warm <- acclimation %>% filter(site_ID %in% warming_sites)
acclimation_compensate <- acclimation %>% filter(TAS_tot < 0 & TAS_totp < 0.05)
acclimation_enhance <- acclimation %>% filter(TAS_tot > 0 & TAS_totp < 0.05)


model <- lm(data=df.warm, trend_ts ~ trend_ta)
summary(model)
residuals(model)
# 7   8  11  14  15  35  40  49  62  63  66  73  76  77  78  79  80  84  96 100 102 107 108 110
# we have ensured R2 of hourly/subhourly TS and TA data > 0.5; interannual growing-season TS and TS are significant correlated (p < 0.01); 

acclimation <- read.csv('data/outcome_temp.csv')
acclimation.direct <- read.csv('data/outcome_temp_water_gpp.csv')
mean(acclimation$R2)
sd(acclimation$R2)

mean(acclimation.direct$R2)
sd(acclimation.direct$R2)

outcome_TA <- read.csv('exploratory/result_SI/outcome_temp_TA.csv')
sum(outcome_TA$TASp < 0.05)
sum(outcome_TA$TASp < 0.05 )


