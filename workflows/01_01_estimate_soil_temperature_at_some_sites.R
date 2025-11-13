# title: Estimate Soil Temperature Based on Air Temp and net Radiation using Random Forests or linear regression for 16 sites

#  Description
#  this script intends to estimate surface soil temperature at sites with lots of missing soil temperature data
#  the method is random forests for 8 sites ("DE-Akm", "FR-Pue", "US-Los", "US-SRG", "CA-Man", "US-Ced", "US-Ho1", "CZ-RAJ"), 
#  but for the remaining sites, we used linear regression due to limited (1-2 years) TS observations or not enough observation of net radiation.
#   
#  Author: Junna Wang, July 2025
#  It will take ~40 min to finish running this script. 

library(librarian)
shelf(randomForest, caret, amerifluxr, lubridate, tidyverse, ranger)

rm(list=ls())

dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'

#--------A function to predict soil temperature
# inputs: a data frame data with 2-3 columns (TS, TA, NETRAD), the last column is optional
# output: the same data frame data, but with one added column TS_pred
predict_soil_temp <- function(data, use_NETRAD) {

  # separate data into trained or tested
  set.seed(222)
  # use a maximum of 60000 data to train and test the model; too many data will cause RF super slow and may not improve accuracy. 
  sampled <- sample(1:nrow(data), size = min(60000, nrow(data)), replace = FALSE)
  ind <- sample(2, length(sampled), replace = TRUE, prob = c(0.7, 0.3))
  train <- data[sampled[ind==1],]        
  test  <- data[sampled[ind==2],]
  # remove NA data
  train <- na.omit(train)
  test  <- na.omit(test)
  #
  if (use_NETRAD) {
    rf <- randomForest(formula = TS ~ TA + NETRAD, data=train)  # this takes lots of time
    lm <- lm(data=data, TS ~ TA + NETRAD)
  } else {
    rf <- randomForest(formula = TS ~ TA, data=train)  # this takes lots of time
    lm <- lm(data=data, TS ~ TA)
  }
  
  y0 <- predict(rf, train)
  y1 <- predict(rf, test)
  
  print('Training data performance by random forests: ')
  print(postResample(pred=y0, obs=train$TS)) 
  
  print('Testing data performance by random forests: ')
  print(postResample(pred=y1, obs=test$TS))

  # compare with linear regression
  print('Performance by linear regression: ')
  print(summary(lm))
  # random forest is much better than lm;
  
  # do the prediction
  data$TS_pred <- predict(rf, data)  
  
  return(data[, c("TIMESTAMP", "TS_pred")])
}

#-------------Predict soil temperature for 9 Ameriflux sites using the function above
site_info <- read.csv('data/site_info.csv')

files_AmeriFlux_BASE <- list.files(file.path(dir_rawdata, "SiteData", "AmeriFlux_BASE"), pattern=".zip$", full.names = T)
#
files_FLUXNET2015  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2015', "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)_", full.names = TRUE)
files_FLUXNET2020  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2020', "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)_", full.names = TRUE)
files_ICOS_after2020 <- list.files(path=file.path(dir_rawdata, 'SiteData', 'ICOS_after2020'), pattern = ".csv$", full.names = TRUE)
#
files_FLUXNET2025 <- list.files(file.path(dir_rawdata, "SiteData", "FLUXNET07202025", "unzip"), pattern=".csv$", full.names = T)
files_ICOS2025 <- list.files(file.path(dir_rawdata, "SiteData", "Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1", "unzip"), pattern=".csv$", full.names = T)

# 
id_estimate_TS <- which(site_info$estimate_Ts == 'YES')


for (id in id_estimate_TS) {
  # id <- id_estimate_TS[4]
  name_site <- site_info$site_ID[id]
  
  if (site_info$source[id] %in% c("AmeriFlux_BASE", "AmeriFlux_FLUXNET")) {
    a <- amf_read_base(files_AmeriFlux_BASE[grepl(name_site, files_AmeriFlux_BASE)], parse_timestamp=TRUE, unzip = T)
    data <- data.frame(TIMESTAMP = a$TIMESTAMP, YEAR = a$YEAR, DOY = a$DOY, HOUR = a$HOUR, MINUTE = a$MINUTE, 
                       TS = a[, trimws(site_info$TS[id])], TA = a[, trimws(site_info$TA[id])])
    if (name_site %in% c('US-Los', "US-Ced")) {
      data$NETRAD <- a$NETRAD_1_1_1
    } else if (name_site %in% c("US-SRG", "CA-Man")) {
      data$NETRAD <- a$NETRAD
    } else if (name_site %in% c("US-Ho1")) {
      data$NETRAD <- a$NETRAD_2_1_1
    } 
  } else {
    data_source <- site_info$source[id]
    data_source <- unlist(strsplit(data_source, "_"))
    for (isource in 1:length(data_source)) {
      if (data_source[isource] == "FLUXNET2015") {
        files <- files_FLUXNET2015
      } else if (data_source[isource] == "FLUXNET2020") {
        files <- files_FLUXNET2020
      } else if (data_source[isource] == "ICOS2020") {
        files <- files_ICOS_after2020
      } else if (data_source[isource] == "FLUXNET2025") {
        files <- files_FLUXNET2025
      } else if (data_source[isource] == "ICOS2025") {
        files <- files_ICOS2025
      }
      file  <- files[grepl(name_site, files)]
      a_tmp <- read.csv(file)
      if (isource == 1) {
        a <- a_tmp
      } else {
        # need to check if they can connect correctly, or use bind_rows. 
        irow_start <- which(a_tmp$TIMESTAMP_START == a$TIMESTAMP_START[nrow(a)]) + 1
        if (is.numeric(irow_start) && length(irow_start) == 0) {
          irow_start = 1
        }
        a <- bind_rows(a, a_tmp[irow_start:nrow(a_tmp), ])
      }
    }
    rm(a_tmp)
    # -9999 means missing information.
    a[a==-9999] <- NA    

    dt <- ymd_hm(a$TIMESTAMP_START[2]) - ymd_hm(a$TIMESTAMP_START[1])
    a$TIMESTAMP  <- ymd_hm(a$TIMESTAMP_START) + dt / 2
    a$DOY <- yday(a$TIMESTAMP)
    data <- data.frame(TIMESTAMP = a$TIMESTAMP, YEAR = year(a$TIMESTAMP), DOY = a$DOY, 
                       HOUR = hour(a$TIMESTAMP), MINUTE = minute(a$TIMESTAMP), TS = a$TS_F_MDS_1, TA = a$TA_F_MDS, NETRAD=a$NETRAD)
    if (name_site == 'DE-Hte') {
      # too few data in TS_F_MDS_1, so use TS_F_MDS_2
      data$TS <- a$TS_F_MDS_2
    } else if (name_site == 'FR-Bil') {
      # abnormal Soil temperature data before 2021
      data$TS[data$YEAR < 2021] <- NA
    } else if (name_site == "FR-Pue") {
      data$TS[data$YEAR < 2016] <- NA
    }
  }

  print(id)
  plot(data$TA)
  plot(data$TS)

  # gap fill TA use DOY method if needed. 
  if (!"TA_gf" %in% colnames(data)) {
    TA_gf <- data %>% group_by(DOY, HOUR, MINUTE) %>% summarise(TA_gf=mean(TA, na.rm=T))
    data   <- data %>% left_join(TA_gf, by=c("DOY", "HOUR", "MINUTE"))
  }
  if (sum(is.na(data$TA)) > 0) {
    data$TA[is.na(data$TA)] <- data$TA_gf[is.na(data$TA)]
  }

  plot(data$TA)
  
  # gap fill NETRAD if using NETRAD method
  if ("NETRAD" %in% names(data)) {
    if (!"NETRAD_gf" %in% colnames(data)) {
      gf <- data %>% group_by(DOY, HOUR, MINUTE) %>% summarise(NETRAD_gf=mean(NETRAD, na.rm=T))
      data   <- data %>% left_join(gf, by=c("DOY", "HOUR", "MINUTE"))
    }
    if (sum(is.na(data$NETRAD)) > 0) {
      data$NETRAD[is.na(data$NETRAD)] <- data$NETRAD_gf[is.na(data$NETRAD)]
    }
    plot(data$NETRAD)    
  }
  
  if (name_site %in% c("DE-Akm", "FR-Pue", "US-Los", "US-SRG", "CA-Man", "US-Ced", "US-Ho1", "CZ-RAJ")) {
    # use NETRAD because it is available
    use_NETRAD = TRUE
    data <- predict_soil_temp(data, use_NETRAD)
  } else if (name_site %in% c("DE-Hte", "FR-FBn")) {
    # use simple linear regression because these sites have 1-year or two-year TS, and NETRAD is incomplete
    lm <- lm(data=data, TS ~ TA + NETRAD)
    data$TS_pred <- predict(lm, data)
  } else {
    # no netrad at this site
    lm <- lm(data=data[data$TA>0, ], TS ~ TA)
    data$TS_pred <- predict(lm, data)
  } 
  plot(data$TS_pred)
  
  write.csv(data, file=file.path(dir_rawdata, 'TS_RandomForest', paste0(name_site, '_TS_rfp.csv')), row.names = F)
}

# all R2 should be higher than 0.83. 
# check TS prediction quality
# data <- read.csv("/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation/TS_RandomForest/US-Ho2_TS_rfp.csv")
# plot(ymd_hms(data$TIMESTAMP[140000:170000]), data$TS_pred[140000:170000])
# data %>% group_by(year(TIMESTAMP)) %>% summarise(TS=mean(TS_pred, na.rm=T))
# 