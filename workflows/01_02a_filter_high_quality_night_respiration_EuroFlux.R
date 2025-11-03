# Purpose of this script: prepare for data for fitting temperature-respiration curves for EuroFlux sites. 
# Output of this script: two .RDS files. One for filtered night NEE, and the other for gap-filled whole time-series.
# Author: Junna Wang
#  It will take ~ 30 mins to finish running this script. 

# Special attention: 
# site FR-Pue - soil water provided by PI
# sites in southern hemisphere ('AU-Tum', 'ZA-Kru') have growing seasons spanning two years. 
# we added 11 European sites in Aug. 2025

library(librarian)
shelf(dplyr, lubridate)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
####################End Attention

files_FLUXNET2015  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2015', "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)_", full.names = TRUE)
files_FLUXNET2020  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2020', "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)_", full.names = TRUE)
files_ICOS_after2020 <- list.files(path=file.path(dir_rawdata, 'SiteData', 'ICOS_after2020'), pattern = ".csv$", full.names = TRUE)
#
files_FLUXNET2025 <- list.files(file.path(dir_rawdata, "SiteData", "FLUXNET07202025", "unzip"), pattern=".csv$", full.names = T)
files_ICOS2025 <- list.files(file.path(dir_rawdata, "SiteData", "Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1", "unzip"), pattern=".csv$", full.names = T)

site_info <- read.csv('data/site_info.csv')

feature_gs <- data.frame(site_ID=character(), gStart=double(), gEnd=double(), tStart=double(), tEnd=double(), nyear=integer())  # growing season feature

# put a for loop here 
for (id in 1:nrow(site_info)) {
  # id = 95
  print(id)
  data_source <- site_info$source[id]
  # This script only works for FLUXNET products
  if (data_source %in% c("AmeriFlux_BASE", "AmeriFlux_FLUXNET")) {
    next
  }
  
  name_site <- site_info$site_ID[id]
  print(name_site)
  
  # connecting multiple data sources
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
  
  # add more columns
  dt <- ymd_hm(a$TIMESTAMP_START[2]) - ymd_hm(a$TIMESTAMP_START[1])
  a$TIMESTAMP  <- ymd_hm(a$TIMESTAMP_START) + dt / 2
  a$YEAR  <- year(a$TIMESTAMP)
  a$MONTH <- month(a$TIMESTAMP)
  a$DAY   <- day(a$TIMESTAMP)
  a$DOY   <- yday(a$TIMESTAMP)
  a$HOUR  <- hour(a$TIMESTAMP)
  a$MINUTE <- minute(a$TIMESTAMP)

  # read soil temperature data, or use specific soil temperature data for some sites. 
  if (name_site %in% c('CZ-Stn')) {
    a$TS_F_MDS_1 <- a$TS_F_MDS_2  # use TS of second layer because the first layer is incomplete
    a$TS_F_MDS_1_QC <- a$TS_F_MDS_2_QC
  } else if (name_site %in% c("FR-Fon", "CH-Dav", "DE-Akm", "DE-Hte", "FR-Bil", "FR-Pue", "FR-FBn", "CZ-RAJ")) {
    df_TS <- readRDS(file=file.path(dir_rawdata, 'TS_RandomForest', paste0(name_site, '_TS_rfp.RDS')))
    df_TS <- left_join(data.frame(TIMESTAMP=a$TIMESTAMP), df_TS, by = "TIMESTAMP")
    a$TS_F_MDS_1 <- df_TS$TS_pred
    a$TS_F_MDS_1_QC[is.na(a$TS_F_MDS_1_QC) | a$TS_F_MDS_1_QC == 3] <- 2
    rm(df_TS)
  } else if (name_site == "FI-Sod") {
    # the first 5 year TS measurements was off, likely because sensor changed locations. 
    # in 2001, measurements in soil profile 2 are correct, so building relations between profiles 1 and 2
    mod1 <- lm(data = a[1:24383,], TS_F_MDS_2 ~ TS_F_MDS_1)
    pred2 <- predict(mod1, data.frame(TS_F_MDS_1 = a$TS_F_MDS_1[a$YEAR <= 2005]))
    mod2 <- lm(data = a[90000:245000,], TS_F_MDS_1 ~ TS_F_MDS_2)
    a$TS_F_MDS_1[a$YEAR <= 2005] <- predict(mod2, data.frame(TS_F_MDS_2 = pred2))
    a$TS_F_MDS_1_QC[a$YEAR <= 2005] <- 2
  }
  
  # FLUXNET data do not include soil data for Site FR-Pue, but site PI provided the SWC data
  if (name_site == 'FR-Pue') {
    swc_2003_2015 <- read.table(file.path(dir_rawdata, "SiteData", "soil_water_content_FR-Pue", "FR-Pue_SWC_Control_Corrected_2003-2015.csv"), sep=';', header=T)
    swc_2016_2020 <- read.table(file.path(dir_rawdata, "SiteData", "soil_water_content_FR-Pue", "FR-Pue_SWC_Control_Corrected_2016-2020.csv"), sep=';', header=T)
    swc_2003_2020 <- bind_rows(swc_2003_2015, swc_2016_2020)
    swc_2003_2020$TIMESTAMP_START <- dmy_hm(swc_2003_2020$date_time)
    
    a <- a %>% mutate(TIMESTAMP_START = ymd_hm(a$TIMESTAMP_START)) %>% 
    left_join(swc_2003_2020[, c("TIMESTAMP_START", "SWC_Mean")], by="TIMESTAMP_START")
    # times 100 to convert SWC into percentage
    a$SWC_Mean <- a$SWC_Mean*100
    a$SWC_F_MDS_1[!is.na(a$SWC_Mean)] <- a$SWC_Mean[!is.na(a$SWC_Mean)]
    a$SWC_F_MDS_1_QC[!is.na(a$SWC_Mean)] <- 1
  }
  
  # Southern Hemisphere, so we need to change DOY values
  southern_hemisphere = FALSE
  if (name_site %in% c('AU-Tum', 'ZA-Kru', "BR-Sa1")) {
    a$DOY[a$DOY < 183] <- a$DOY[a$DOY < 183] + 366
    southern_hemisphere = TRUE
  }
  #
  # check if SWC was missed. 
  SWC <- TRUE
  if (!'SWC_F_MDS_1' %in% colnames(a)) {
    a$SWC_F_MDS_1 <- NA
    SWC <- FALSE
  }
  if (site_info$SWC_use[id] == 'NO') {
    SWC <- FALSE
  }
  #
  # The four important variables I need to look for: TS, TA, SWC, NEE, and their corresponding variable names. 
  # plot(a$TA_F_MDS)
  # plot(a$TS_F_MDS_1)   # quality control == -9999 is missing data, 3 = poor, so filtering out these data
  # if (SWC) {
  #   plot(a$SWC_F_MDS_1)     
  # }
  # plot(a$NEE_VUT_REF)   # I should use QC=0, measured data
  
  # we decided to relax the SWC constrains
  a_measure_night_complete <- a %>% filter(!is.na(TA_F_MDS)) %>% filter(TS_F_MDS_1_QC %in% c(0, 1, 2)) %>% filter(NEE_VUT_REF_QC == 0) %>% filter(NIGHT == 1) %>% filter(NEE_VUT_REF > -5 & NEE_VUT_REF < 30) 
  # %>% filter(SWC_F_MDS_1_QC %in% c(0, 1, 2))
  
  # plot the data of each year
  yStart <- a_measure_night_complete$YEAR[1]
  yEnd   <- a_measure_night_complete$YEAR[nrow(a_measure_night_complete)]
  for (i in yStart:yEnd) {
    data <- a_measure_night_complete %>% filter(YEAR==i)
    if (nrow(data) > 100) {
      plot(data$TS_F_MDS_1, data$NEE_VUT_REF, main=i)
      plot(data$TIMESTAMP, data$NEE_VUT_REF, main=i, xlim=c(as.POSIXct(paste0(i, "-01-01 00:00:00")), as.POSIXct(paste0(i, "-12-31 00:00:00"))))
    } else {
      print(i)
    }
  }
  
  # decide which years data to exclude: rule 1: maximum data gaps during growing season < 1 month; rule2: large disturbance years due to fire, insect breakout
  # estimate growing season DOY
  NEE_yearly  <- a %>% group_by(DOY) %>% summarise(NEE=mean(NEE_VUT_REF, na.rm=T), TS=mean(TS_F_MDS_1, na.rm=T))
  NEE_yearly_night  <- a %>% filter(NIGHT==1) %>% group_by(DOY) %>% summarise(NEE=mean(NEE_VUT_REF, na.rm=T), TS=mean(TS_F_MDS_1, na.rm=T))
  
  
  # save TS~NEE curve
  png(paste0("graphs/", name_site, "TS_NEEnight_relation.png"), width = 1200, height = 1200)
  par(mfrow = c(2, 2))
  plot(a_measure_night_complete$TS_F_MDS_1, a_measure_night_complete$NEE_VUT_REF, main=name_site)
  plot(NEE_yearly$DOY, NEE_yearly$NEE)
  plot(NEE_yearly$DOY, NEE_yearly$TS)
  plot(NEE_yearly_night$DOY, NEE_yearly_night$NEE)
  dev.off()
  
  
  # I will use the mean of the first three values, and the mean of the last three values
  if (name_site %in% c('FI-Sod', 'DE-RuC')) {
    tmp <- NEE_yearly %>% filter(NEE < 0.0)
  } else {
    tmp <- NEE_yearly %>% filter(NEE < max(min(NEE_yearly$NEE, na.rm=T) * 0.2, -0.8))
  }
  
  gStart <- as.integer(mean(tmp$DOY[7])) - 4
  gEnd   <- as.integer(mean(tmp$DOY[(nrow(tmp)-6)])) + 4
  
  # In general, growing-season temperature should be > 0 C
  gStart <- max(gStart, min(which(NEE_yearly$TS >= 0)))
  if (!is.na(site_info$gStart[id])) {
    gStart <- as.numeric(site_info$gStart[id])
  }
  if (!is.na(site_info$gEnd[id])) {
    gEnd <- as.numeric(site_info$gEnd[id])
  }
  
  # temperature range of growing season
  tStart <- quantile(tmp$TS, c(0.025), na.rm=T)
  tEnd   <- quantile(tmp$TS, c(0.975), na.rm=T) 
  
  # some sites NEE measurements under TS < 2C was inaccurate, so for these sites tStart should be greater than 2C
  if (name_site %in% c('CH-Dav')) {
    tStart <- max(tStart, 2)
    a_measure_night_complete <- a_measure_night_complete %>% filter(TS_F_MDS_1 >= 2.0)
  }
  
  # max gap and total large gap threshold
  if (name_site %in% c('ZA-Kru', 'FI-Sod', 'GF-Guy')) {
    # use larger gap threshold because lack of data at these Tropical and Tundra sites. 
    gap_max_thresh <- 60
    gap_total_thresh <- 0.7   
  } else {
    gap_max_thresh <- max(31, (gEnd-gStart+1) * 0.225)
    gap_total_thresh <- max(1/3, gap_max_thresh / (gEnd-gStart+1))
  }
  
  # find long gap years
  # I need to check if this works for the southern hemisphere.
  a_check_gaps <- a_measure_night_complete[, c("TIMESTAMP", "DOY")]
  
  # adding the start and end of growing season date
  if (southern_hemisphere) {
    gStart_original <- ifelse(gStart > 366, gStart-366, gStart)
    gEnd_original <- ifelse(gEnd > 366, gEnd-366, gEnd)
    a_gs_dates <- data.frame(TIMESTAMP = c(as.POSIXct((gStart_original - 1) * 86400 + as.numeric(dt/2, units = "secs"), origin = paste0(yStart:yEnd, "-01-01"), tz = "UTC"), 
                                           as.POSIXct((gEnd_original - 1) * 86400 + as.numeric(dt/2, units = "secs"), origin = paste0(yStart:yEnd, "-01-01"), tz = "UTC")), 
                             DOY = c(rep(gStart, yEnd - yStart + 1), rep(gEnd, yEnd - yStart + 1)))
  } else {
    a_gs_dates <- data.frame(TIMESTAMP = c(as.POSIXct((gStart - 1) * 86400 + as.numeric(dt/2, units = "secs"), origin = paste0(yStart:yEnd, "-01-01"), tz = "UTC"), 
                                           as.POSIXct((gEnd - 1) * 86400 + as.numeric(dt/2, units = "secs"), origin = paste0(yStart:yEnd, "-01-01"), tz = "UTC")), 
                             DOY = c(rep(gStart, yEnd - yStart + 1), rep(gEnd, yEnd - yStart + 1)))
  }
  a_check_gaps <- a_check_gaps %>% rbind(a_gs_dates) %>% unique %>% arrange(TIMESTAMP)

  good_years <- a_check_gaps %>% 
    mutate(growing_year = case_when(DOY <= 366 ~ year(TIMESTAMP), 
                                    TRUE ~ year(TIMESTAMP) - 1)) %>% 
    arrange(growing_year, TIMESTAMP) %>%  # <-- Sort before lag
    group_by(growing_year) %>%
    filter(DOY >= gStart & DOY <= gEnd) %>% 
    mutate(lag_date = dplyr::lag(TIMESTAMP)) %>% 
    mutate(gap = as.numeric(difftime(TIMESTAMP, lag_date, units = "days"))) %>% 
    mutate(gap_large = ifelse(gap > 14, gap, 0)) %>%
    filter(!is.na(gap)) %>%
    summarise(gap_max = max(gap, na.rm=T), gap_total = sum(gap_large, na.rm=T) / (gEnd - gStart + 1)) %>% 
    filter(gap_max < gap_max_thresh &  gap_total < gap_total_thresh) %>%
    distinct(growing_year) %>% pull()
  
  # # ensure every degree of TS have measurements
  # bad_TS_years <- a_measure_night_complete %>% filter(between(TS_F_MDS_1, tStart, tEnd)) %>% 
  #   arrange(YEAR, TS_F_MDS_1) %>% group_by(YEAR) %>%
  #   mutate(gap = TS_F_MDS_1 - lag(TS_F_MDS_1)) %>%
  #   summarise(max_gap = max(gap, na.rm = TRUE)) %>%
  #   filter(max_gap > max(0.2 * (tEnd - tStart), 1.0)) %>% distinct(YEAR) %>% pull()
  # if (length(bad_TS_years) >= 1) {
  #   good_years <- setdiff(good_years, bad_TS_years)    
  # }
  
  years2remove_automation <- setdiff(yStart:yEnd, good_years)
  
  # years to remove due to PI reported disturbances and bad data quality
  year_str <- site_info$year_removed[id]
  year_parts <- strsplit(year_str, ",")[[1]]
  years2remove <- unlist(lapply(year_parts, function(part) {
    if (grepl(":", part)) {
      rng <- as.numeric(strsplit(part, ":")[[1]])
      seq(rng[1], rng[2])
    } else {
      as.numeric(part)
    }
  }))
  
  if (length(years2remove) >= 1) {
    good_years <- setdiff(good_years, years2remove)
  }
  
  
  # select good_years based on growing year
  a_measure_night_complete <- a_measure_night_complete  %>% 
    mutate(growing_year = case_when(DOY <= 366 ~ YEAR, 
                                    TRUE ~ YEAR - 1)) %>% 
    filter(growing_year %in% good_years)
  # change column names
  a_measure_night_complete$TA <- a_measure_night_complete$TA_F_MDS
  a_measure_night_complete$TS <- a_measure_night_complete$TS_F_MDS_1
  a_measure_night_complete$SWC <- a_measure_night_complete$SWC_F_MDS_1
  a_measure_night_complete$NEE <- a_measure_night_complete$NEE_VUT_REF
  
  a_measure_night_complete <- a_measure_night_complete %>% dplyr::select(c(YEAR, MONTH, DAY, DOY, HOUR, MINUTE, NEE, TA, TS, SWC))
  
  # select useful column and year to output
  ac <- a[, c("YEAR", "MONTH", "DAY", "DOY", "HOUR", "MINUTE")]
  ac$NEE <- a$NEE_VUT_REF
  ac$NEE_QC  <- a$NEE_VUT_REF_QC
  ac$TA  <- a$TA_F_MDS
  ac$TS  <- a$TS_F_MDS_1
  ac$SWC <- a$SWC_F_MDS_1
  ac$NEE_uStar_f <- a$NEE_VUT_REF
  ac$daytime <- a$NIGHT!=1
  ac$SW_IN <- a$SW_IN_F_MDS
  ac$GPP_DT <- a$GPP_DT_VUT_REF
  
  # only keep years between start year and end year
  iStart = a_measure_night_complete$YEAR[1]
  iEnd   = a_measure_night_complete$YEAR[nrow(a_measure_night_complete)]
  ac <- ac %>% filter(between(YEAR, iStart, iEnd))
  
  # save the ac and a_measure_night_complete data:
  saveRDS(ac, file=file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.RDS')))
  saveRDS(a_measure_night_complete, file=file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.RDS')))
  
  feature_gs <- bind_rows(feature_gs, data.frame(site_ID=name_site, gStart=gStart, gEnd=gEnd, tStart=max(tStart, 0.0), tEnd=tEnd, nyear=length(good_years)))
}

# output growing season features
write.csv(feature_gs, file='data/growing_season_feature_EuropFlux.csv', row.names = F)

# Sys.time()