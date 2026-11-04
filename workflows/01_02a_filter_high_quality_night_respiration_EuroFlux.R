# Purpose of this script: prepare data for fitting temperature-respiration curves for EuroFlux sites. 
# Output of this script: two .csv files. One for filtered night NEE, and the other for gap-filled whole time-series.
# Author: Junna Wang
#  It will take ~ 30 mins to finish running this script. 

# Special attention: 
# site FR-Pue - soil water provided by PI
# sites in southern hemisphere ('AU-Tum', 'ZA-Kru') have growing seasons spanning two years. 
# we added 11 European sites in Aug. 2025

library(librarian)
shelf(dplyr, lubridate)
rm(list=ls())

source(file.path("workflows", "00_functions.R"))
####################Attention: change this directory based on your own directory of raw data
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
dir_rawdata <- '/Volumes/WZZ_disk/Thermal_Acclimation'
####################End Attention
shuttle_dir = file.path(dir_rawdata, 'SiteData', 'FLUXNET_Shuttle09102026', "unzip")

files_FLUXNET2015  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2015', "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)_", full.names = TRUE)
files_FLUXNET2020  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2020', "unzip"), pattern = "_FLUXNET2015_FULLSET_(HH|HR)_", full.names = TRUE)
#
files_FLUXNET2025 <- list.files(file.path(dir_rawdata, "SiteData", "FLUXNET07202025", "unzip"), pattern=".csv$", full.names = T)
files_ICOS2025 <- list.files(file.path(dir_rawdata, "SiteData", "Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1", "unzip"), pattern=".csv$", full.names = T)

site_info <- read.csv(file.path('data', 'site_info.csv'))

feature_gs <- data.frame(site_ID=character(), gStart=double(), gEnd=double(), taStart=double(), taEnd=double(), tsStart=double(), tsEnd=double(), nyear=integer())  # growing season feature

# put a for loop here 
for (id in 1:nrow(site_info)) {
  # id = 95
  print(id)
  data_source <- site_info$source[id]
  # This script only works for FLUXNET products
  if (data_source %in% c("AmeriFlux_BASE", "FLUXNET_Shuttle", "self_compile")) {
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

  # deal with soil temperature data
  if (name_site %in% no_TS_sites) {
    ts_result <- reconstruct_missing_TS(a = a, name_site = name_site, ts_config = ts_config, long_mean_days = 8, shuttle_dir = shuttle_dir)
    a$TS_F_MDS_1 <- ts_result$data$TS_F_MDS_1
    
  } else if (name_site %in% ts_config$site) {
    ts_result <- reconstruct_TS(a=a, name_site=name_site, ts_config=ts_config, long_mean_days=8)
    a$TS_F_MDS_1 <- ts_result$data$TS_F_MDS_1
    # summary(ts_result$model)
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
  
  #
  NEE_yearly  <- a %>% group_by(DOY) %>% summarise(NEE=mean(NEE_VUT_REF, na.rm=T), TS=mean(TS_F_MDS_1, na.rm=T), TA=mean(TA_F_MDS, na.rm=T))
  NEE_yearly_night  <- a %>% filter(NIGHT==1) %>% group_by(DOY) %>% summarise(NEE=mean(NEE_VUT_REF, na.rm=T), TS=mean(TS_F_MDS_1, na.rm=T), TA=mean(TA_F_MDS, na.rm=T))
  NEE_yearly_day  <- a %>% filter(NIGHT==0) %>% group_by(DOY) %>% summarise(NEE=mean(NEE_VUT_REF, na.rm=T), TS=mean(TS_F_MDS_1, na.rm=T), TA=mean(TA_F_MDS, na.rm=T))
  
  # Estimate start and end date of the ER active period. 
  # Tropical sites, all year growing season; check if we need to add more sites. 
  if (abs(site_info$LAT[id]) < 23.5 || (site_info$Climate_class[id] %in% c("Af", "Am", "Aw", "As", "Cfa", "Cfb", "Cwa") && site_info$IGBP[id] == "EBF")) {
    gStart <- min(NEE_yearly_night$DOY, na.rm = T)
    gEnd <- max(NEE_yearly_night$DOY, na.rm = T)
  } else {
    df_ER_period <- estimate_active_ER_period(a = a, site_id = name_site, ta_var = "TA_F_MDS", ts_var = "TS_F_MDS_1", 
                                              southern_hemisphere = southern_hemisphere, threshold_fraction = 0.20,
                                              temp_threshold = 2, min_observed_doy = 180L)
    if (df_ER_period$status == "ok") {
      gStart <- df_ER_period$gStart
      gEnd   <- df_ER_period$gEnd
      if (df_ER_period$trimmed_days > 0) {
        print(paste0(name_site, " has trimmed days at boundaries ", df_ER_period$trimmed_days))
      }
    } else if (df_ER_period$status == "skipped") {
      # jump out of the loop
      next
    }
  }
  # if user inputs the start and end dates, using them.  
  if (!is.na(site_info$gStart[id])) {
    gStart <- as.numeric(site_info$gStart[id])
  }
  if (!is.na(site_info$gEnd[id])) {
    gEnd <- as.numeric(site_info$gEnd[id])
  }
  
  
  # save TS~NEE curve
  png(paste0("graphs/", name_site, "TS_NEEnight_relation.png"), width = 1200, height = 1200)
  par(mfrow = c(3, 2))
  plot(a_measure_night_complete$TS_F_MDS_1, a_measure_night_complete$NEE_VUT_REF, main=name_site)
  plot(NEE_yearly$DOY, NEE_yearly$NEE)
  mark_period(NEE_yearly$DOY, NEE_yearly$NEE, gStart, gEnd)
  
  plot(NEE_yearly$DOY, NEE_yearly$TS)
  mark_period(NEE_yearly$DOY, NEE_yearly$TS, gStart, gEnd)
  
  plot(NEE_yearly_night$DOY, NEE_yearly_night$NEE)
  mark_period(NEE_yearly_night$DOY, NEE_yearly_night$NEE, gStart, gEnd)
  
  plot(NEE_yearly_day$DOY, NEE_yearly_day$TS)
  mark_period(NEE_yearly_day$DOY, NEE_yearly_day$TS, gStart, gEnd)
  
  plot(NEE_yearly_day$DOY, NEE_yearly_day$NEE)
  mark_period(NEE_yearly_day$DOY, NEE_yearly_day$NEE, gStart, gEnd)
  dev.off()
  
  # temperature range of growing season
  tmp <- NEE_yearly %>% filter(between(DOY, gStart, gEnd))
  tsStart <- quantile(tmp$TS, c(0.025), na.rm=T)
  tsEnd   <- quantile(tmp$TS, c(0.975), na.rm=T)
  taStart <- quantile(tmp$TA, c(0.025), na.rm=T)
  taEnd   <- quantile(tmp$TA, c(0.975), na.rm=T)
  
  # # max gap and total large gap threshold
  # if (name_site %in% c('ZA-Kru', 'FI-Sod', 'GF-Guy')) {
  #   # use larger gap threshold because lack of data at these Tropical and Tundra sites. 
  #   gap_max_thresh <- 60
  #   gap_total_thresh <- 0.7   
  # } else {
  #   gap_max_thresh <- max(31, (gEnd-gStart+1) * 0.225)
  #   gap_total_thresh <- max(1/3, gap_max_thresh / (gEnd-gStart+1))
  # }

  gap_max_thresh <- (gEnd-gStart+1) * 0.5
  gap_total_thresh <- 0.7
  
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
  
  # only keep years in a_measure_night_complete
  ac <- ac %>% filter(YEAR %in% unique(a_measure_night_complete$YEAR))
  
  # save the ac and a_measure_night_complete data:
  write.csv(ac, file=file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.csv')), row.names = F)
  write.csv(a_measure_night_complete, file=file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.csv')), row.names = F)
  
  feature_gs <- bind_rows(feature_gs, data.frame(site_ID=name_site, gStart=gStart, gEnd=gEnd, taStart=taStart, taEnd=taEnd, tsStart=tsStart, tsEnd=tsEnd, nyear=length(good_years)))
}

# output growing season features
write.csv(feature_gs, file=file.path('data', 'growing_season_feature_EuropFlux.csv'), row.names = F)

# Sys.time()