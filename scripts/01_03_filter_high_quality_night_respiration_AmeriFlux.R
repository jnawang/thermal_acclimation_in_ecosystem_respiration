# Purpose of this script: prepare for data for fitting temperature-respiration curves for AmeriFlux sites. 
# Output of this script: two .RDS files. One for selected night NEE, and the other for gap-filled whole time-series.
# Author: Junna Wang

# step by step notes: 
# step 1: double check removed years. [some sites cannot be eliminated]
# step 1: add statement of two special sites.
# step 2: ensure the previously used sites all worked. 
# step 3: clearly indicate data sources (routine of data sources)
# step 4: ensure other sites also worked. 

# to do list
# Sites we want to remove soil temperature: TS < 2C.
# US-Uaf: first 5 year TS has problem. 
# a few stations miss RH data
# US-Los, I should only use one section of TS data?
# US-Ha1: remove data when TS < 2C
# 
# check each comment

library(librarian)
shelf(dplyr, lubridate, amerifluxr, suncalc, REddyProc, lutz, zoo)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
####################End Attention

files_AmeriFlux_BASE <- list.files(file.path(dir_rawdata, "SiteData", "AmeriFlux_BASE"), pattern=".zip$", full.names = T)

site_info <- read.csv('data/site_info.csv')

feature_gs <- data.frame(site_ID=character(), gStart=double(), gEnd=double(), tStart=double(), tEnd=double(), nyear=integer())  # growing season feature

# put a for loop here; nrow(site_info)
for (id in 1:nrow(site_info)) {
  # id = 42
  print(id)
  data_source <- site_info$source[id]
  # This script only works for FLUXNET products
  if (!data_source %in% c("AmeriFlux_BASE", "AmeriFlux_FLUXNET")) {
    next
  }
  
  name_site <- site_info$site_ID[id]
  print(name_site)
  
  # read data
  a <- amf_read_base(files_AmeriFlux_BASE[grepl(name_site, files_AmeriFlux_BASE)], parse_timestamp=TRUE, unzip = T)
  a[a==-9999] <- NA
  if (name_site == "US-Myb") {
    # use data from the second year, because lots of missing NEE in the first year. 
    a <- a[17521:245424, ]
    # combine TS data at two depth
    a$TS_2_1_1[is.na(a$TS_2_1_1)] <- a$TS_2_2_1[is.na(a$TS_2_1_1)] * 0.9324892 + 0.9077066
  }
  ######################################decide day or night for every TIMESTAMP#######################################
  # get site long and lat
  sites <- amf_site_info()
  long_site  <- sites$LOCATION_LONG[sites$SITE_ID==name_site]
  lat_site   <- sites$LOCATION_LAT[sites$SITE_ID==name_site]
  
  # get the exact time zone, not using daylight saving so choose a winter date
  tz <- tz_offset(as.Date("2000-01-01"), tz_lookup_coords(lat=lat_site, lon=long_site, method='accurate'))
  if (tz$utc_offset_h > 0) {
    tz_site  <- paste0('GMT-', tz$utc_offset_h)
  } else {
    tz_site  <- paste0('GMT+', -tz$utc_offset_h)
  }
  
  # get sunrise and sunset time
  sunrise_set <- getSunlightTimes(
    date = seq.Date(as.Date(a$TIMESTAMP[1]), as.Date(a$TIMESTAMP[nrow(a)]), by = 1),       # Date format: YYYY-MM-DD
    keep = c("sunrise", "sunset"),
    lat = lat_site,
    lon = long_site,
    tz  = tz_site
  )
  
  # Force the time zone to GMT in order to be consistent with the date format in AmeriFlux database
  sunrise_set$sunrise <- as.POSIXct(as.character(sunrise_set$sunrise), tz='GMT')
  sunrise_set$sunset <- as.POSIXct(as.character(sunrise_set$sunset), tz='GMT')
  
  # create ancillary df
  ac <- a[, 1:7]  
  ac$DATE <- as.Date(ac$TIMESTAMP) 
  ac <- ac %>% left_join(sunrise_set, by = c('DATE' = 'date'))
  
  # determine: day-time or night-time
  dt <- a$TIMESTAMP[2] - a$TIMESTAMP[1]
  a$daytime <- a$TIMESTAMP + dt/2.0 >= ac$sunrise & a$TIMESTAMP - dt/2.0 <= ac$sunset
  
  if (name_site == "US-CMW") {
    # only use data after 2000
    a <- a %>% filter (YEAR > 2000)
    ac <- ac %>% filter (YEAR > 2000)
  }
  if (name_site == "CA-Mer") {
    # only use data after 1999, the first two-year data have some problems. 
    a <- a %>% filter (YEAR > 1999)
    ac <- ac %>% filter (YEAR > 1999)
  }
  
  ######################################decide which variables to use#######################################
  # modify the dataframe ac
  if ("DATE" %in% colnames(ac)) {
    ac = subset(ac, select=-c(DATE, lat, lon, sunrise, sunset))
  }
  
  if (name_site == 'US-Ha2') {
  # combine data from two locations
    a$FC_A_1_1 <- a$FC_1_1_1
    a$FC_A_1_1[is.na(a$FC_A_1_1)] <- a$FC_2_1_1[is.na(a$FC_A_1_1)]
    a$TA_A_1_1 <- a$TA_1_1_1
    a$TA_A_1_1[is.na(a$TA_A_1_1)] <- a$TA_2_1_1[is.na(a$TA_A_1_1)]
    a$PPFD_IN_A_1_1 <- a$PPFD_IN_1_1_1
    a$PPFD_IN_A_1_1[is.na(a$PPFD_IN_A_1_1)] <- a$PPFD_IN_2_1_1[is.na(a$PPFD_IN_A_1_1)]
    a$RH_A_1_1 <- a$RH_1_1_1
    a$RH_A_1_1[is.na(a$RH_A_1_1)] <- a$RH_2_1_1[is.na(a$RH_A_1_1)]
    a$USTAR_A_1_1 <- a$USTAR_1_1_1
    a$USTAR_A_1_1[is.na(a$USTAR_A_1_1)] <- a$USTAR_2_1_1[is.na(a$USTAR_A_1_1)]
  } 
  
  NEEvars <- trimws(unlist(strsplit(site_info$NEE[id], "\\+")))
  for (iNEEvar in 1:length(NEEvars)) {
    if (iNEEvar == 1) {
      ac$NEE <- a[, NEEvars[iNEEvar]]
    } else {
      ac$NEE <- ac$NEE + a[, NEEvars[iNEEvar]]
    }
  }
  
  if (site_info$FC[id] != "") {
    # remove the NEE data without FC measurements
    ac$NEE[is.na(a[, trimws(site_info$FC[id])])]  <- NA
  }
  
  # combine NEE time series and FC time series; because either one is incomplete
  if (name_site == "CA-Man") {
    ac$NEE[is.na(ac$NEE)] <- a$FC[is.na(ac$NEE)] + a$SC[is.na(ac$NEE)]
  }
  
  # air temperature
  if (site_info$TA[id] != "") {
    ac$TA  <- a[, trimws(site_info$TA[id])]
  }
  
  # soil temperature TS
  if (site_info$estimate_Ts[id] == 'YES') {
    df_TS <- readRDS(file=file.path(dir_rawdata, 'TS_RandomForest', paste0(name_site, '_TS_rfp.RDS')))
    df_TS <- left_join(data.frame(TIMESTAMP=ac$TIMESTAMP), df_TS, by = "TIMESTAMP")
    ac$TS <- df_TS$TS_pred
    rm(df_TS)
  } else {
    if (site_info$TS[id] != "") {
      ac$TS  <- a[, trimws(site_info$TS[id])]
    } 
    # deal with special cases
    if (name_site %in% c('US-NR1', 'US-ICh', "US-ICs")) {
      # use PI gap-filled data
      ac$TS[is.na(ac$TS)] <- a$TS_PI_1[is.na(ac$TS)]
    } else if (name_site == "US-Cwt") {
      # this site has no TS measurements, so we used TS-TA relationships from nearby US-xGB of the same DBF category. 
      ac$TS  <- ac$TA * 0.64718 + 5.13873
    } else if (name_site == "US-MBP") {
      # this site only missed a few TS data, so only estimate these missing data.
      ac$TS[is.na(ac$TS)] <- ac$TA[is.na(ac$TS)] * 0.3688005 + 5.8670273
    }
  }
  
  # soil water SWC
  if (site_info$SWC_use[id] == 'YES') {
    ac$SWC <- a[, trimws(site_info$SWC[id])]
    # deal with special cases
    if (name_site == 'US-NR1') {
      # use PI gap-filled data
      ac$SWC[is.na(ac$SWC)] <- a$SWC_PI_1[is.na(ac$SWC)]
    } 
    # interpolate SWC data if necessary because TS sensor has different frequency with other data
    x <- ac$SWC
    rle_x <- rle(is.na(x))
    min_gap <- min(rle_x$lengths[rle_x$values == TRUE])
    print(paste0('minimum soil gaps: ', min_gap))
    if (min_gap <= 8) {
      # Perform interpolation on eligible gaps only
      ac$SWC <- na.approx(x, na.rm = FALSE, maxgap = 8)
    }
  } else {
    ac$SWC <- NA
  }
  
  # SW_IN
  ac$SW_IN  <- a[, trimws(site_info$SW_IN[id])]
  if (substr(site_info$SW_IN[id], 1, 4) == 'PPFD') {
    # convert PPFD to SW_IN
    ac$SW_IN  <- ac$SW_IN / 2.3
  }
  
  # Special cases: US-Jo2: interpolating SW_IN, otherwise ReddyProc does not work for years before 2013 because of no SW_IN data. 
  if (name_site == "US-Jo2") {
    ac$SW_IN[which(is.na(ac$SW_IN))] <- ac$SW_IN[which(is.na(ac$SW_IN)) + 365*48*4]
  } else if (name_site == 'US-KM4') {
    # use PI gap-filled data
    ac$SW_IN[is.na(ac$SW_IN)] <- a$SW_IN_PI_F[is.na(ac$SW_IN)]
  }
  
  # USTAR
  ac$USTAR  <- a[, trimws(site_info$USTAR[id])]
  
  # RH or VPD
  if (site_info$RH[id] != '') {
    ac$RH  <- a[, trimws(site_info$RH[id])]
    convert_RH_VPD <- TRUE
  } else {
    ac$VPD <- a[, trimws(site_info$VPD[id])]
    convert_RH_VPD <- FALSE
  }
  
  ac$daytime <- a$daytime
  
  ######################################plot and check selected variables###################################
  # plot(ac$NEE)
  # plot(ac$TA)
  # plot(ac$TS)
  # if (site_info$SWC_use[id] == 'YES') {
  #   plot(ac$SWC)
  # }
  # plot(ac$SW_IN)
  # plot(ac$USTAR)
  # if (convert_RH_VPD) {
  #   plot(ac$RH)
  # } else {
  #   plot(ac$VPD)
  # }
  
  ##########################################do Ustar filtering#############################################
  NEE_yearly <- ac %>% group_by(DOY) %>% summarise(NEE=mean(NEE, na.rm=T), TS=mean(TS, na.rm=T))
  NEE_yearly_night <- ac %>% filter(!daytime) %>% group_by(DOY) %>% summarise(NEE=mean(NEE, na.rm=T), TS=mean(TS, na.rm=T))

  # I will use the mean of the first three values, and the mean of the last three values
  tmp    <- NEE_yearly %>% filter(NEE < min(NEE_yearly$NEE, na.rm=T) * 0.2)
  #
  gStart <- as.integer(mean(tmp$DOY[7])) - 4
  gEnd   <- as.integer(mean(tmp$DOY[(nrow(tmp)-6)])) + 4
  
  # In general, growing-season temperature should be > 0 C. 
  gStart <- max(gStart, min(which(NEE_yearly$TS >= 0)))
  if (!is.na(site_info$gStart[id])) {
    gStart <- as.numeric(site_info$gStart[id])
  }
  if (!is.na(site_info$gEnd[id])) {
    gEnd <- as.numeric(site_info$gEnd[id])
  }
  #
  # temperature range of growing season
  tStart <- quantile(tmp$TS, c(0.025), na.rm=T)
  tEnd   <- quantile(tmp$TS, c(0.975), na.rm=T)
  #
  # Create seasonStarts with two columns: one is day, and the other is year!
  years  <- unique(ac$YEAR)
  seasonStarts <- data.frame(DOY = rep(c(gStart, gEnd), length(years)), YEAR=rep(years, each=2))
  #
  # create the data frame
  ac_u       <- data.frame(DateTime = ymd_hm(a$TIMESTAMP_END))   # The first dataset has to be 00:30
  ac_u$NEE   <- ac$NEE
  ac_u$Rg    <- ac$SW_IN
  ac_u$Tair  <- ac$TA
  ac_u$Ustar <- ac$USTAR

  if (convert_RH_VPD) {
    ac_u$rH    <- ac$RH
    ac_u$VPD <- fCalcVPDfromRHandTair(ac_u$rH, ac_u$Tair)
  } else {
    ac_u$VPD <- ac$VPD
  }

  if (dt == 1) {
    EProc <- sEddyProc$new(name_site, ac_u, c('NEE', 'Rg', 'Tair', 'VPD', 'Ustar'), DTS.n=24)
  } else {
    EProc <- sEddyProc$new(name_site, ac_u, c('NEE', 'Rg', 'Tair', 'VPD', 'Ustar'))
  }
  # Rg should not have negative values.
  EProc$sSetLocationInfo(LatDeg = lat_site, LongDeg = long_site, TimeZoneHour = tz$utc_offset_h)
  #
  ac_u$season <- usCreateSeasonFactorYdayYear(
    ac_u$DateTime - 15*60, starts = seasonStarts)  # it sets back 15 min.
  uStarTh <- EProc$sEstUstarThold(seasonFactor = ac_u$season)
  # gap fill NEE, air temperature and soil temperature
  # By default the gap-filling uses annually aggregated estimates of uStar-Threshold.
  # we can also use a different threshold for each of the defined seasons, by calling the two functions
  # EProc$useSeaonsalUStarThresholds()
  # EProc$sGetUstarScenarios()
  # I only want annually aggregated estimated, because some sites have no seasonal estimates
  EProc$sMDSGapFillAfterUstar('NEE', FillAll = FALSE, isVerbose = FALSE)
  ac$NEE_uStar_f <- EProc$sExportResults()$NEE_uStar_f
  #
  ac_u <- ac_u %>% left_join(uStarTh[,3:4], by="season")
  ac$uStarTh <- ac_u$uStar
  #
  #########################################end with Ustar filtering#######################################
  if (site_info$SWC_use[id] == 'YES') {
    a_measure     <- ac %>% filter(!is.na(NEE) & !is.na(TS) & !is.na(SWC) & !is.na(USTAR)) %>% filter(USTAR >= uStarTh)
  } else {
    a_measure     <- ac %>% filter(!is.na(NEE) & !is.na(TS) & !is.na(USTAR)) %>% filter(USTAR >= uStarTh)
  }
  #
  if (name_site %in% c("US-ICh", "US-ICt", "US-ICs")) {
    # tundra sites take advantage of data with little light
    a_measure_night_complete <- a_measure %>% filter(SW_IN < 10 | !daytime) %>% filter(NEE > -5 & NEE < 30)
  } else {
    a_measure_night_complete <- a_measure %>% filter(!daytime) %>% filter(NEE > -5 & NEE < 30)
  }
  
  #
  for (i in years) {
    a_measure_night_year <- a_measure_night_complete %>% filter(YEAR==i)
    # check length of gaps and do a distribution
    if (nrow(a_measure_night_year) > 1) {
      gap <- as.integer(a_measure_night_year$TIMESTAMP[2:nrow(a_measure_night_year)] - a_measure_night_year$TIMESTAMP[1:(nrow(a_measure_night_year)-1)])
      plot(a_measure_night_year$TS, a_measure_night_year$NEE, main=as.character(i), xlim=c(-10, 40), ylim=c(-10, 40))      # TA_1_1_1
      # check air temperature and soil temperature relationship
      plot(a_measure_night_year$TA, a_measure_night_year$TS, main=as.character(i), xlim=c(-10, 40), ylim=c(-5, 40))
      # see gaps each year
      plot(a_measure_night_year$TIMESTAMP, a_measure_night_year$NEE, main=as.character(i), xlim=c(as.POSIXct(paste0(i, "-01-01 00:00:00")), as.POSIXct(paste0(i, "-12-31 00:00:00"))))
    } else {
      print(i)
    }
  }

  # save site-specific TS~NEE curve
  png(paste0("graphs/", name_site, "TS_NEEnight_relation.png"), width = 1200, height = 1200)
  par(mfrow = c(2, 2))
  plot(a_measure_night_complete$TS, a_measure_night_complete$NEE, main=name_site)
  plot(NEE_yearly$DOY, NEE_yearly$NEE)
  plot(NEE_yearly$DOY, NEE_yearly$TS)
  plot(NEE_yearly_night$DOY, NEE_yearly_night$NEE)
  dev.off()
  
  #########################################remove years with big gaps#######################################
  # max gap and total large gap threshold
  if (name_site %in% c("US-ICt", "US-ICh", "US-ICs")) {
    gap_max_thresh <- 60
    gap_total_thresh <- 0.8
  } else {
    gap_max_thresh <- max(31, (gEnd-gStart+1) * 0.2)
    gap_total_thresh <- max(1/3, gap_max_thresh / (gEnd-gStart+1))
  }

  # find long gap years
  a_check_gaps <- a_measure_night_complete[, c("TIMESTAMP", "DOY")]
  yStart <- a_measure_night_complete$YEAR[1]
  yEnd   <- a_measure_night_complete$YEAR[nrow(a_measure_night_complete)]

  # adding the start and end of growing season date
  a_gs_dates <- data.frame(TIMESTAMP = c(as.POSIXct((gStart - 1) * 86400 + as.numeric(dt/2, units = "secs"), origin = paste0(yStart:yEnd, "-01-01"), tz = "UTC"),
                                         as.POSIXct((gEnd - 1) * 86400 + as.numeric(dt/2, units = "secs"), origin = paste0(yStart:yEnd, "-01-01"), tz = "UTC")),
                           DOY = c(rep(gStart, yEnd - yStart + 1), rep(gEnd, yEnd - yStart + 1)))
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
  # bad_TS_years <- a_measure_night_complete %>% filter(between(TS, tStart, tEnd)) %>%
  #   arrange(YEAR, TS) %>% group_by(YEAR) %>%
  #   mutate(gap = TS - lag(TS)) %>%
  #   summarise(max_gap = max(gap, na.rm = TRUE)) %>%
  #   filter(max_gap > max(0.2 * (tEnd - tStart), 1.0)) %>% distinct(YEAR) %>% pull()
  # if (length(bad_TS_years) >= 1) {
  #   good_years <- setdiff(good_years, bad_TS_years)
  # }

  years2remove_automation <- setdiff(ac$YEAR[1]:ac$YEAR[nrow(ac)], good_years)

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

  # Junna temperory output.
  if (setequal(years2remove_automation, years2remove)) {
    print('perfectly matched year2remove')
  } else {
    print('automation does not match year2remove')
    print(years2remove_automation)
    print(years2remove)
  }

  ############################prepare to output two datasets: a_measure_night_complete and ac#######################################
  # we can only use TS > 2C for some sites, because data quality issue below TS < 2C 
  if (name_site %in% c("US-Ha1", "US-GLE")) {
    tStart <- max(tStart, 2.0)
    a_measure_night_complete <- a_measure_night_complete %>% filter(TS >= 2.0)
  } 
  
  a_measure_night_complete <- a_measure_night_complete %>% 
    filter(YEAR %in% good_years) %>%
    dplyr::select(c(YEAR, MONTH, DAY, DOY, HOUR, MINUTE, NEE, TA, TS, SWC))

  iStart = a_measure_night_complete$YEAR[1]
  iEnd   = a_measure_night_complete$YEAR[nrow(a_measure_night_complete)]
  ac <- ac %>% filter(between(YEAR, iStart, iEnd)) %>% dplyr::select(c(YEAR, MONTH, DAY, DOY, HOUR, MINUTE, NEE_uStar_f, TA, TS, SWC, SW_IN, daytime))

  # I have to gap fill TA, TS, NEE;
  T_gf <- ac %>% group_by(DOY, HOUR, MINUTE) %>% summarise(TA_gf=mean(TA, na.rm=T), TS_gf=mean(TS, na.rm=T), NEE_gf=mean(NEE_uStar_f, na.rm=T))  # for gap fill only.
  if (!"TA_gf" %in% colnames(ac)) {
    ac   <- ac %>% left_join(T_gf, by=c("DOY", "HOUR", "MINUTE"))
  }

  if (sum(is.na(ac$TA)) > 0) {
    ac$TA[is.na(ac$TA)] <- ac$TA_gf[is.na(ac$TA)]
  }
  if (sum(is.na(ac$TS)) > 0) {
    ac$TS[is.na(ac$TS)] <- ac$TS_gf[is.na(ac$TS)]
  }
  if (sum(is.na(ac$NEE_uStar_f)) > 0) {
    ac$NEE_uStar_f[is.na(ac$NEE_uStar_f)] <- ac$NEE_gf[is.na(ac$NEE_uStar_f)]
  }

  # save the ac and a_measure_night_complete data:
  saveRDS(ac, file=file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.RDS')))
  saveRDS(a_measure_night_complete, file=file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.RDS')))

  feature_gs <- bind_rows(feature_gs, data.frame(site_ID=name_site, gStart=gStart, gEnd=gEnd, tStart=max(tStart, 0.0), tEnd=tEnd, nyear=length(good_years)))
}
# output growing season features
# write.csv(feature_gs, file='data/growing_season_feature_AmeriFlux.csv', row.names = F)

