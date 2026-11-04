library(dplyr)
library(tibble)
library(lubridate)

#----------------------------Two functions to calculate soil temperature-----------------
# first one: estimate TS at sites without TS measurements. 
# 19 sites
no_TS_sites <- c(
  "AT-Zoe", "BR-Ma2", "BR-Sa1", "CH-Cha", "FR-FBn", "ID-JOP", "ID-PaD", "ID-Pag", "JP-NsC","JP-NsM",
  "MY-LHP", "US-Cwt", "US-HBK", "US-Mpj", "US-Seg", "US-Ses", "US-Vcm", "US-Vcp", "US-Wjs"
)

reconstruct_missing_TS <- function(
    a, name_site, ts_config, long_mean_days = 8,
    shuttle_dir) {
  
  candidates <- list(
    "AT-Zoe" = c("CA-Gro", "US-Syv", "JP-MMF"),
    "BR-Ma2" = c("CR-AqC", "AU-Cow", "AU-Ctr"),
    "BR-Sa1" = c("CR-AqC", "AU-Cow", "AU-Ctr"),
    "CH-Cha" = c("CH-Oe3", "CH-Oe1", "DE-RuR"),
    "FR-FBn" = "IT-Cp2",
    "ID-JOP" = c("AU-Ctr", "AU-Cow", "CR-AqC"),
    "ID-PaD" = c("AU-Ctr", "AU-Cow", "CR-AqC"),
    "ID-Pag" = c("AU-Ctr", "AU-Cow", "CR-AqC"),
    "JP-NsC" = c("US-Kon", "US-KFS"),
    "JP-NsM" = c("US-Kon", "US-KFS"),
    "MY-LHP" = c("AU-Ctr", "AU-Cow", "CR-AqC"),
    "US-Cwt" = c("CA-TPD", "CA-Cbo", "US-Ha1"),
    "US-HBK" = c("US-Bar", "US-Ha1", "CA-Cbo"),
    "US-Mpj" = c("US-SRS", "US-SRM"),
    "US-Seg" = c("US-Wkg", "US-SRG"),
    "US-Ses" = c("US-Whs", "US-Rws"),
    "US-Vcm" = c("CA-TP4", "CA-TP3", "CA-LP1"),
    "US-Vcp" = c("CA-TP4", "CA-TP3", "CA-LP1"),
    "US-Wjs" = c("US-SRS", "US-SRM")
  )
  
  if (!name_site %in% names(candidates))
    stop(name_site, " is not a site without soil-temperature measurements.")
  
  if (!all(c("TIMESTAMP", "TA_F_MDS") %in% names(a)))
    stop(name_site, " requires TIMESTAMP and TA_F_MDS.")
  
  files <- list.files(
    shuttle_dir, "_(HH|HR)_.*[.]csv$",
    full.names = TRUE
  )
  files <- files[!startsWith(basename(files), "._")]
  
  prepare_donor <- function(site) {
    
    file <- files[
      grepl(paste0("_", site, "_"), basename(files), fixed = TRUE)
    ]
    
    if (!length(file)) stop("No Shuttle file found for ", site)
    
    hh_file <- file[grepl("_HH_", basename(file), fixed = TRUE)]
    file <- if (length(hh_file)) hh_file[1] else file[1]
    
    cfg <- ts_config[ts_config$site == site, , drop = FALSE]
    
    if (nrow(cfg) > 1) stop(site, " has multiple configuration rows.")
    
    if (!nrow(cfg)) {
      cfg <- data.frame(
        site = site,
        ts_variable = "TS_F_MDS_1",
        training_start = as.Date(NA),
        training_end = as.Date(NA),
        ta_min = 0
      )
    }
    
    ts_var <- cfg$ts_variable[1]
    required <- c("TIMESTAMP_START", "TA_F_MDS", ts_var)
    
    if (!all(required %in% names(data.table::fread(file, nrows = 0))))
      stop(site, " is missing required variables.")
    
    x <- data.table::fread(
      file, select = required,
      na.strings = c("NA", "NaN", "-9999", "-9999.0")
    )
    
    x[, TIMESTAMP := as.POSIXct(
      sprintf("%012s", as.character(TIMESTAMP_START)),
      format = "%Y%m%d%H%M", tz = "UTC"
    )]
    
    x <- x[!is.na(TIMESTAMP)]
    data.table::setorder(x, TIMESTAMP)
    
    dt <- median(
      as.numeric(diff(x$TIMESTAMP), units = "hours"),
      na.rm = TRUE
    )
    
    if (!is.finite(dt) || dt <= 0)
      stop("Invalid time interval for ", site)
    
    n_day <- round(24 / dt)
    
    x[, `:=`(
      donor_site = site,
      TS_observed = get(ts_var),
      TA_current = TA_F_MDS,
      TA_previous_1day = previous_mean(TA_F_MDS, n_day),
      TA_previous_8day = previous_mean(
        TA_F_MDS, n_day * long_mean_days
      )
    )]
    
    x[is.na(TA_previous_1day),
      TA_previous_1day := TA_current]
    
    x[is.na(TA_previous_8day),
      TA_previous_8day := TA_current]
    
    keep <- rep(TRUE, nrow(x))
    date <- as.Date(x$TIMESTAMP, tz = "UTC")
    
    if (!is.na(cfg$training_start[1]))
      keep <- keep & date >= cfg$training_start[1]
    
    if (!is.na(cfg$training_end[1]))
      keep <- keep & date <= cfg$training_end[1]
    
    if (!is.na(cfg$ta_min[1]))
      keep <- keep & x$TA_F_MDS > cfg$ta_min[1]
    
    vars <- c(
      "TS_observed", "TA_current",
      "TA_previous_1day", "TA_previous_8day"
    )
    
    x[
      keep & complete.cases(x[, ..vars]),
      c("donor_site", vars),
      with = FALSE
    ]
  }
  
  donor_data <- lapply(candidates[[name_site]], prepare_donor)
  names(donor_data) <- candidates[[name_site]]
  
  training_data <- data.table::rbindlist(donor_data)
  
  model <- lm(
    TS_observed ~ TA_current +
      TA_previous_1day + TA_previous_8day,
    data = training_data
  )
  
  dt <- median(
    as.numeric(diff(a$TIMESTAMP), units = "hours"),
    na.rm = TRUE
  )
  
  if (!is.finite(dt) || dt <= 0)
    stop("Invalid time interval for ", name_site)
  
  n_day <- round(24 / dt)
  
  predictors <- data.frame(
    TA_current = a$TA_F_MDS,
    TA_previous_1day = previous_mean(a$TA_F_MDS, n_day),
    TA_previous_8day = previous_mean(
      a$TA_F_MDS, n_day * long_mean_days
    )
  )
  
  predictors$TA_previous_1day[
    is.na(predictors$TA_previous_1day)
  ] <- predictors$TA_current[
    is.na(predictors$TA_previous_1day)
  ]
  
  predictors$TA_previous_8day[
    is.na(predictors$TA_previous_8day)
  ] <- predictors$TA_current[
    is.na(predictors$TA_previous_8day)
  ]
  
  a$TS_F_MDS_1 <- predict(model, newdata = predictors)
  
  if ("TA_F_MDS_QC" %in% names(a))
    a$TS_F_MDS_1_QC <- a$TA_F_MDS_QC
  
  list(
    data = a,
    model = model,
    training_n = nrow(training_data),
    n_per_day = n_day,
    donors_used = candidates[[name_site]],
    donor_training_n = vapply(donor_data, nrow, integer(1)),
    coefficients = coef(model)
  )
}


# Second one: estimate TS at sites with some TS measurements

# Sites using all TS_F_MDS_1 observations with TA_F_MDS > 0
all_data_sites <- c(
  "AU-DaS","AU-Dry","AU-Lit","AU-Ync","US-Dk2","DE-HoH","FR-Hes","JP-Yms",
  "BE-Bra","CA-Cbo","CA-Gro","CA-Mer","CA-Obs","CA-TP3","CH-Lae","DE-RuC","DE-SfS",
  "IT-Ren","NL-Loo","US-Bar","US-BZB","US-BZF","US-BZS","US-CMW","US-GLE","US-Ha2","US-IB2",
  "US-Jo2","US-KL2","US-Kon","US-LL1","US-MBP","US-Myb","US-NC4","US-Tw1","US-ICt","BE-Dor",
  "CA-TP4","UK-AMo","RU-Fyo","ZA-Kru","IT-Tor","AT-Neu","BJ-Nhu","CH-Oe3","CN-YaS",
  "DE-GsB","ES-Amo","FI-Var","FR-LBr","JP-Fjy","RU-Ege","SN-Dhr","US-TKs","CH-Rtb","CL-SDP",
  "CN-HaM","IT-Cpz","NZ-ADw","US-Dk3","DE-Akm","US-Los","US-SRG","CA-Man","US-Ced","US-Ho1",
  "CZ-RAJ","US-Ha1","US-BZo","US-SRS","US-Oho","CA-LP1","CA-Ca3","CA-ARB","CA-ARF"
)

config_all <- tibble(
  site=unique(all_data_sites),
  ts_variable="TS_F_MDS_1",
  training_start=as.Date(NA),
  training_end=as.Date(NA),
  ta_min=0
)

# Sites using restricted training periods
config_period <- tribble(
  ~site,   ~ts_variable, ~training_start, ~training_end, ~ta_min,
  "AU-Rob","TS_F_MDS_1","2013-01-01","2022-09-12",NA_real_,
  "AU-TTE","TS_F_MDS_1","2012-01-01","2019-12-26",NA_real_,
  "AU-Whr","TS_F_MDS_1","2011-01-01","2016-12-27",NA_real_,
  "AU-Wrr","TS_F_MDS_1","2013-01-01","2018-12-28",NA_real_,
  "BJ-Bfg","TS_F_MDS_1","2008-01-01","2013-12-27",NA_real_,
  "ES-Agu","TS_F_MDS_1","2006-01-01","2011-12-28",NA_real_,
  "US-Jo1","TS_F_MDS_1","2010-01-01","2015-12-28",NA_real_,
  "AU-Wom","TS_F_MDS_1","2024-10-30","2025-12-31",NA_real_,
  "TH-Mae","TS_F_MDS_1","2014-02-16","2016-12-31",NA_real_,
  "CN-HeM","TS_F_MDS_1","2018-09-15","2024-12-31",NA_real_,
  "CN-Zha","TS_F_MDS_1","2017-09-14","2023-12-31",NA_real_,
  "EE-Stg","TS_F_MDS_1","2023-06-01","2024-12-31",NA_real_,
  "IT-Cp2","TS_F_MDS_1","2012-01-01","2023-12-23",NA_real_,
  "JP-Khw","TS_F_MDS_1","2017-02-09","2021-12-31",NA_real_,
  "CA-SCC","TS_F_MDS_1","2023-11-03","2024-12-31",0,
  "RU-SkP","TS_F_MDS_1","2009-09-14","2014-12-31",0,
  "US-Elm","TS_F_MDS_1","2008-01-01","2013-09-14",0,
  "US-Ho2","TS_F_MDS_1","2008-01-01","2022-12-31",0,
  "AU-Cow","TS_F_MDS_1","2011-01-01","2020-12-31",0,
  "FR-Pue","TS_F_MDS_1","2022-01-01",NA_character_,0,
  "FR-Bil","TS_F_MDS_1","2020-01-01",NA_character_,0,
  "FI-Sod","TS_F_MDS_1","2006-01-01",NA_character_,0,
  "DE-RuW","TS_F_MDS_1","2019-01-01",NA_character_,0,
  "US-Ro4","TS_F_MDS_1","2015-01-01",NA_character_,0,
  "AU-How","TS_F_MDS_1","2004-01-01","2024-12-31",0,
  "US-PFa","TS_F_MDS_1","2023-01-01","2024-12-31",0,
  "CR-AqC","TS_F_MDS_1","2016-01-01","2016-12-31",0,
  "IT-SRo","TS_F_MDS_1","2000-01-01","2012-12-31",0,
  "CN-Dsh","TS_F_MDS_1","2014-01-01","2021-12-31",0,
  "US-Prr","TS_F_MDS_1","2014-01-01","2022-12-31",0,
  "US-ICs","TS_F_MDS_1","2010-01-01","2025-12-31",0
) %>%
  mutate(
    training_start=as.Date(training_start),
    training_end=as.Date(training_end)
  )

# Sites using different TS layer
config_layer <- tribble(
  ~site,    ~ts_variable,  ~training_start, ~training_end, ~ta_min,
  "US-EML","TS_F_MDS_2",NA_character_,NA_character_, 0,
  "CH-Dav","TS_F_MDS_2",NA_character_,NA_character_,0,
  "CZ-Stn","TS_F_MDS_2",NA_character_,NA_character_,0,
  "DE-Hte","TS_F_MDS_2","2012-01-01","2017-12-31",0
) %>%
  mutate(
    training_start=as.Date(training_start),
    training_end=as.Date(training_end)
  )

ts_config <- bind_rows(config_all,config_period,config_layer)

if(anyDuplicated(ts_config$site)){
  stop("Some sites occur more than once in ts_config.")
}

# Mean of the preceding n observations, excluding the current observation
previous_mean <- function(x,n){
  xlag <- c(NA_real_,head(x,-1))
  cs <- c(0,cumsum(ifelse(is.na(xlag),0,xlag)))
  cn <- c(0,cumsum(!is.na(xlag)))
  i <- seq_along(xlag)
  left <- pmax(1L,i-n+1L)
  total <- cs[i+1L]-cs[left]
  count <- cn[i+1L]-cn[left]
  result <- total/count
  result[count==0] <- NA_real_
  result
}

reconstruct_TS <- function(a,name_site,ts_config,long_mean_days=8){
  cfg <- ts_config %>% filter(site==name_site)
  
  if(nrow(cfg)==0){
    return(list(data=a,model=NULL,training_n=0))
  }
  
  if(nrow(cfg)>1){
    stop(paste(name_site,"has multiple configuration rows."))
  }
  
  required <- c("TIMESTAMP","TA_F_MDS","TA_F_MDS_QC",cfg$ts_variable)
  missing_columns <- setdiff(required,names(a))
  
  if(length(missing_columns)>0){
    stop(paste(name_site,"is missing:",paste(missing_columns,collapse=", ")))
  }
  
  dt_hours <- median(as.numeric(diff(a$TIMESTAMP),units="hours"),na.rm=TRUE)
  n_per_day <- round(24/dt_hours)
  n_long <- n_per_day*long_mean_days
  
  d <- a %>%
    mutate(
      .date=as.Date(TIMESTAMP),
      TA_current=TA_F_MDS,
      TA_previous_1day=previous_mean(TA_F_MDS,n_per_day),
      TA_previous_long=previous_mean(TA_F_MDS,n_long)
    )
  
  # Allow prediction at the beginning of the record
  d$TA_previous_1day[is.na(d$TA_previous_1day)] <-
    d$TA_current[is.na(d$TA_previous_1day)]
  
  d$TA_previous_long[is.na(d$TA_previous_long)] <-
    d$TA_current[is.na(d$TA_previous_long)]
  
  keep <- rep(TRUE,nrow(d))
  
  if(!is.na(cfg$training_start)){
    keep <- keep & d$.date>=cfg$training_start
  }
  
  if(!is.na(cfg$training_end)){
    keep <- keep & d$.date<=cfg$training_end
  }
  
  if(!is.na(cfg$ta_min)){
    keep <- keep & d$TA_F_MDS>cfg$ta_min
  }
  
  model_variables <- c(
    cfg$ts_variable,
    "TA_current",
    "TA_previous_1day",
    "TA_previous_long"
  )
  
  keep <- keep & complete.cases(d[,model_variables])
  training_data <- d[keep,]
  
  if(nrow(training_data)<30){
    stop(paste(name_site,"has too few training observations."))
  }
  
  model_formula <- reformulate(
    c("TA_current","TA_previous_1day","TA_previous_long"),
    response=cfg$ts_variable
  )
  
  mod_TS <- lm(model_formula,data=training_data)
  
  # Replace the complete TS_F_MDS_1 time series
  a$TS_F_MDS_1 <- predict(mod_TS,newdata=d)
  a$TS_F_MDS_1_QC <- a$TA_F_MDS_QC
  
  list(
    data=a,
    model=mod_TS,
    training_n=nrow(training_data),
    n_per_day=n_per_day
  )
}


#-----------------A function to estimate the start and end of ER active period----------------------------
# A function to estimate the start and end of ER active period
estimate_active_ER_period <- function(a, site_id, ta_var,
                                      ts_var = NULL,
                                      southern_hemisphere = FALSE,
                                      threshold_fraction = 0.20,
                                      temp_threshold = 2,
                                      min_observed_doy = 300L) {
  # This function expects a to already contain YEAR, DOY, and NIGHT, as in
  # 01_02c_filter_high_quality_night_respiration_FluxShuttle.R. For southern
  # sites, DOY must already have been shifted there: DOY < 183 was increased
  # by 366. Returned gStart/gEnd use exactly that same DOY convention.
  stopifnot(is.data.frame(a), length(site_id) == 1L,
            length(ta_var) == 1L, is.character(ta_var),
            length(southern_hemisphere) == 1L,
            !is.na(southern_hemisphere),
            is.finite(threshold_fraction), threshold_fraction > 0,
            threshold_fraction < 1,
            is.numeric(temp_threshold), length(temp_threshold) == 1L,
            is.finite(temp_threshold),
            is.finite(min_observed_doy), min_observed_doy >= 30)
  
  skip <- function(reason, nee_product = NA_character_,
                   temperature_source = NA_character_,
                   trimmed_days = NA_integer_) {
    message(sprintf('%s: skipped active ER period (%s).', site_id, reason))
    list(status = 'skipped', reason = reason, site_id = site_id,
         gStart = NA_integer_, gEnd = NA_integer_,
         trimmed_days = trimmed_days, nee_product = nee_product,
         temperature_source = temperature_source,
         threshold_fraction = threshold_fraction,
         temp_threshold = temp_threshold)
  }
  
  needed <- c('YEAR', 'DOY', 'NIGHT', ta_var)
  if (!all(needed %in% names(a))) {
    return(skip(paste('missing columns:',
                      paste(setdiff(needed, names(a)), collapse = ', '))))
  }
  
  # Prefer the VUT product when its REF and QC columns exist. CUT is used
  # only when the VUT pair is unavailable, matching the requested rule.
  if (all(c('NEE_VUT_REF', 'NEE_VUT_REF_QC') %in% names(a))) {
    nee_var <- 'NEE_VUT_REF'
    qc_var <- 'NEE_VUT_REF_QC'
    nee_product <- 'VUT'
  } else if (all(c('NEE_CUT_REF', 'NEE_CUT_REF_QC') %in% names(a))) {
    nee_var <- 'NEE_CUT_REF'
    qc_var <- 'NEE_CUT_REF_QC'
    nee_product <- 'CUT'
  } else {
    return(skip('neither VUT nor CUT has both REF and QC columns'))
  }
  
  # Undo the script's +366 transformation for fitting a cyclic annual GAM.
  # We apply that same transformation again to predictions before finding
  # gStart/gEnd. This deliberately retains the script's leap-year DOY usage.
  doy <- as.integer(a$DOY)
  if (southern_hemisphere) doy[doy > 366L] <- doy[doy > 366L] - 366L
  year <- as.integer(a$YEAR)
  valid_day <- is.finite(doy) & doy >= 1L & doy <= 366L & is.finite(year)
  
  # Make a site-year daily mean, then a climatological DOY mean. Each year
  # contributes at most one daily value for a DOY, regardless of whether
  # its source has half-hourly or hourly records.
  daily_climatology <- function(values, keep) {
    values <- suppressWarnings(as.numeric(values))
    values[values == -9999] <- NA_real_
    good <- valid_day & keep & is.finite(values)
    if (!any(good)) return(NULL)
    dat <- data.frame(year = year[good], doy = doy[good], value = values[good])
    daily <- stats::aggregate(value ~ year + doy, dat, mean)
    mean_by_doy <- stats::aggregate(value ~ doy, daily, mean)
    n_years <- stats::aggregate(value ~ doy, daily, length)
    mean_by_doy$n_years <- n_years$value[match(mean_by_doy$doy, n_years$doy)]
    mean_by_doy
  }
  
  qc <- suppressWarnings(as.numeric(a[[qc_var]]))
  night <- suppressWarnings(as.numeric(a$NIGHT))
  nee_clim <- daily_climatology(a[[nee_var]], night == 1 & qc == 0)
  if (is.null(nee_clim) || nrow(nee_clim) < min_observed_doy) {
    return(skip(sprintf('fewer than %d DOYs with QC=0 nighttime NEE',
                        min_observed_doy), nee_product))
  }
  
  # Use soil temperature when provided and sufficiently complete. Otherwise
  # use air temperature. Temperature climatology uses all hours, not only
  # records with observed nighttime NEE.
  temperature_source <- ta_var
  temp_clim <- NULL
  if (!is.null(ts_var) && length(ts_var) == 1L && ts_var %in% names(a)) {
    soil_clim <- daily_climatology(a[[ts_var]], rep(TRUE, nrow(a)))
    if (!is.null(soil_clim) && nrow(soil_clim) >= min_observed_doy) {
      temp_clim <- soil_clim
      temperature_source <- ts_var
    }
  }
  if (is.null(temp_clim)) {
    temp_clim <- daily_climatology(a[[ta_var]], rep(TRUE, nrow(a)))
  }
  if (is.null(temp_clim) || nrow(temp_clim) < min_observed_doy) {
    return(skip(sprintf('fewer than %d DOYs with usable temperature',
                        min_observed_doy), nee_product, temperature_source))
  }
  
  if (!requireNamespace('mgcv', quietly = TRUE)) {
    return(skip('R package mgcv is not installed', nee_product,
                temperature_source))
  }
  # Cyclic cubic smooths join DOY 366 back to DOY 1. Year counts weight
  # climatological DOYs, and REML selects the amount of smoothing.
  fit_annual <- function(clim) {
    # A constant series has no seasonal smooth to estimate; mgcv's REML
    # optimization can fail in that case, so preserve the constant directly.
    if (stats::sd(clim$value) < 1e-10) {
      return(rep(mean(clim$value), 366L))
    }
    s <- mgcv::s
    k <- min(30L, nrow(clim) - 1L)
    model <- mgcv::gam(value ~ s(doy, bs = 'cc', k = k), data = clim,
                       weights = n_years, method = 'REML',
                       knots = list(doy = c(0.5, 366.5)))
    as.numeric(stats::predict(model, newdata = data.frame(doy = 1:366)))
  }
  fitted <- tryCatch(list(nee = fit_annual(nee_clim),
                          temp = fit_annual(temp_clim)),
                     error = function(e) e)
  if (inherits(fitted, 'error')) {
    return(skip(paste('cyclic GAM failed:', conditionMessage(fitted)),
                nee_product, temperature_source))
  }
  
  # Order predictions on the exact axis used by the downstream script.
  # North: 1...366. South: 183...548 (Jan-Jun appear at 367...548).
  axis <- if (southern_hemisphere) 183:548 else 1:366
  calendar_doy <- ifelse(axis > 366L, axis - 366L, axis)
  nee_fit <- fitted$nee[calendar_doy]
  temp_fit <- fitted$temp[calendar_doy]
  threshold <- min(nee_fit) + threshold_fraction *
    (max(nee_fit) - min(nee_fit))
  eligible <- nee_fit > threshold & temp_fit > temp_threshold
  n <- length(axis)
  
  if (!any(eligible)) {
    return(skip(sprintf('no DOY exceeds both the NEE and %.1f C thresholds',
                        temp_threshold),
                nee_product, temperature_source))
  }
  if (all(eligible)) {
    return(list(status = 'ok', reason = NA_character_, site_id = site_id,
                gStart = min(axis), gEnd = max(axis), trimmed_days = 0L,
                offStart = NA_integer_, offEnd = NA_integer_,
                nee_product = nee_product,
                temperature_source = temperature_source,
                threshold_fraction = threshold_fraction,
                temp_threshold = temp_threshold,
                nee_threshold = threshold))
  }
  
  # Find each circular run where at least one criterion is not met. Select
  # the run with the lowest MEAN fitted temperature as the true off-season.
  # Other, warmer dips remain inside one continuous active-season envelope;
  # this is what keeps multiple NEE peaks together.
  inactive <- !eligible
  starts <- which(inactive & !inactive[c(n, seq_len(n - 1L))])
  runs <- lapply(starts, function(start) {
    indices <- integer(0)
    j <- start
    while (inactive[j]) {
      indices <- c(indices, j)
      j <- if (j == n) 1L else j + 1L
    }
    indices
  })
  mean_temp <- vapply(runs, function(indices) mean(temp_fit[indices]), 0.0)
  off <- runs[[which.min(mean_temp)]]
  active <- rep(TRUE, n)
  active[off] <- FALSE
  
  # A circular active season can appear as two pieces on a fixed DOY axis.
  # Retain its longer linear piece; do not return a cross-year interval.
  # Always return the site, but report how many active days were excluded.
  pieces <- rle(active)
  piece_starts <- cumsum(c(1L, head(pieces$lengths, -1L)))
  active_pieces <- which(pieces$values)
  chosen <- active_pieces[which.max(pieces$lengths[active_pieces])]
  trimmed_days <- as.integer(sum(active) - pieces$lengths[chosen])
  first <- piece_starts[chosen]
  last <- first + pieces$lengths[chosen] - 1L
  if (trimmed_days > 0L) {
    message(sprintf('%s: trimmed %d active days at the DOY boundary.',
                    site_id, trimmed_days))
  }
  
  list(status = 'ok', reason = NA_character_, site_id = site_id,
       gStart = axis[first], gEnd = axis[last], trimmed_days = trimmed_days,
       offStart = axis[off[1L]], offEnd = axis[tail(off, 1L)],
       nee_product = nee_product, temperature_source = temperature_source,
       threshold_fraction = threshold_fraction,
       temp_threshold = temp_threshold, nee_threshold = threshold)
}


#------------a function to add red point
mark_period <- function(doy, value, gStart, gEnd) {
  boundaries <- c(gStart, gEnd)
  if (!all(is.finite(boundaries))) return(invisible(NULL))
  
  for (i in 1:2) {
    abline(v = boundaries[i], col = "red", lty = c(2, 3)[i], lwd = 2)
  }
  
  valid <- is.finite(doy) & is.finite(value)
  if (sum(valid) >= 2) {
    y <- approx(doy[valid], value[valid],
                xout = boundaries, ties = mean)$y
    points(boundaries, y, col = "red", pch = c(16, 17), cex = 1.5)
  }
}


#------------------------a function to identify control year----------------------
select_control_year <- function(
    ac,
    night_nee,
    gStart,
    gEnd,
    years = NULL,
    long_gap_threshold = 14,
    alternative_gap_weight = 4
) {
  
  # Required input columns
  required_ac <- c(
    "growing_year",
    "DOY",
    "TS",
    "daytime"
  )
  
  required_night <- c(
    "growing_year",
    "DOY",
    "NEE"
  )
  
  missing_ac <- setdiff(
    required_ac,
    names(ac)
  )
  
  missing_night <- setdiff(
    required_night,
    names(night_nee)
  )
  
  if (length(missing_ac) > 0) {
    stop(
      "ac is missing: ",
      paste(missing_ac, collapse = ", ")
    )
  }
  
  if (length(missing_night) > 0) {
    stop(
      "night_nee is missing: ",
      paste(missing_night, collapse = ", ")
    )
  }
  
  if (
    !is.finite(gStart) ||
    !is.finite(gEnd) ||
    gStart > gEnd
  ) {
    stop(
      "gStart and gEnd must be finite, ",
      "and gStart must be <= gEnd."
    )
  }
  
  # Use growing years available in both data frames
  if (is.null(years)) {
    years <- intersect(
      unique(ac$growing_year),
      unique(night_nee$growing_year)
    )
  }
  
  years <- sort(
    unique(years[is.finite(years)])
  )
  
  if (length(years) == 0) {
    stop(
      "No common growing years are available."
    )
  }
  
  # ==================================================
  # 1. Mean growing-season soil temperature
  # ==================================================
  
  yearly_temperature <- ac |>
    dplyr::filter(
      growing_year %in% years,
      dplyr::between(DOY, gStart, gEnd)
    ) |>
    dplyr::group_by(growing_year) |>
    dplyr::summarise(
      # mean_TS = mean(
      #   TS,
      #   na.rm = TRUE
      # ),
      mean_TS = mean(TS), 
      .groups = "drop"
    ) |>
    dplyr::filter(
      is.finite(mean_TS)
    )
  
  if (nrow(yearly_temperature) == 0) {
    stop(
      "No finite growing-season soil temperature is available."
    )
  }
  
  years <- intersect(
    years,
    yearly_temperature$growing_year
  )
  
  # ==================================================
  # 2. Obtain every actual growing-season day from AC
  # ==================================================
  
  expected_nights <- ac |>
    dplyr::filter(
      growing_year %in% years,
      dplyr::between(DOY, gStart, gEnd)
    ) |>
    dplyr::group_by(
      growing_year,
      DOY
    ) |>
    dplyr::summarise(
      # Check whether the daytime classification exists
      daytime_known = any(
        !is.na(daytime)
      ),
      
      # Does this date contain any nighttime period?
      has_night = any(
        !daytime,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      unknown_daytime = !daytime_known,
      
      # If daytime is completely unknown, retain the
      # date as an expected night. This prevents an
      # unknown classification from hiding a gap.
      has_night = dplyr::if_else(
        unknown_daytime,
        TRUE,
        has_night
      )
    )
  
  # ==================================================
  # 3. Determine whether each night has observed NEE
  # ==================================================
  
  observed_nights <- night_nee |>
    dplyr::filter(
      growing_year %in% years,
      dplyr::between(DOY, gStart, gEnd)
    ) |>
    dplyr::group_by(
      growing_year,
      DOY
    ) |>
    dplyr::summarise(
      # A night is observed if it contains at least
      # one non-missing NEE value.
      observed = any(
        !is.na(NEE)
      ),
      .groups = "drop"
    )
  
  # Join observed NEE to every date available in AC
  night_status <- expected_nights |>
    dplyr::left_join(
      observed_nights,
      by = c(
        "growing_year",
        "DOY"
      )
    ) |>
    dplyr::mutate(
      observed = tidyr::replace_na(
        observed,
        FALSE
      ),
      
      # A date is a gap only if:
      # 1. the date contains nighttime, and
      # 2. no NEE observation is available.
      is_gap = has_night & !observed
    ) |>
    dplyr::arrange(
      growing_year,
      DOY
    ) |>
    dplyr::group_by(growing_year) |>
    dplyr::mutate(
      # Every non-gap date breaks a continuous gap.
      # This includes polar-day dates.
      gap_run = cumsum(!is_gap)
    ) |>
    dplyr::ungroup()
  
  # ==================================================
  # 4. Calculate individual continuous gaps
  # ==================================================
  
  gap_runs <- night_status |>
    dplyr::filter(is_gap) |>
    dplyr::group_by(
      growing_year,
      gap_run
    ) |>
    dplyr::summarise(
      gap_start_DOY = min(DOY),
      gap_end_DOY = max(DOY),
      gap_days = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      growing_year,
      gap_start_DOY
    )
  
  # ==================================================
  # 5. Summarize nighttime coverage by growing year
  # ==================================================
  
  yearly_counts <- night_status |>
    dplyr::group_by(growing_year) |>
    dplyr::summarise(
      growing_season_days = dplyr::n(),
      
      # Days without any nighttime period
      polar_day_days = sum(
        !has_night
      ),
      
      unknown_daytime_days = sum(
        unknown_daytime
      ),
      
      # Number of dates that actually contain nighttime
      expected_nights = sum(
        has_night
      ),
      
      # Expected nights containing observed NEE
      observed_nights = sum(
        has_night & observed
      ),
      
      # Expected nights without observed NEE
      total_gap_days = sum(
        is_gap
      ),
      
      .groups = "drop"
    )
  
  # ==================================================
  # 6. Calculate gap length and gap burden
  # ==================================================
  
  yearly_gaps <- gap_runs |>
    dplyr::group_by(growing_year) |>
    dplyr::summarise(
      gap_lengths_days = paste(
        gap_days,
        collapse = ", "
      ),
      
      max_gap_days = max(
        gap_days
      ),
      
      # Squaring the gap length penalizes long,
      # continuous gaps more strongly.
      gap_burden = sum(
        gap_days^2
      ),
      
      .groups = "drop"
    )
  
  # ==================================================
  # 7. Combine temperature and gap information
  # ==================================================
  
  yearly_table <- yearly_temperature |>
    dplyr::left_join(
      yearly_counts,
      by = "growing_year"
    ) |>
    dplyr::left_join(
      yearly_gaps,
      by = "growing_year"
    ) |>
    dplyr::mutate(
      gap_lengths_days = tidyr::replace_na(
        gap_lengths_days,
        "none"
      ),
      
      max_gap_days = tidyr::replace_na(
        max_gap_days,
        0L
      ),
      
      gap_burden = tidyr::replace_na(
        gap_burden,
        0
      )
    )
  
  # ==================================================
  # 8. Initial selection using equal weights
  # ==================================================
  
  target_TS <- median(
    yearly_table$mean_TS,
    na.rm = TRUE
  )
  
  yearly_table <- yearly_table |>
    dplyr::mutate(
      target_TS = target_TS,
      
      temperature_distance = abs(
        mean_TS - target_TS
      ),
      
      temperature_rank = rank(
        temperature_distance,
        ties.method = "average"
      ),
      
      gap_burden_rank = rank(
        gap_burden,
        ties.method = "average"
      ),
      
      max_gap_rank = rank(
        max_gap_days,
        ties.method = "average"
      ),
      
      # Equal weights for:
      # 1. temperature representativeness,
      # 2. gap burden,
      # 3. maximum gap.
      initial_score =
        temperature_rank +
        gap_burden_rank +
        max_gap_rank
    )
  
  initial_result <- yearly_table |>
    dplyr::arrange(
      initial_score,
      gap_burden,
      temperature_distance,
      growing_year
    ) |>
    dplyr::slice(1)
  
  initial_year <- initial_result$growing_year
  final_year <- initial_year
  
  selection_reason <-
    "initial_equal_weight_selection"
  
  # ==================================================
  # 9. Second-stage long-gap adjustment
  # ==================================================
  
  if (
    initial_result$max_gap_days >
    long_gap_threshold
  ) {
    
    alternatives <- yearly_table |>
      dplyr::filter(
        # Do not select the same initial year
        growing_year != initial_year,
        
        # The alternative must satisfy the
        # maximum-gap requirement.
        max_gap_days <= long_gap_threshold,
        
        # Avoid dates with unknown nighttime status.
        unknown_daytime_days == 0
      ) |>
      dplyr::mutate(
        # Apply the larger gap-burden weight only
        # when choosing among acceptable alternatives.
        alternative_score =
          temperature_rank +
          alternative_gap_weight *
          gap_burden_rank +
          max_gap_rank
      ) |>
      dplyr::arrange(
        alternative_score,
        gap_burden,
        temperature_distance,
        growing_year
      )
    
    if (nrow(alternatives) > 0) {
      
      final_year <-
        alternatives$growing_year[1]
      
      selection_reason <-
        "replacement_for_gap_longer_than_threshold"
      
    } else {
      
      # Retain the initial year if no alternative
      # satisfies the maximum-gap requirement.
      selection_reason <-
        "initial_year_kept_no_acceptable_alternative"
    }
  }
  
  # ==================================================
  # 10. Mark the initial and final control years
  # ==================================================
  
  yearly_table <- yearly_table |>
    dplyr::mutate(
      initial_selected =
        growing_year == initial_year,
      
      final_selected =
        growing_year == final_year
    ) |>
    dplyr::arrange(growing_year)
  
  # ==================================================
  # 11. Return results
  # ==================================================
  
  list(
    control_year = final_year,
    initial_control_year = initial_year,
    selection_changed = final_year != initial_year,
    selection_reason = selection_reason,
    long_gap_threshold = long_gap_threshold,
    alternative_gap_weight = alternative_gap_weight,
    yearly_table = yearly_table,
    gap_runs = gap_runs
  )
}



repair_TS_periods <- function(a, bad_periods) {
  
  #------------------------------------------------------------
  # 1. Determine how many observations occur each day
  #
  # Half-hourly data: n_per_day = 48
  # Hourly data:      n_per_day = 24
  #------------------------------------------------------------
  
  interval_hours <- median(
    as.numeric(diff(a$TIMESTAMP), units = "hours"),
    na.rm = TRUE
  )
  
  n_per_day <- round(24 / interval_hours)
  
  
  #------------------------------------------------------------
  # 2. Calculate the trailing air-temperature means
  #
  # TA_24h uses the current and previous 24 hours.
  # TA_8day uses the current and previous eight days.
  #
  # TA_F is used because it is the consolidated FLUXNET
  # air-temperature product.
  #------------------------------------------------------------
  
  a$TA_24h <- data.table::frollmean(
    a$TA_F,
    n = n_per_day,
    align = "right",
    na.rm = TRUE
  )
  
  a$TA_8day <- data.table::frollmean(
    a$TA_F,
    n = 8 * n_per_day,
    align = "right",
    na.rm = TRUE
  )
  
  
  #------------------------------------------------------------
  # 3. Preserve the original soil-temperature information
  #
  # TS_repair_flag:
  #   0 = original soil temperature
  #   1 = modeled soil temperature
  #------------------------------------------------------------
  
  a$TS_F_MDS_1_original <- a$TS_F_MDS_1
  a$TS_F_MDS_1_QC_original <- a$TS_F_MDS_1_QC
  a$TS_repair_flag <- 0L
  
  
  #------------------------------------------------------------
  # 4. Identify all known problematic periods
  #
  # These observations must not be used to train any model.
  #------------------------------------------------------------
  
  all_bad <- rep(FALSE, nrow(a))
  
  for (i in seq_len(nrow(bad_periods))) {
    
    all_bad <- all_bad |
      (
        a$YEAR == bad_periods$YEAR[i] &
          a$DOY >= bad_periods$DOY_start[i] &
          a$DOY <= bad_periods$DOY_end[i]
      )
  }
  
  
  #------------------------------------------------------------
  # 5. Model each problematic period separately
  #------------------------------------------------------------
  
  for (i in seq_len(nrow(bad_periods))) {
    
    bad_year  <- bad_periods$YEAR[i]
    bad_start <- bad_periods$DOY_start[i]
    bad_end   <- bad_periods$DOY_end[i]
    
    
    #----------------------------------------------------------
    # Select the training data
    #
    # Use the same seasonal period from all years, extending
    # the interval by 30 days before and after.
    #
    # For the problematic year, the actual bad period is
    # removed by !all_bad. The clean 30-day periods before and
    # after remain available for training.
    #
    # No fixed DOY limits are used, so shifted Southern
    # Hemisphere DOY values are supported.
    #----------------------------------------------------------
    
    train_index <-
      a$DOY >= bad_start - 30 &
      a$DOY <= bad_end + 30 &
      !all_bad &
      a$TS_F_MDS_1_QC %in% c(0, 1) &
      complete.cases(
        a$TS_F_MDS_1,
        a$TA_F,
        a$TA_24h,
        a$TA_8day
      )
    
    
    # Stop if there are too few observations to fit the model
    if (sum(train_index) < 100) {
      stop(
        "Not enough good training observations for ",
        bad_year, ", DOY ", bad_start, "-", bad_end
      )
    }
    
    
    #----------------------------------------------------------
    # Fit the linear soil-temperature model
    #----------------------------------------------------------
    
    TS_model <- lm(
      TS_F_MDS_1 ~ TA_F + TA_24h + TA_8day,
      data = a[train_index, ]
    )
    
    
    #----------------------------------------------------------
    # Select observations inside the problematic period
    #
    # Air-temperature predictors must be available.
    #----------------------------------------------------------
    
    predict_index <-
      a$YEAR == bad_year &
      a$DOY >= bad_start &
      a$DOY <= bad_end &
      complete.cases(
        a$TA_F,
        a$TA_24h,
        a$TA_8day
      )
    
    
    #----------------------------------------------------------
    # Predict and replace soil temperature only during the
    # problematic period
    #----------------------------------------------------------
    
    a$TS_F_MDS_1[predict_index] <- predict(
      TS_model,
      newdata = a[predict_index, ]
    )
    
    
    #----------------------------------------------------------
    # Mark the modeled soil-temperature observations
    #
    # QC = 2 identifies these values as modeled rather than
    # original measured or good-quality observations.
    #----------------------------------------------------------
    
    a$TS_F_MDS_1_QC[predict_index] <- 1
    a$TS_repair_flag[predict_index] <- 1L
    
    
    # Print a short summary
    message(
      "Repaired TS for ", bad_year,
      ", DOY ", bad_start, "-", bad_end,
      "; training observations = ", sum(train_index),
      "; predicted observations = ", sum(predict_index)
    )
  }
  
  
  # Print the detected temporal resolution
  message(
    "Data interval = ", interval_hours,
    " hours; observations per day = ", n_per_day
  )
  
  return(a)
}