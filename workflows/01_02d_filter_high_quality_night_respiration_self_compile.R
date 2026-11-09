# Purpose: prepare standard ac.csv and nightNEE.csv files for self-compiled
# flux sites that do not have standard FLUXNET/AmeriFlux products.
#
# The input NEE products are already quality controlled, u* filtered, and
# gap-filled by their data providers. Therefore, this script does not rerun
# REddyProc or estimate u* thresholds. It retains the remaining workflow from
# 01_02b: fixed local time, sunrise/sunset classification, TS reconstruction,
# growing-season features, good-year gap screening, climatological filling,
# and standard output formatting.

library(data.table)
library(dplyr)
library(lubridate)
library(suncalc)

rm(list = ls())

dir_rawdata <- "/Volumes/WZZ_disk/Thermal_Acclimation"
additional_dir <- file.path(dir_rawdata, "additional_data")
output_dir <- file.path(dir_rawdata, "RespirationData")
graph_dir <- "graphs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(graph_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("workflows", "00_functions.R"))
site_info <- read.csv(
  file.path("data", "site_info.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

site_files <- c(
  "BR-ATT" = file.path(additional_dir, "ATTO_617_with_TA_TS_SWC_local_time.csv"),
  "BR-Tgf" = file.path(additional_dir, "EC_Fire_2014-2024_GAPFILLED_Final.csv"),
  "BR-Tgc" = file.path(additional_dir, "EC_Control_2014-2024_GAPFILLED_Final.csv")
)

missing_files <- site_files[!file.exists(site_files)]
if (length(missing_files)) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}

# Use the existing reconstruct_TS() configuration from 00_functions.R.
ts_config_self <- ts_config
missing_ts_config <- setdiff(names(site_files), ts_config_self$site)
if (length(missing_ts_config)) {
  stop("Missing TS configuration for: ", paste(missing_ts_config, collapse = ", "))
}

as_numeric <- function(x) suppressWarnings(as.numeric(x))

# Fill residual gaps from the mean seasonal-diurnal cycle. This is the same
# basic approach used near the end of 01_02b, with an hour-of-day fallback.
fill_climatology <- function(x, doy, hour, minute) {
  out <- as.numeric(x)
  if (!anyNA(out)) return(out)

  d <- data.table(value = out, DOY = doy, HOUR = hour, MINUTE = minute)
  clim <- d[, .(clim = mean(value, na.rm = TRUE)), by = .(DOY, HOUR, MINUTE)]
  clim[!is.finite(clim), clim := NA_real_]
  key <- paste(d$DOY, d$HOUR, d$MINUTE)
  clim_key <- paste(clim$DOY, clim$HOUR, clim$MINUTE)
  replacement <- clim$clim[match(key, clim_key)]
  out[is.na(out)] <- replacement[is.na(out)]

  if (anyNA(out)) {
    d[, value := out]
    hour_clim <- d[, .(clim = mean(value, na.rm = TRUE)), by = .(HOUR, MINUTE)]
    hour_clim[!is.finite(clim), clim := NA_real_]
    key <- paste(d$HOUR, d$MINUTE)
    hour_key <- paste(hour_clim$HOUR, hour_clim$MINUTE)
    replacement <- hour_clim$clim[match(key, hour_key)]
    out[is.na(out)] <- replacement[is.na(out)]
  }

  if (anyNA(out)) {
    fallback <- mean(out, na.rm = TRUE)
    if (is.finite(fallback)) out[is.na(out)] <- fallback
  }
  out
}

parse_year_removed <- function(x) {
  if (length(x) == 0L || is.na(x) || trimws(x) == "") return(integer())
  parts <- trimws(strsplit(x, ",")[[1]])
  unique(unlist(lapply(parts, function(part) {
    if (grepl(":", part, fixed = TRUE)) {
      limits <- as.integer(strsplit(part, ":", fixed = TRUE)[[1]])
      seq(limits[1], limits[2])
    } else {
      as.integer(part)
    }
  })))
}

# Convert an adjusted Southern-Hemisphere DOY back to the correct local date.
adjusted_doy_timestamp <- function(growing_year, adjusted_doy,
                                   southern_hemisphere, tz_site,
                                   half_interval_seconds) {
  if (southern_hemisphere && adjusted_doy > 366) {
    calendar_year <- growing_year + 1L
    calendar_doy <- adjusted_doy - 366
  } else {
    calendar_year <- growing_year
    calendar_doy <- adjusted_doy
  }
  as.POSIXct(
    as.Date(sprintf("%04d-01-01", calendar_year)) + calendar_doy - 1,
    tz = tz_site
  ) + half_interval_seconds
}

read_self_compiled_site <- function(name_site, file_path) {
  if (name_site == "BR-ATT") {
    timestamp_col <- "Local Time (start time each 30 min)"
    x <- fread(file_path, colClasses = setNames("character", timestamp_col))
    timestamp <- as.POSIXct(x[[timestamp_col]], format = "%Y-%m-%d %H:%M:%S",
                            tz = "Etc/GMT+4")

    ta <- as_numeric(x$TA_F)
    ta_measured <- as_numeric(x$TA_26m)
    ta_qc <- ifelse(is.na(ta), NA_real_, ifelse(!is.na(ta_measured), 0, 1))

    data.table(
      TIMESTAMP = timestamp,
      NEE_measured = as_numeric(x$NEE_nf_umol_m2_s),
      NEE_continuous = as_numeric(x$NEE_f_umol_m2_s),
      TA_F_MDS = ta,
      TA_F_MDS_QC = ta_qc,
      TS_F_MDS_1 = as_numeric(x$TS_5cm),
      TS_F_MDS_1_QC = ifelse(!is.na(as_numeric(x$TS_5cm)), 0, NA_real_),
      SWC = as_numeric(x$SWC_5cm),
      SW_IN = NA_real_
    )
  } else {
    x <- fread(file_path)
    timestamp_utc <- as.POSIXct(
      x$DateTime,
      format = "%Y-%m-%dT%H:%M:%OSZ",
      tz = "UTC"
    )
    # Convert UTC to a fixed UTC-4 local clock without daylight-saving changes.
    timestamp_local <- as.POSIXct(
      format(timestamp_utc, tz = "Etc/GMT+4", format = "%Y-%m-%d %H:%M:%S"),
      format = "%Y-%m-%d %H:%M:%S",
      tz = "Etc/GMT+4"
    )

    ta <- as_numeric(x$Tair_f)
    ta_source_qc <- as_numeric(x$Tair_fqc)
    ta_qc <- ifelse(is.na(ta), NA_real_, ifelse(ta_source_qc == 0, 0, 1))
    ts_measured <- as_numeric(x$Tsoil)

    data.table(
      TIMESTAMP = timestamp_local,
      NEE_measured = as_numeric(x$NEE_AfterU_orig),
      NEE_continuous = as_numeric(x$NEE_AfterU_f),
      TA_F_MDS = ta,
      TA_F_MDS_QC = ta_qc,
      TS_F_MDS_1 = ts_measured,
      TS_F_MDS_1_QC = ifelse(!is.na(ts_measured), 0, NA_real_),
      SWC = as_numeric(x$SWC_1_1_1) / 100,
      SW_IN = as_numeric(x$Rg_f)
    )
  }
}

save_site_diagnostic_plot <- function(name_site, measured_night,
                                      NEE_yearly, NEE_yearly_night,
                                      NEE_yearly_day, gStart, gEnd) {
  png(
    file.path(graph_dir, paste0(name_site, "TS_NEEnight_relation.png")),
    width = 1200,
    height = 1200
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(3, 2))

  plot(measured_night$TS, measured_night$NEE,
       main = name_site, xlab = "Soil temperature", ylab = "Nighttime NEE")
  plot(NEE_yearly$DOY, NEE_yearly$NEE,
       xlab = "DOY", ylab = "Mean NEE", main = "All observations")
  mark_period(NEE_yearly$DOY, NEE_yearly$NEE, gStart, gEnd)
  plot(NEE_yearly$DOY, NEE_yearly$TS,
       xlab = "DOY", ylab = "Mean soil temperature", main = "All observations")
  mark_period(NEE_yearly$DOY, NEE_yearly$TS, gStart, gEnd)
  plot(NEE_yearly_night$DOY, NEE_yearly_night$NEE,
       xlab = "DOY", ylab = "Mean nighttime NEE", main = "Nighttime")
  mark_period(NEE_yearly_night$DOY, NEE_yearly_night$NEE, gStart, gEnd)
  plot(NEE_yearly_day$DOY, NEE_yearly_day$TS,
       xlab = "DOY", ylab = "Mean daytime soil temperature", main = "Daytime")
  mark_period(NEE_yearly_day$DOY, NEE_yearly_day$TS, gStart, gEnd)
  plot(NEE_yearly_day$DOY, NEE_yearly_day$NEE,
       xlab = "DOY", ylab = "Mean daytime NEE", main = "Daytime")
  mark_period(NEE_yearly_day$DOY, NEE_yearly_day$NEE, gStart, gEnd)
}

feature_gs <- data.frame(
  site_ID = character(), gStart = double(), gEnd = double(),
  taStart = double(), taEnd = double(), tsStart = double(), tsEnd = double(),
  nyear = integer()
)
gap_diagnostics <- list()
model_diagnostics <- list()
processing_errors <- list()

for (name_site in names(site_files)) {
  message("Processing ", name_site)

  tryCatch({
    id <- match(name_site, site_info$site_ID)
    if (is.na(id)) stop(name_site, " is not present in site_info.csv")

    a <- read_self_compiled_site(name_site, site_files[[name_site]])
    setorder(a, TIMESTAMP)
    if (anyNA(a$TIMESTAMP)) stop(name_site, " has invalid timestamps")
    if (anyDuplicated(a$TIMESTAMP)) stop(name_site, " has duplicate timestamps")

    dt_hours <- median(as.numeric(diff(a$TIMESTAMP), units = "hours"), na.rm = TRUE)
    if (!is.finite(dt_hours) || abs(dt_hours - 0.5) > 1e-6) {
      stop(name_site, " does not have a regular half-hour interval")
    }

    tz_site <- "Etc/GMT+4"
    a[, `:=`(
      YEAR = year(TIMESTAMP),
      MONTH = month(TIMESTAMP),
      DAY = day(TIMESTAMP),
      DOY = yday(TIMESTAMP),
      HOUR = hour(TIMESTAMP),
      MINUTE = minute(TIMESTAMP),
      DATE = as.Date(TIMESTAMP, tz = tz_site)
    )]

    # Updated 01_02b rule: identify Southern-Hemisphere sites by latitude.
    southern_hemisphere <- site_info$LAT[id] < -10
    if (southern_hemisphere) {
      a[DOY < 183, DOY := DOY + 366]
    }
    a[, growing_year := ifelse(DOY <= 366, YEAR, YEAR - 1L)]

    # Calculate sunrise and sunset using fixed local standard time.
    sunlight <- as.data.table(getSunlightTimes(
      date = seq(min(a$DATE), max(a$DATE), by = "day"),
      keep = c("sunrise", "sunset"),
      lat = site_info$LAT[id],
      lon = site_info$LONG[id],
      tz = tz_site
    ))
    a[, sunrise := sunlight$sunrise[match(DATE, sunlight$date)]]
    a[, sunset := sunlight$sunset[match(DATE, sunlight$date)]]
    half_interval_seconds <- dt_hours * 3600 / 2
    a[, daytime := TIMESTAMP + half_interval_seconds >= sunrise &
                      TIMESTAMP - half_interval_seconds <= sunset]
    if (anyNA(a$daytime)) {
      stop(name_site, " has undefined daytime values")
    }

    # Fill any remaining TA gaps before TS reconstruction. Preserve QC=0 only
    # for originally measured TA; all newly filled records receive QC=1.
    missing_ta <- is.na(a$TA_F_MDS)
    a[, TA_F_MDS := fill_climatology(TA_F_MDS, DOY, HOUR, MINUTE)]
    a[missing_ta & !is.na(TA_F_MDS), TA_F_MDS_QC := 1]
    if (anyNA(a$TA_F_MDS)) stop(name_site, " still has missing TA after filling")

    ts_result <- reconstruct_TS(
      a = as.data.frame(a),
      name_site = name_site,
      ts_config = ts_config_self,
      long_mean_days = 8
    )
    a <- as.data.table(ts_result$data)
    a[, TS := TS_F_MDS_1]
    a[, TA := TA_F_MDS]
    a[, NEE := NEE_measured]

    # Interactive checks retained from 01_02c. These appear in the RStudio
    # Plots pane when the script is sourced interactively.
    if (interactive()) {
      plot(a$TA, type = "l", main = paste(name_site, "air temperature"),
           xlab = "Record", ylab = "TA")
      plot(a$TS, type = "l", main = paste(name_site, "soil temperature"),
           xlab = "Record", ylab = "TS")
      if (site_info$SWC_use[id] == "YES") {
        plot(a$SWC, type = "l", main = paste(name_site, "soil water content"),
             xlab = "Record", ylab = "SWC")
      }
      plot(a$NEE, type = "l", main = paste(name_site, "measured NEE"),
           xlab = "Record", ylab = "NEE")
    }

    NEE_yearly <- a[, .(
      NEE = mean(NEE, na.rm = TRUE),
      TS = mean(TS, na.rm = TRUE),
      TA = mean(TA, na.rm = TRUE)
    ), by = DOY]
    NEE_yearly_night <- a[daytime == FALSE, .(
      NEE = mean(NEE, na.rm = TRUE),
      TS = mean(TS, na.rm = TRUE),
      TA = mean(TA, na.rm = TRUE)
    ), by = DOY]
    NEE_yearly_day <- a[daytime == TRUE, .(
      NEE = mean(NEE, na.rm = TRUE),
      TS = mean(TS, na.rm = TRUE),
      TA = mean(TA, na.rm = TRUE)
    ), by = DOY]

    model_summary <- summary(ts_result$model)
    model_diagnostics[[name_site]] <- data.frame(
      site_ID = name_site,
      training_n = ts_result$training_n,
      r_squared = model_summary$r.squared,
      adjusted_r_squared = model_summary$adj.r.squared,
      rmse = sqrt(mean(residuals(ts_result$model)^2)),
      stringsAsFactors = FALSE
    )

    # Tropical-site growing-season logic retained from 01_02b.
    if (abs(site_info$LAT[id]) < 23.5 ||
        (site_info$Climate_class[id] %in% c("Af", "Am", "Aw", "As", "Cfa", "Cfb", "Cwa") &&
         site_info$IGBP[id] == "EBF")) {
      gStart <- min(a[daytime == FALSE, DOY], na.rm = TRUE)
      gEnd <- max(a[daytime == FALSE, DOY], na.rm = TRUE)
    } else {
      ac_tmp <- as.data.frame(a)
      ac_tmp$NIGHT <- !ac_tmp$daytime
      ac_tmp$NEE_VUT_REF <- ac_tmp$NEE
      ac_tmp$NEE_VUT_REF_QC <- as.numeric(is.na(ac_tmp$NEE))
      er_period <- estimate_active_ER_period(
        a = ac_tmp,
        site_id = name_site,
        ta_var = "TA",
        ts_var = "TS",
        southern_hemisphere = southern_hemisphere,
        threshold_fraction = 0.20,
        temp_threshold = 2,
        min_observed_doy = 180L
      )
      if (er_period$status != "ok") stop(name_site, " growing season was skipped")
      gStart <- er_period$gStart
      gEnd <- er_period$gEnd
    }

    if (!is.na(site_info$gStart[id])) gStart <- as.numeric(site_info$gStart[id])
    if (!is.na(site_info$gEnd[id])) gEnd <- as.numeric(site_info$gEnd[id])

    growing_data <- a[DOY >= gStart & DOY <= gEnd]
    tsStart <- quantile(growing_data$TS, 0.025, na.rm = TRUE)
    tsEnd <- quantile(growing_data$TS, 0.975, na.rm = TRUE)
    taStart <- quantile(growing_data$TA, 0.025, na.rm = TRUE)
    taEnd <- quantile(growing_data$TA, 0.975, na.rm = TRUE)

    # Measured, provider-filtered nighttime NEE. SWC is not required because
    # all three sites have SWC_use == "NO" in site_info.csv.
    measured_night <- a[
      daytime == FALSE &
        !is.na(NEE) & !is.na(TA) & !is.na(TS) &
        NEE > -5 & NEE < 30
    ]
    if (!nrow(measured_night)) stop(name_site, " has no usable nighttime NEE")

    # Match the two annual diagnostic plots displayed by 01_02c.
    if (interactive()) {
      for (iyear in sort(unique(measured_night$YEAR))) {
        data_year <- measured_night[YEAR == iyear]
        if (nrow(data_year) > 100) {
          plot(data_year$TS, data_year$NEE, main = iyear,
               xlab = "Soil temperature", ylab = "Nighttime NEE")
          plot(data_year$TIMESTAMP, data_year$NEE, main = iyear,
               xlim = as.POSIXct(c(
                 paste0(iyear, "-01-01 00:00:00"),
                 paste0(iyear, "-12-31 23:59:59")
               ), tz = tz_site),
               xlab = "Local time", ylab = "Nighttime NEE")
        } else {
          message(name_site, " has <=100 nighttime observations in ", iyear)
        }
      }
    }

    # Always save the same six-panel diagnostic PNG produced by 01_02c.
    save_site_diagnostic_plot(
      name_site = name_site,
      measured_night = measured_night,
      NEE_yearly = NEE_yearly,
      NEE_yearly_night = NEE_yearly_night,
      NEE_yearly_day = NEE_yearly_day,
      gStart = gStart,
      gEnd = gEnd
    )

    # Good/bad growing-year calculation retained from 01_02b.
    season_length <- gEnd - gStart + 1
    gap_max_thresh <- season_length * 0.5
    gap_total_thresh <- 0.7
    all_growing_years <- sort(unique(a$growing_year))

    boundary_rows <- rbindlist(lapply(all_growing_years, function(gy) {
      data.table(
        growing_year = gy,
        TIMESTAMP = c(
          adjusted_doy_timestamp(gy, gStart, southern_hemisphere, tz_site,
                                 half_interval_seconds),
          adjusted_doy_timestamp(gy, gEnd, southern_hemisphere, tz_site,
                                 half_interval_seconds)
        ),
        DOY = c(gStart, gEnd)
      )
    }))

    gap_input <- unique(rbindlist(list(
      measured_night[, .(growing_year, TIMESTAMP, DOY)],
      boundary_rows
    )))
    setorder(gap_input, growing_year, TIMESTAMP)
    gap_result <- gap_input[
      DOY >= gStart & DOY <= gEnd,
      {
        gaps <- as.numeric(diff(TIMESTAMP), units = "days")
        if (!length(gaps)) {
          list(gap_max = Inf, gap_total = Inf)
        } else {
          list(
            gap_max = max(gaps, na.rm = TRUE),
            gap_total = sum(ifelse(gaps > 14, gaps, 0), na.rm = TRUE) /
              season_length
          )
        }
      },
      by = growing_year
    ]
    gap_result[, `:=`(
      site_ID = name_site,
      gap_max_thresh = gap_max_thresh,
      gap_total_thresh = gap_total_thresh,
      accepted = gap_max < gap_max_thresh & gap_total < gap_total_thresh
    )]

    good_years <- gap_result[accepted == TRUE, growing_year]
    good_years <- setdiff(good_years, parse_year_removed(site_info$year_removed[id]))
    gap_result[, accepted := growing_year %in% good_years]
    gap_diagnostics[[name_site]] <- as.data.frame(gap_result)
    if (!length(good_years)) stop(name_site, " has no good growing years")

    measured_night <- measured_night[growing_year %in% good_years]
    ac <- a[growing_year %in% good_years, .(
      YEAR, MONTH, DAY, DOY, HOUR, MINUTE,
      NEE_uStar_f = NEE_continuous,
      TA, TS, SWC, SW_IN, daytime
    )]

    # Fill residual gaps in continuous ac variables using the same seasonal-
    # diurnal approach. Measured nightNEE is never filled.
    ac[, TA := fill_climatology(TA, DOY, HOUR, MINUTE)]
    ac[, TS := fill_climatology(TS, DOY, HOUR, MINUTE)]
    ac[, NEE_uStar_f := fill_climatology(NEE_uStar_f, DOY, HOUR, MINUTE)]

    night_output <- measured_night[, .(
      YEAR, MONTH, DAY, DOY, HOUR, MINUTE,
      NEE, TA, TS, SWC
    )]

    fwrite(ac, file.path(output_dir, paste0(name_site, "_ac.csv")), na = "NA")
    fwrite(
      night_output,
      file.path(output_dir, paste0(name_site, "_nightNEE.csv")),
      na = "NA"
    )

    feature_gs <- bind_rows(
      feature_gs,
      data.frame(
        site_ID = name_site,
        gStart = gStart,
        gEnd = gEnd,
        taStart = as.numeric(taStart),
        taEnd = as.numeric(taEnd),
        tsStart = as.numeric(tsStart),
        tsEnd = as.numeric(tsEnd),
        nyear = length(good_years)
      )
    )

    message(
      name_site, ": wrote ", nrow(ac), " ac rows and ",
      nrow(night_output), " nighttime NEE rows; good years = ",
      paste(good_years, collapse = ", ")
    )
  }, error = function(e) {
    processing_errors[[name_site]] <<- data.frame(
      site_ID = name_site,
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    message(name_site, " failed: ", conditionMessage(e))
  })
}

write.csv(
  feature_gs,
  file.path("data", "growing_season_feature_self_compile.csv"),
  row.names = FALSE
)
# if (length(gap_diagnostics)) {
#   fwrite(
#     rbindlist(gap_diagnostics, fill = TRUE),
#     file.path(output_dir, "self_compile_year_gap_diagnostics.csv")
#   )
# }
# if (length(model_diagnostics)) {
#   fwrite(
#     rbindlist(model_diagnostics, fill = TRUE),
#     file.path(output_dir, "self_compile_TS_model_diagnostics.csv")
#   )
# }
# if (length(processing_errors)) {
#   fwrite(
#     rbindlist(processing_errors, fill = TRUE),
#     file.path(output_dir, "self_compile_processing_errors.csv")
#   )
# }

message("Finished. Outputs saved in: ", output_dir)
