library(librarian)
shelf(dplyr, lubridate)
rm(list=ls())

site_info <- read.csv(file.path('data', 'site_info.csv'))
site_sparkle <- read.csv(file.path('data', 'FLUXNET_missing_sites_TS.csv'))

site_sparkle <- site_sparkle %>% filter(!IGBP %in% c("CRO", "SNO", "BSV", "CVM")) # %>% filter(Checked != "bad")
# 112 sites
# 21 sites has no top soil temperature.

site_shuttle <- read.csv(file.path('data', 'fluxnet_shuttle_snapshot_20260908T163200_selectedsites.csv'))

site_shuttle <- site_shuttle %>% filter(!igbp %in% c("CRO", "SNO", "BSV", "CVM"))   
site_wesort <- unique(c(site_info$site_ID, site_sparkle$Site))

site_shuttle$site_id[!site_shuttle$site_id %in% site_wesort]


# "DE-Zrk" "FI-Ken" "FR-Tou" "IT-Cpz" "IT-Ro1" "IT-Ro2" "SE-St1" "US-Aud" "US-Ho3" "US-KLS" "US-NR3"
# we have downloaded all the sites. 

# ICOS's data: "IT-Cp2", IT-Ro1, IT-Ro2

# not enough: SE-St1; "US-Ho3" harvested, "US-NR3"

# why i did not use this one:  US-Aud; US-KLS; "DE-Zrk"; US-KFS(too bad data)?

# explanation of why excluding sites with >= 8-year data and satisfy our site selection criteria
# AU-Wrr	A wildfire burned through the tower site in January 2019, within its 2013–2021 record.
# CA-Ca2	Douglas-fir stand clear-cut in winter 1999/2000; measurements began in 2000.
# FI-Let	Partially harvested in March 2016, within its measurement record.
# IT-Ro1	Oak forest coppiced in December 1999; its 2000–2008 record follows early regrowth.
# US-Aud	Burned in 2002 and reported to be recovering through late 2005.
# US-NC1	A 70-ha forest clear-cut in 2004, then replanted; measurements began in 2005.
# US-Rpf	Burned black-spruce forest undergoing post-fire succession during its 2008–2024 record.
# US-UMd	Experimental girdling in 2008 killed a substantial fraction of the canopy and initiated succession.
# US-Tw4  cropland restored to wetland at the start of monitoring. 
# CA-TP1  plantation established in 2002, with stand development measured from its first year. 
# US-Me6  post-fire pine stand, but the fire preceded its 2010 record by many years
# IT-Lav  site's ER has known problems. 
# IT-Tor  site's ER has known problems. Galvagno et al. (2017), Contribution of advection to nighttime ecosystem respiration at a mountain grassland in complex terrain. Agricultural and Forest Meteorology, 237–238, 270–281.
# IT-TrF  nighttime NEEs have problems. seems not balanced.
# DE-Zrk  previously drained peatland rewetted in 2004–2005, before its 2013–2022 record; data has clear trend
# CN-Aro  data are super noisy. 
# FR-LGt  bad data
# CA-DBB  I removed this site because of its negative NEE during summer nighttime (abnormal). 
# CN-Ash  I removed this site because of the same reason. 


#-----------------------------------------
# US-esm  weird higly fluctuating data
# US-elm  weird higly fluctuating data


# sites we have to exclude years
# CA-DBB: 

# check these sites
# AT-Neu: this site is a managed agricultural and grassland. NEE is highly positive, distinct from other sites; it may have problems


# unzip files
# Extract only the FLUXMET HH, DD, and YY CSV files from a FLUXNET zip archive.

#----------------------------------------------
# Unzip selected FLUXNET files:
# For each zip, extract only FLUXMET HH, DD, and YY CSV files.

site_dir <- "/Volumes/WZZ_disk/Thermal_Acclimation/SiteData/FLUXNET_Shuttle09102026"
out_dir <- file.path(site_dir, "unzip")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zip_files <- list.files(
  site_dir,
  pattern = "\\.zip$",
  full.names = TRUE
)

for (zip_file in zip_files) {
  zip_contents <- unzip(zip_file, list = TRUE)
  
  files_to_extract <- zip_contents$Name[
    grepl("_FLUXMET_(HH|DD|YY)_.*\\.csv$", zip_contents$Name)
  ]
  
  if (length(files_to_extract) != 3) {
    warning(
      "Skipping ", basename(zip_file),
      ": expected 3 FLUXMET files (HH, DD, YY), found ",
      length(files_to_extract)
    )
    next
  }
  
  output_files <- file.path(out_dir, basename(files_to_extract))
  
  missing_files <- files_to_extract[!file.exists(output_files)]
  
  if (length(missing_files) == 0) {
    message("Already complete, skipping: ", basename(zip_file))
    next
  }
  
  unzip(
    zipfile = zip_file,
    files = missing_files,
    exdir = out_dir
  )
  
  message(
    "Extracted ", length(missing_files),
    " missing file(s) from: ", basename(zip_file)
  )
}

message("Done. Files saved to: ", out_dir)

# file processing is done!
#################################

TS_scrren_table <- read.csv("/Users/junnawang/Downloads/soil-temperature-screentable.csv")
bad_sites <- TS_scrren_table$Site[TS_scrren_table$verdict == "BAD"] 


