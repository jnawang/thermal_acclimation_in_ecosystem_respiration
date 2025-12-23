# This script obtains environmental and biological conditions potentially affecting thermal response strength
# Authors: Junna Wang, October, 2025

# This script gets soil carbon, climate conditions, and spectral information (GPP, LAI, EVI, NDVI). 
# Spectral information is from NASA earth data. 
# Elevation of DE-SfS is from 5-m global elevation map and elevation of other sites are from PI.

# this script takes 5 mins to run. 

library(librarian)
shelf(dplyr, lubridate, ggpubr, ggplot2, zoo, terra, corrplot)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
dir_rawdata <- '/Users/junnawang/YaleLab/data_server/'
####################End Attention

site_info <- read.csv(file.path('data', 'site_info.csv'))

#--------------------------------------------SOIL DATA-------------------------------------
# get measured soil carbon data from AmeriFlux BIF data
BIF <- read.csv(file.path(dir_rawdata, 'AMF_AA-Net_BIF_CCBY4_20251017.csv'))
BIF.site <- BIF %>% filter(SITE_ID %in% site_info$site_ID) %>% filter(VARIABLE=='SOIL_CHEM_C_ORG') %>% group_by(SITE_ID) %>% summarise(soc_obs = mean(as.numeric(DATAVALUE), na.rm=T), .groups = 'drop')     # 26 sites
# the unit of SOIL_CHEM_C_ORG: g C kg soil-1

# if no measured data available, we use GSOCmap1.5.0.tif
# the unit of this product: total SOC stocks t/ha
stat.soil <- data.frame(site_ID=site_info$site_ID)
#
GSOCmap <- terra::rast(file.path(dir_rawdata, 'GSOCmap1.5.0.tif'))
# plot(GSOCmap)
#
xy <- data.frame(x=site_info$LONG, y=site_info$LAT)
#
stat.soil$GSOC <- terra::extract(GSOCmap, xy)$GSOCmap1.5.0

# check how GSOC data are good
# Compare these field data with map data, why I did not think of this method? # 27 sites
tmp <- BIF.site %>% left_join(stat.soil, by=c('SITE_ID'='site_ID'))
cor.test(tmp$soc_obs, tmp$GSOC)  # p-value = 0.02817, r = 0.4224124
# plot(tmp$soc_obs, tmp$GSOC)

# I need to do some corrections for GSOC
model <- lm(data=tmp, GSOC ~ soc_obs)
BIF.site$GSOC_obs <- predict(model, tmp[c("soc_obs")])

stat.soil <- stat.soil %>% left_join(BIF.site, by=c('site_ID'='SITE_ID'))
stat.soil <- stat.soil %>% mutate(SOC = case_when(is.na(GSOC_obs) ~ GSOC,
                                                  !is.na(GSOC_obs) ~ GSOC_obs))

#---------------------------------------------------CLIMATE DATA---------------------------------
# calculate the following climatic variables
# mean annual NEE, mean annual GPP, and mean annual night NEE [I have to differentiate daytime and nighttime].
stat.climate <- data.frame(site_ID=character(), NEE=double(), NEE_day=double(), NEE_night=double(),                      
                           MATA=double(), SSTA=double(), IATA=double(), DRTA=double(), RSTA=double(),
                           MATS=double(), SSTS=double(), IATS=double(), DRTS=double(), RSTS=double(),
                           NEE2=double(), NEE_day2=double(), NEE_night2=double(),
                           NEE_FLUXNET=double(), MAT_FLUXNET=double(), P_FLUXNET=double())
# notation of variable names
# MA: mean annual; SS: seasonal variation; IA: inter-annual variation; RD: daily range; RS: seasonal range
#
# Read ac data
files <- list.files(file.path(dir_rawdata, 'RespirationData'), pattern = '_ac.csv$', full.names = TRUE)
for (i in 1:length(files)) {
  # i = 67
  name_site <- substring(files[i], nchar(files[i])-12, nchar(files[i])-7)
  print(paste0(i, name_site))
  stat.climate[i, 1] <- name_site        # the last value
  ac <- read.csv(files[i])
  
  # use only the years with qualified data
  a_measure_night_complete <- read.csv(file.path(dir_rawdata, 'RespirationData', paste0(name_site, '_nightNEE.csv')))
  good_years <- unique(a_measure_night_complete$YEAR)
  
  # I have to gap fill TA, TS, and NEE, if needed.
  if (!"TA_gf" %in% colnames(ac)) {
    T_gf <- ac %>% group_by(DOY, HOUR, MINUTE) %>% summarise(TA_gf=mean(TA, na.rm=T), TS_gf=mean(TS, na.rm=T), NEE_gf=mean(NEE, na.rm=T), .groups = 'drop')  # for gap fill only.    
    ac   <- ac %>% left_join(T_gf, by=c("DOY", "HOUR", "MINUTE"))
    #
    # I still need a standard gap-fill methods. 
    if (sum(is.na(ac$TA)) > 0) {
      ac$TA[is.na(ac$TA)] <- ac$TA_gf[is.na(ac$TA)]
    }
    if (sum(is.na(ac$TS)) > 0) {
      ac$TS[is.na(ac$TS)] <- ac$TS_gf[is.na(ac$TS)]
    }
    if (sum(is.na(ac$NEE)) > 0) {
      ac$NEE[is.na(ac$NEE)] <- ac$NEE_gf[is.na(ac$NEE)]
    } 
  }
  #
  # Do NEE separately
  # only use the selected years data; find good years from other datasets!
  # use interpolated values
  annual <- ac %>% filter(YEAR %in% good_years) %>% group_by(YEAR) %>% summarise(NEE=sum(NEE_uStar_f, na.rm=T), .groups = 'drop')
  annual_day   <- ac %>% filter(YEAR %in% good_years & daytime) %>% group_by(YEAR) %>% summarise(NEE=sum(NEE_uStar_f, na.rm=T), .groups = 'drop')
  annual_night <- ac %>% filter(YEAR %in% good_years & !daytime) %>% group_by(YEAR) %>% summarise(NEE=sum(NEE_uStar_f, na.rm=T), .groups = 'drop')
  # get time interval
  # dt <- as.numeric(ac$TIMESTAMP[2] - ac$TIMESTAMP[1])     this is a wrong way!
  dt <- abs(ac$MINUTE[2] + ac$HOUR[2] * 60 - ac$MINUTE[1] - ac$HOUR[1] * 60)
  #  
  ######alternative methods
  stat.climate[i, 2] <- mean(annual$NEE) * dt * 60 / 1000000 * 12        ### g C / m2;
  stat.climate[i, 3] <- mean(annual_day$NEE) * dt * 60 / 1000000 * 12    ### g C / m2;
  stat.climate[i, 4] <- mean(annual_night$NEE) * dt * 60 / 1000000 * 12  ### g C / m2;
  #
  # I have to interpolate missing values
  data_yearly  <- ac %>% filter(YEAR %in% good_years) %>% group_by(DOY, HOUR, MINUTE, daytime) %>% summarise(NEE=mean(NEE_uStar_f, na.rm=T), .groups = 'drop')
  # remove duplicate rows
  data_yearly  <- data_yearly[!duplicated(data_yearly[,c('DOY', 'HOUR', 'MINUTE')]), ]
  #
  data_yearly$time <- data_yearly$DOY + data_yearly$HOUR / 24 + data_yearly$MINUTE / 60 / 24
  # re-order the data
  data_yearly <- data_yearly %>% arrange(time)
  # then linear interpolation
  data_yearly$NEE.appr <- na.approx(data_yearly$NEE, data_yearly$time, na.rm = FALSE)
  # dealing with missing head and end data
  num_min <- min(which(!is.na(data_yearly$NEE.appr)))
  num_max <- max(which(!is.na(data_yearly$NEE.appr)))
  if (num_min > 1) {
    data_yearly$NEE.appr[1:(num_min-1)] <- data_yearly$NEE.appr[num_min]
  }
  if (num_max < nrow(data_yearly)) {
    data_yearly$NEE.appr[(num_max+1):nrow(data_yearly)] <- data_yearly$NEE.appr[num_max]
  } 
  # then sum
  stat.climate[i, 15] <- sum(data_yearly$NEE.appr) * dt * 60 / 1000000 * 12  ### g CO2 / m2  
  stat.climate[i, 16] <- sum(data_yearly$NEE.appr[data_yearly$daytime], na.rm=T) * dt * 60 / 1000000 * 12   ### g CO2 / m2; na.rm for tundra sites
  stat.climate[i, 17] <- sum(data_yearly$NEE.appr[!data_yearly$daytime], na.rm=T) * dt * 60 / 1000000 * 12  ### g CO2 / m2
  ##############################air temperature and soil temperature############################
  data_yearly  <- ac %>% filter(YEAR %in% good_years) %>% group_by(DOY) %>% summarise(TS=mean(TS, na.rm=T), TA=mean(TA, na.rm=T), .groups = 'drop')
  # Mean annual
  stat.climate[i, 5] <- mean(data_yearly$TA)
  stat.climate[i, 10] <- mean(data_yearly$TS)
  #
  # Seasonal variations
  data_1day      <- ac %>% filter(YEAR %in% good_years) %>% group_by(YEAR, DOY) %>% summarise(TSrange=max(TS, na.rm=T)-min(TS, na.rm=T), TArange=max(TA, na.rm=T)-min(TA, na.rm=T), TS=mean(TS, na.rm=T), TA=mean(TA, na.rm=T), .groups = 'drop')
  data_1year_sd  <- data_1day %>% group_by(YEAR) %>% summarise(TS=sd(TS, na.rm=T), TA=sd(TA, na.rm=T), .groups = 'drop')
  stat.climate[i, 6]  <- mean(data_1year_sd$TA)
  stat.climate[i, 11] <- mean(data_1year_sd$TS)
  #
  # Inter-annual variations
  data_1year_mn  <- data_1day %>% group_by(YEAR) %>% summarise(TS=mean(TS, na.rm=T), TA=mean(TA, na.rm=T), .groups = 'drop')
  stat.climate[i, 7]  <- sd(data_1year_mn$TA)
  stat.climate[i, 12] <- sd(data_1year_mn$TS)
  # Daily range
  stat.climate[i, 8]   <- mean(data_1day$TArange)
  stat.climate[i, 13]  <- mean(data_1day$TSrange, na.rm=T)
  # Seasonal range
  data_1year_rg  <- data_1day %>% group_by(YEAR) %>% summarise(TS=max(TS, na.rm=T)-min(TS, na.rm=T), TA=max(TA, na.rm=T)-min(TA, na.rm=T), .groups = 'drop')
  stat.climate[i, 9]  <- mean(data_1year_rg$TA)
  stat.climate[i, 14] <- mean(data_1year_rg$TS)
}

# correlation of these climatic variables?
cor(stat.climate[, 2:14])

#---------------------------------------------------SPECTRAL DATA------------------------------------------------
# Take a look at spectral data
file1    <- read.csv(file.path(dir_rawdata, 'towers-MOD13A2-061-results.csv'))  ##EVI and NDVI

# only use good quality data
file1$MOD13A2_061__1_km_16_days_EVI[!file1$MOD13A2_061__1_km_16_days_VI_Quality_MODLAND_Description %in% c("VI produced, good quality", "VI produced, but check other QA")]  <- NA
file1$MOD13A2_061__1_km_16_days_NDVI[!file1$MOD13A2_061__1_km_16_days_VI_Quality_MODLAND_Description %in% c("VI produced, good quality", "VI produced, but check other QA")]  <- NA
file1$Date <- as.Date(file1$Date, format = "%m/%d/%y")

# Alternative way to calculate annual average values
file1_interp <- file1 %>% mutate(month = month(Date)) %>% group_by(ID, month) %>% 
  summarise(NDVI = max(mean(MOD13A2_061__1_km_16_days_NDVI, na.rm=T), 0.0), n_NDVI = sum(!is.na(MOD13A2_061__1_km_16_days_NDVI)), 
            EVI = max(mean(MOD13A2_061__1_km_16_days_EVI, na.rm=T), 0.0), n_EVI = sum(!is.na(MOD13A2_061__1_km_16_days_EVI)), .groups = 'drop')
# interpolate months with missed values
file1_interp <- file1_interp %>% group_by(ID) %>% mutate(NDVI_interp = na.approx(NDVI, rule = 2), EVI_interp = na.approx(EVI, rule = 2))


# file1_interp <- data.frame()
# IDs <- unique(file1$ID)
# for (id in IDs) {
#   # print(id)
#   data.site <- file1 %>% filter(ID == id)
#   data.site$NDVI <- na.approx(data.site$MOD13A2_061__1_km_16_days_NDVI, data.site$Date, na.rm = FALSE)
#   data.site$EVI <- na.approx(data.site$MOD13A2_061__1_km_16_days_EVI, data.site$Date, na.rm = FALSE)
#   file1_interp <- rbind(file1_interp, data.site[, c('ID', 'Date', 'NDVI', 'EVI')])
#   #
#   # p <- ggplot(data = data.site, aes(x=Date, y=MOD13A2_061__1_km_16_days_NDVI)) +
#   #   geom_point() +
#   #   labs(title = id)
#   # plot(p)
# }

#######
file2    <- read.csv(file.path(dir_rawdata, 'towers-MOD15A2H-061-results.csv'))  ## Fpar and LAI

# remove low quality data
file2$MOD15A2H_061_Fpar_500m[file2$MOD15A2H_061_FparLai_QC_MODLAND_Description == "Other Quality (back-up algorithm or fill values)"] <- NA
file2$MOD15A2H_061_Lai_500m[file2$MOD15A2H_061_FparLai_QC_MODLAND_Description == "Other Quality (back-up algorithm or fill values)"] <- NA
file2$Date <- as.Date(file2$Date, format = "%m/%d/%y")

# Alternative way to calculate annual average values
file2_interp <- file2 %>% mutate(month = month(Date)) %>% group_by(ID, month) %>% 
  summarise(Fpar = max(mean(MOD15A2H_061_Fpar_500m, na.rm=T), 0.0), n_Fpar = sum(!is.na(MOD15A2H_061_Fpar_500m)), 
            Lai = max(mean(MOD15A2H_061_Lai_500m, na.rm=T), 0.0), n_Lai = sum(!is.na(MOD15A2H_061_Lai_500m)), .groups = 'drop')

# interpolate months with missed values
file2_interp <- file2_interp %>% group_by(ID) %>% mutate(Fpar_interp = na.approx(Fpar, rule = 2), Lai_interp = na.approx(Lai, rule = 2))

#####
file3  <- read.csv(file.path(dir_rawdata, 'towers-MYD17A2HGF-061-results.csv'))  ##GPP

# remove low quality data
file3$MYD17A2HGF_061_Gpp_500m[file3$MYD17A2HGF_061_Psn_QC_500m_MODLAND_Description == "Other quality (back-up algorithm or fill values)"] <- NA
file3$Date <- as.Date(file3$Date, format = "%m/%d/%y")

# Alternative way to calculate annual average values
file3_interp <- file3 %>% mutate(month = month(Date)) %>% group_by(ID, month) %>% 
  summarise(GPP = max(mean(MYD17A2HGF_061_Gpp_500m, na.rm=T), 0.0), n_GPP = sum(!is.na(MYD17A2HGF_061_Gpp_500m)), .groups = 'drop')

file3_interp <- file3_interp %>% group_by(ID) %>% mutate(GPP_interp = na.approx(GPP, rule = 2))

# site average data: vegetation index
VI <- file1_interp %>% group_by(ID) %>% summarise(EVI=mean(EVI_interp, na.rm=T), NDVI=mean(NDVI_interp, na.rm=T)) 
# LAI: unitless
LAI <- file2_interp %>% group_by(ID) %>% summarise(Fpar=mean(Fpar_interp, na.rm=T), LAI=mean(Lai_interp, na.rm=T)) 
# GPP: kgC/m2/8d -> convert to kgC/m2/year
GPP <- file3_interp %>% group_by(ID) %>% summarise(GPP=mean(GPP_interp, na.rm=T)/8*365.25) 
data.spectral <- VI %>% left_join(LAI, by='ID') %>% left_join(GPP, by="ID")
# site year average data

#### the correlations between GPP, LAI, EVI, and NDVI of these sites; need to take a look at this first!  
data.spectral.tower <- data.spectral %>% left_join(stat.climate, by=c("ID"="site_ID" ))
corrplot(cor(data.spectral.tower[, 2:9]),
         method = "number",
         sig.level = 0.05,
         type = "upper" # show only upper side
)
# LAI, GPP are strongly related to NEE_day; try to use this. 

#--------------------------------combine soil, climate, spectral, and thermal response strength data together-------------
data.TAS_tot <- read.csv(file.path('data', 'outcome_temp.csv'))
data.TAS_tot <- data.TAS_tot %>% rename("TAS_tot" = "TAS", "TAS_totp" = "TASp")
data.TAS <- read.csv(file.path('data', 'outcome_temp_water_gpp.csv'))

acclimation <- site_info[, 1:7] %>% left_join(stat.climate[, 1:8], by = "site_ID") %>% 
  left_join(data.spectral[, c("ID", "EVI", "NDVI", "LAI", "GPP")], by=c("site_ID" = "ID")) %>% 
  left_join(stat.soil[, c('site_ID', 'SOC')], by = "site_ID") %>% 
  left_join(data.TAS_tot[, c('site_ID', 'TAS_tot', 'TAS_totp')], by = "site_ID") %>% 
  left_join(data.TAS[, c('site_ID', 'TAS', 'TASp')], by = "site_ID")

#--------------------------------
# output the data used for identify drivers of thermal acclimation
write.csv(acclimation, file=file.path('data', 'acclimation_data.csv'), row.names = F)
