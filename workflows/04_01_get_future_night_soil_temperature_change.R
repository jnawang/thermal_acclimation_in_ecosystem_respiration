# This script estimates monthly nighttime air temperature change from the current (2000-2020) to the future (2041-2060). 
# Authors: Junna Wang, October, 2025

# This script needs to read data stored on data server. 
# In order to reproduce this script, you have to download data from https://www.worldclim.org/
# The raw data includes wc2.1_2.5m_tmin in 2000-2020 and wc2.1_2.5m_tmin_2041-2060 of 13 global circulation models. 

# it takes ~30min to run through this script. 

library(librarian)
shelf(dplyr, terra)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
# dir_rawdata <- '/Users/junnawang/YaleLab/data_server/'
####################End Attention

site_info <- read.csv('data/site_info.csv')
xy <- data.frame(x=site_info$LONG, y=site_info$LAT)  

#----------------Step 1: get monthly temperature of the current period
# dir <- '/Users/jw2946/Documents/data/climate/wc2.1_2.5m_tmin/'
dir <- file.path(dir_rawdata, 'Climate', 'wc2.1_2.5m_tmin')
#
files <- list.files(dir)
base_line <- terra::rast()
for (i in 1:12) {
  if (i < 10) {
    month <- paste0(0, i)
  } else {
    month <- i
  }
  print(month)
  id <- grepl(pattern=paste0(month, '.tif$'), files)
  base_line <- c(base_line, terra::app(terra::rast(file.path(dir, files[id])), fun=mean))
}

names(base_line) <- 1:12
plot(base_line)

#----------------Step 2: get monthly temperature of the future period under emissions scenario (SSP245)
files <- list.files(path=file.path(dir_rawdata, 'Climate', 'wc2.1_2.5m_tmin_2041-2060'), pattern='.tif$', full.names = T)  # 12 future models
nfiles <- length(files)
#
ssps <- c('ssp245')            #  c('ssp126', 'ssp245', 'ssp370', 'ssp585')
for (ssp in ssps) {
  id.select <- grepl(pattern=ssp, files)
  #
  raster_list <- lapply(files[id.select], rast)
  GCM <- list()
  #
  for (i in 1:nlyr(raster_list[[1]])) {
    #  i = 1
    print(i)
    # Extract the layer from each raster
    layers <- lapply(raster_list, function(x) x[[i]])
    
    # Compute the average of layers
    GCM[[i]] <- terra::app(terra::rast(layers), fun=mean)
  }
  #
  GCM <- terra::rast(GCM)
  #
  Tmin_change2010_2050 <- GCM - base_line
  # # take a look at change in summer
  # mean(terra::values(Tmin_change2010_2050[[6]]), na.rm=T)
  #
  #----------------extract temperature change at our study sites----------------
  icell <- adjacent(Tmin_change2010_2050, cellFromXY(Tmin_change2010_2050, xy), include=TRUE)
  tmp <- terra::extract(Tmin_change2010_2050, c(icell))
  Tmin_month <- data.frame(site_ID=site_info$site_ID)
  for (i in 1:ncol(tmp)) {
    Tmin_month[,i+1] <- as.vector(tapply(tmp[,i], rep(1:nrow(xy), times=5), mean, na.rm=T))
  }
  colnames(Tmin_month)[2:13] <- paste0('Tmin', 1:12)
  write.csv(Tmin_month, file=paste0('data/Tmin_month_', ssp, '_wc.csv'), row.names=FALSE)
  #
}

