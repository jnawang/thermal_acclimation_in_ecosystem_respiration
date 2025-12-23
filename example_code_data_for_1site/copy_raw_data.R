library(librarian)
shelf(dplyr, lubridate)

####################Attention: change this directory based on your own directory of raw data
dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
####################End Attention

files_FLUXNET2015  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2015'), pattern = ".zip$", full.names = TRUE)
files_FLUXNET2020  <- list.files(path=file.path(dir_rawdata, 'SiteData', 'FLUXNET2020'), pattern = ".zip$", full.names = TRUE)
files_ICOS_after2020 <- list.files(path=file.path(dir_rawdata, 'SiteData', 'ICOS_after2020'), pattern = ".zip$", full.names = TRUE)
#
files_FLUXNET2025 <- list.files(file.path(dir_rawdata, "SiteData", "FLUXNET07202025"), pattern=".zip$", full.names = T)
files_ICOS2025 <- list.files(file.path(dir_rawdata, "SiteData", "Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1"), pattern="FLUXNET_(HH|HR)_L2.zip$", full.names = T)

files_AmeriFlux_BASE <- list.files(file.path(dir_rawdata, "SiteData", "AmeriFlux_BASE"), pattern=".zip$", full.names = T)

site_info <- read.csv(file.path('data', 'site_info.csv'))

for (id in 1:nrow(site_info)) {
  # id = 57
  print(id)
  name_site <- site_info$site_ID[id]
  print(name_site)
  
  data_source <- site_info$source[id]
  if (data_source %in% c("AmeriFlux_BASE", "AmeriFlux_FLUXNET")) {
    # files <- files_AmeriFlux_BASE
    # file  <- files[grepl(name_site, files)]    
    # file.copy(file, "/Users/jw2946/Downloads/SiteData2Zenodo/AmeriFlux_BASE/")
  } else {
    data_source <- unlist(strsplit(data_source, "_"))
    for (isource in 1:length(data_source)) {
      if (data_source[isource] == "FLUXNET2015") {
        # files <- files_FLUXNET2015
        # file  <- files[grepl(name_site, files)]
        # file.copy(file, "/Users/jw2946/Downloads/SiteData2Zenodo/FLUXNET2015/")
      } else if (data_source[isource] == "FLUXNET2020") {
        # files <- files_FLUXNET2020
        # file  <- files[grepl(name_site, files)]
        # file.copy(file, "/Users/jw2946/Downloads/SiteData2Zenodo/FLUXNET2020/")
      } else if (data_source[isource] == "ICOS2020") {
        files <- files_ICOS_after2020
        file  <- files[grepl(name_site, files)]
        file.copy(file, "/Users/jw2946/Downloads/SiteData2Zenodo/ICOS_after2020/")
      } else if (data_source[isource] == "FLUXNET2025") {
      #   files <- files_FLUXNET2025
      #   file  <- files[grepl(name_site, files)]
      #   file.copy(file, "/Users/jw2946/Downloads/SiteData2Zenodo/FLUXNET07202025/")
      } else if (data_source[isource] == "ICOS2025") {
        files <- files_ICOS2025
        file  <- files[grepl(name_site, files)]
        file.copy(file, "/Users/jw2946/Downloads/SiteData2Zenodo/Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1/")
      } 
    }
  }
}