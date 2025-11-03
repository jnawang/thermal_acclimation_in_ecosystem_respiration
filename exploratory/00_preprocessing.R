rm(list=ls())

library(librarian)
shelf(ggplot2, dplyr, tidyr, amerifluxr, data.table, pander, lubridate)


data_acclimation <- read.csv('data/acclimation_data.csv')
sites_AmeriFlux <- read.csv('data/AmeriFlux-site-search-results-all-sites-202411111059.csv')

sites_AmeriFlux_BASE <- sites_AmeriFlux %>% filter(Site.ID %in% data_acclimation$site_ID & AmeriFlux.FLUXNET.Data == 'No')
sites_AmeriFlux_FLUXNET <- sites_AmeriFlux %>% filter(Site.ID %in% data_acclimation$site_ID & AmeriFlux.FLUXNET.Data == 'Yes')

#
data_acclimation <- data_acclimation %>% mutate(source = case_when(site_ID %in% sites_AmeriFlux_BASE$Site.ID ~ "AmeriFlux_BASE",
                                                                   site_ID %in% sites_AmeriFlux_FLUXNET$Site.ID ~ "AmeriFlux_FLUXNET",
                    !site_ID %in% sites_AmeriFlux_BASE$Site.ID & !site_ID %in% sites_AmeriFlux_FLUXNET$Site.ID ~ "FLUXNET"))


# some European sites
files <- list.files('/Users/jw2946/Documents/stability_project/code_ICOS2020/', pattern='.Rmd$')
site_ID_ICOS2020 <- substr(files, start = 26, stop = 31)

data_acclimation$source[data_acclimation$site_ID %in% site_ID_ICOS2020] <- "FLUXNET_ICOS2020"

write.csv(data_acclimation, 'data/acclimation_data_new.csv')


# Can I have more sites?
ameri_10yr <- read.csv('/Users/jw2946/Downloads/AmeriFlux-site-search-results-202507191557_FC_TA_10yr.csv')
#
ameri_omit <- ameri_10yr %>% filter(!Site.ID %in% data_acclimation$site_ID) %>% 
  filter(!Vegetation.Abbreviation..IGBP. %in% c('CRO'))

# how many FLUXNET2025 sites > 10 year data
files_FLUXNET2025 <- list.files(file.path(dir_rawdata, "SiteData", "FLUXNET07202025"), pattern=".zip$", full.names = F)
df <- data.frame(site_ID=integer(), yStart=integer(), yEnd=integer())
for (i in 1:length(files_FLUXNET2025)) {
  df[i,1] <- substr(files_FLUXNET2025[i], 5, 10)
  df[i,2] <- as.integer(substr(files_FLUXNET2025[i], 35, 38))
  df[i,3] <- as.integer(substr(files_FLUXNET2025[i], 40, 43))
}

df$nyear <- df$yEnd - df$yStart + 1
# df$nyear <- 2024 - df$yStart + 1
df_sub <- df %>% filter(nyear>=10)

# unzip files
files_FLUXNET2025 <- list.files(file.path(dir_rawdata, "SiteData", "FLUXNET07202025"), pattern=".zip$", full.names = T)
for (i in 1:length(files_FLUXNET2025)) {
  print(i)
  unzip(files_FLUXNET2025[i], exdir = file.path(dir_rawdata, "SiteData", "FLUXNET07202025", "unzip"))
}




# unzip files
# files <- list.files(file.path(dir_rawdata, "SiteData", "Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1"), pattern="_FLUXNET_HH_L2.zip$", full.names = T)
# for (file in files) {
#   unzip(file, exdir = file.path(dir_rawdata, "SiteData", "Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1", "unzip"))
# }

df_ICOS <- data.frame(site_ID=character(), yStart=integer(), yEnd=integer())
files <- list.files(file.path(dir_rawdata, "SiteData", "Ecosystem final quality (L2) product in ETC-Archive format - release 2025-1", "unzip"), pattern=".csv$", full.names = T)
for (i in 1:length(files)) {
  print(i)
  df_ICOS[i, 1] <- substr(files[i], 166, 171)
  a <- read.csv(files[i])
  df_ICOS[i, 2] <- year(ymd_hm(a$TIMESTAMP_START[1]))
  df_ICOS[i, 3] <- year(ymd_hm(a$TIMESTAMP_START[nrow(a)]))
}

df_FLUXNET_ICOS <- left_join(df, df_ICOS, by="site_ID")
df_FLUXNET_ICOS$yStart_ICOS_corrected <- pmax(df_FLUXNET_ICOS$yStart.y, df_FLUXNET_ICOS$yEnd.x + 1) 
df_FLUXNET_ICOS$nyear_ICOS <- df_FLUXNET_ICOS$yEnd.y - df_FLUXNET_ICOS$yStart_ICOS_corrected + 1
df_FLUXNET_ICOS$nyear_ICOS[is.na(df_FLUXNET_ICOS$nyear_ICOS)] <- 0
df_FLUXNET_ICOS$nyear_total <- df_FLUXNET_ICOS$nyear + df_FLUXNET_ICOS$nyear_ICOS

df_FLUXNET_ICOS$site_ID[df_FLUXNET_ICOS$nyear_total>=10]
setdiff(df_FLUXNET_ICOS$site_ID[df_FLUXNET_ICOS$nyear_total>=10], site_info$site_ID[57:93])


data_HH <- read.csv('/Users/jw2946/Documents/stability_project/SiteData/Ameriflux_FLUXNET/unzip/AMF_US-xGR_FLUXNET_FULLSET_HH_2017-2021_3-5.csv')
data_HH[data_HH == -9999] <- NA
plot(data_HH$TA_F)

summary(lm(data = data_HH[30000:nrow(data_HH),], TS_F_MDS_1 ~ TA_F))

  
  

