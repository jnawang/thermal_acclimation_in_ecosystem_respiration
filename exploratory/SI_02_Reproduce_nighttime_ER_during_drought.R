# Test model performance for drought events to address reviewer's concern. 
# Questions to address
# 1) what's the best way of adding water?


library(librarian)
shelf(dplyr, lubridate, gslnls, caret, performance, ggpubr, ggplot2, zoo, brms, cowplot)
rm(list=ls())

####################Attention: change this directory based on your own directory of raw data
# dir_rawdata <- '/Volumes/MaloneLab/Research/Stability_Project/Thermal_Acclimation'
dir_rawdata <- '/Users/junnawang/YaleLab/data_server/'
####################End Attention

name_site <- "DK-Sor"    # 'CH-LAE'
path <- file.path(dir_rawdata, "RespirationData", paste0(name_site, '_nightNEE.csv'))
if (file.exists(path)) {
  a_measure_night_complete <- read.csv(path)
  ac <- read.csv(file.path(dir_rawdata, "RespirationData", paste0(name_site, '_ac.csv')))
} else {
  next
}

a_measure_night_complete <- a_measure_night_complete %>% mutate(DATETIME=ISOdatetime(YEAR, MONTH, DAY, HOUR, MINUTE, 0, tz = "UTC"))
ac <- ac %>% mutate(DATETIME=ISOdatetime(YEAR, MONTH, DAY, HOUR, MINUTE, 0, tz = "UTC"))

# get daily GPP
ac_day <- ac %>% filter(daytime) %>% group_by(YEAR, DOY) %>% 
  summarise(GPP_DT1 = sum(GPP_DT, na.rm=T) * 0.5 / 24, .groups = "drop")      # average GPP across a year

a_measure_night_complete <- a_measure_night_complete %>% mutate(DOY_gpp = case_when(HOUR >= 12 ~ DOY,
                                                                                    HOUR < 12 ~ DOY - 1)) %>%
  left_join(ac_day, by = c("YEAR", "DOY_gpp"="DOY")) %>% select(-"DOY_gpp")

# find the drought year and NEE characteristics
plot(a_measure_night_complete$DATETIME, a_measure_night_complete$SWC)
plot(a_measure_night_complete$DATETIME[a_measure_night_complete$YEAR %in% 2018:2020], a_measure_night_complete$SWC[a_measure_night_complete$YEAR %in% 2018:2020])
plot(a_measure_night_complete$DATETIME[a_measure_night_complete$YEAR %in% 2018:2020], a_measure_night_complete$NEE[a_measure_night_complete$YEAR %in% 2018:2020])

# 2015:2020
p1 <- a_measure_night_complete %>% filter(YEAR %in% 2019 & between(DOY, 121, 240)) %>% mutate(YEAR=as.factor(YEAR)) %>%
  ggplot(aes(x=DOY, y=SWC)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 170, linetype = "dotted") +
  geom_vline(xintercept = 205, linetype = "dotted") +
  labs(x = 'Day of year', y="Soil water content (%)", tag = 'a') +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))


p2 <- a_measure_night_complete %>% filter(YEAR %in% 2019 & between(DOY, 121, 240)) %>% mutate(YEAR=as.factor(YEAR)) %>% 
  ggplot(aes(x=DOY, y=NEE)) + 
  geom_point(alpha = 0.5) + 
  geom_vline(xintercept = 170, linetype = "dotted") +
  geom_vline(xintercept = 205, linetype = "dotted") +  
  labs(x = 'Day of year', y="Night NEE (umolCO2 m-2 s-1)", tag = 'b') + 
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))


p3 <- a_measure_night_complete %>% filter(YEAR %in% 2019 & between(DOY, 170, 205)) %>% mutate(YEAR=as.factor(YEAR)) %>% 
  ggplot(aes(x=TS, y=NEE)) + 
  geom_point(alpha = 0.5) +
  labs(x = "Soil temperature (°C)", y="Night NEE (umolCO2 m-2 s-1)", tag = 'c') + 
  theme_bw()

p4 <- a_measure_night_complete %>% filter(YEAR %in% 2019 & between(DOY, 170, 205)) %>% mutate(YEAR=as.factor(YEAR)) %>% 
  ggplot(aes(x=SWC, y=NEE)) + 
  geom_point(alpha = 0.5) + 
  labs(x = "Soil water content (%)", y="Night NEE (umolCO2 m-2 s-1)", tag = 'd') + 
  theme_bw()

plot_grid(p1, p2, p3, p4, nrow=2, ncol=2)
ggsave('exploratory/drought_example.png')

# 2008 and 2020 both suffered from drought; VPD can also represent dryness. 

# work on data in July and August
frmu_temp <- NEE ~ exp(alpha * TS + beta*TS^2) * C0
frmu_temp_water <- NEE ~ exp(alpha * TS + beta*TS^2) * C0 * SWC / (Hs + SWC)
frmu_temp_water_gpp <- NEE ~ exp(alpha * TS + beta*TS^2) * C0 * SWC / (Hs + SWC) + k2 * GPP_DT1 * exp(gamma * TS)

priors_temp <- brms::prior("normal(0.1, 0.5)", nlpar = "alpha", lb = 0, ub = 0.2) +
  brms::prior("normal(-0.001, 0.1)", nlpar = "beta", lb = -0.01, ub = 0.01) +
  brms::prior("normal(2, 5)", nlpar = "C0", lb = 0, ub = 20) 

priors_temp_water <- priors_temp +
  brms::prior("normal(10, 10)", nlpar = "Hs", lb = 0, ub = 1000) 

priors_temp_water_gpp <- priors_temp_water +
  brms::prior("normal(0.2, 1)", nlpar = "k2", lb = 0, ub = 1) + 
  brms::prior("normal(0.05, 0.5)", nlpar = "gamma", lb = 0, ub = 0.2)

# drought occurred in 2020: 2019:2020; 210, 290
data <- a_measure_night_complete %>% filter(YEAR %in% 2019 & between(DOY, 170, 205)) %>% filter(!is.na(SWC))

model_temp <- brms::brm(brms::bf(frmu_temp, alpha+beta+C0 ~ 1, nl = TRUE),
                        prior = priors_temp, data = data, iter = 2000, cores = 4, chains = 4, backend = "cmdstanr", 
                        control = list(adapt_delta = 0.95, max_treedepth = 15))
print(summary(model_temp), digits = 4)
data$NEE_temp <- fitted(model_temp)[, "Estimate"]
postResample(pred=data$NEE_temp, obs=data$NEE)

data$residual_NEE_temp <- data$NEE - data$NEE_temp
cor.test(data$residual_NEE_temp, data$SWC)


model_temp_water <- brms::brm(brms::bf(frmu_temp_water, alpha+beta+C0+Hs ~ 1, nl = TRUE),
                        prior = priors_temp_water, data = data, iter = 2000, cores = 4, chains = 4, backend = "cmdstanr", 
                        control = list(adapt_delta = 0.95, max_treedepth = 15))
print(summary(model_temp_water), digits = 4)
data$NEE_temp_water <- fitted(model_temp_water)[, "Estimate"]
postResample(pred=data$NEE_temp_water, obs=data$NEE)

data$residual_NEE_temp_water <- data$NEE - data$NEE_temp_water
cor.test(data$residual_NEE_temp_water, data$SWC)


model_temp_water_gpp <- brms::brm(brms::bf(frmu_temp_water_gpp, alpha+beta+C0+Hs+k2+gamma ~ 1, nl = TRUE),
                              prior = priors_temp_water_gpp, data = data, iter = 2000, cores = 4, chains = 4, backend = "cmdstanr", 
                              control = list(adapt_delta = 0.95, max_treedepth = 15))
print(summary(model_temp_water_gpp), digits = 4)
data$NEE_temp_water_gpp <- fitted(model_temp_water_gpp)[, "Estimate"]
postResample(pred=data$NEE_temp_water_gpp, obs=data$NEE)


data %>% filter(YEAR==2019) %>% 
  ggplot(aes(x=DOY)) +
  geom_point(aes(y=NEE)) + 
  geom_line(aes(y = NEE_temp, color = "Temperature")) +
  geom_line(aes(y = NEE_temp_water, color = "Temp + Water")) +
  scale_color_manual(
    name = "Model",                     # legend title
    values = c("Temperature" = "red", 
               "Temp + Water" = "blue")) +
  labs(y="Nighttime NEE (umolCO2 m-2 s-1)", tag = 'a', title = "2019") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))

data %>% filter(YEAR==2020) %>% 
  ggplot(aes(x=DOY)) +
  geom_point(aes(y=NEE)) + 
  geom_line(aes(y = NEE_temp, color = "Temperature")) +
  geom_line(aes(y = NEE_temp_water, color = "Temp + Water")) +
  scale_color_manual(
    name = "Model",                     # legend title
    values = c("Temperature" = "red", 
               "Temp + Water" = "blue")) +
  labs(y="Nighttime NEE (umolCO2 m-2 s-1)", tag = 'b', title = "2020") +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8))

# the first example:
# "Dk-Sor", 2019 (wet) and 2020 (dry); using data in July and August (DOY: 210~290). 
# Adding GPP almost has no effect for this site. 


