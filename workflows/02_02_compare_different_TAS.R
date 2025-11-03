# This scripts compare different thermal response strength (total, direct, indirect via water, and indirect via GPP)
# We also compare our estimated total thermal response strength with these from field experiments
# Authors: Junna Wang, October, 2025

library(librarian)
shelf(dplyr, ggplot2, tidyr)
rm(list=ls())

outcome_temp <- read.csv('data/outcome_temp.csv')

outcome_temp_water <- read.csv('data/outcome_temp_water.csv')

outcome_temp_gpp <- read.csv('data/outcome_temp_gpp.csv')

outcome_temp_water_gpp <- read.csv('data/outcome_final.csv')

outcome <- data.frame(site_ID = outcome_temp$site_ID, TAS_tot = outcome_temp$TAS, TAS_totp = outcome_temp$TASp)

outcome$TAS <- outcome_temp_water_gpp$TAS
outcome$TASp <- outcome_temp_water_gpp$TASp

outcome$TAS_water <- outcome_temp_gpp$TAS - outcome$TAS
outcome$TAS_gpp <- outcome_temp_water$TAS - outcome$TAS
outcome$TAS_water_gpp <- outcome$TAS_tot - outcome$TAS - outcome$TAS_water - outcome$TAS_gpp

outcome$TAS_app <- outcome$TAS_tot - outcome$TAS
# outcome$TAS_water <- outcome$TAS_tot - outcome_temp_water$TAS
# outcome$TAS_gpp <- outcome$TAS_tot - outcome_temp_gpp$TAS


# comparison of different thermal response strengths
outcome %>% tidyr::pivot_longer(col = c("TAS_tot", "TAS", "TAS_app"), names_to = 'TAS_type', values_to = "TAS_value") %>%
  ggplot(aes(x = TAS_type, y = TAS_value)) + 
  geom_boxplot() +
  theme_bw()
  

outcome %>% tidyr::pivot_longer(col = c("TAS_tot", "TAS", "TAS_water", "TAS_gpp", "TAS_water_gpp"), names_to = 'TAS_type', values_to = "TAS_value") %>%
  mutate(TAS_type = factor(TAS_type, levels = c("TAS_tot", "TAS", "TAS_gpp", "TAS_water", "TAS_water_gpp"))) %>%
  ggplot(aes(x = TAS_type, y = TAS_value)) + 
  geom_boxplot() +
  theme_bw()

site_info <- read.csv('data/site_info.csv')
if (!"IGBP" %in% names(outcome)) {
  outcome <- outcome %>% left_join(site_info[, c('site_ID', 'IGBP', 'Climate_class')], by = 'site_ID')
}

outcome %>% tidyr::pivot_longer(col = c("TAS_tot", "TAS", "TAS_app"), names_to = 'TAS_type', values_to = "TAS_value") %>% 
  ggplot(aes(x = Climate_class, y = TAS_value, col = TAS_type)) + 
  geom_boxplot() +
  theme_bw()

outcome %>% tidyr::pivot_longer(col = c("TAS_tot", "TAS", "TAS_app"), names_to = 'TAS_type', values_to = "TAS_value") %>% 
  mutate(IGBP = factor(IGBP, levels = c('OSH', 'SAV', 'WSA', 'CSH', 'DBF', 'DNF', 'EBF', 'MF', 'ENF', 'GRA', 'WET'))) %>%
  ggplot(aes(x = IGBP, y = TAS_value, col = TAS_type)) + 
  geom_boxplot() +
  theme_bw()

#-------------------------------------
# compare with field experiment data
r3_field <- read.csv('data/r3_field.csv')

TAS_estimate_field <- data.frame(r3 = outcome$TAS_tot, respiration = 'ecosystem', result = ifelse(outcome$TAS_totp < 0.05, 'YES', 'NO'))
TAS_estimate_field <- rbind(TAS_estimate_field, r3_field[, c('r3', 'respiration', 'result')])

TAS_estimate_field %>% filter(result == 'YES' & r3 < 0) %>%
  ggplot(aes(x = respiration, y = r3)) +
  geom_boxplot() + 
  theme_bw()

boxplot(r3_field$r3)


t.test(r3_field$r3[r3_field$result == 'YES'], outcome$TAS_tot[outcome$TAS_tot < 0.0 & outcome$TAS_totp < 0.01])
