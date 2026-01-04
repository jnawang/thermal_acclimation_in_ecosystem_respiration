# This scripts compare different thermal response strength (total, direct, indirect via water and GPP)
# We also compare our estimated total thermal response strength with these from field experiments
# Authors: Junna Wang, October, 2025

# it takes 1 min to run. 

library(librarian)
shelf(dplyr, ggplot2, tidyr)
rm(list=ls())

outcome_temp <- read.csv(file.path('data', 'outcome_temp.csv'))

outcome_temp_water_gpp <- read.csv(file.path('data', 'outcome_temp_water_gpp.csv'))

outcome <- data.frame(site_ID = outcome_temp$site_ID, TAS_tot = outcome_temp$TAS, TAS_totp = outcome_temp$TASp)

outcome$TAS <- outcome_temp_water_gpp$TAS
outcome$TASp <- outcome_temp_water_gpp$TASp

# outcome$TAS_water <- outcome_temp_gpp$TAS - outcome$TAS
# outcome$TAS_gpp <- outcome_temp_water$TAS - outcome$TAS
# outcome$TAS_water_gpp <- outcome$TAS_tot - outcome$TAS - outcome$TAS_water - outcome$TAS_gpp

outcome$TAS_app <- outcome$TAS_tot - outcome$TAS

# comparison of different thermal response strengths
outcome %>% tidyr::pivot_longer(col = c("TAS_tot", "TAS", "TAS_app"), names_to = 'TAS_type', values_to = "TAS_value") %>%
  ggplot(aes(x = TAS_type, y = TAS_value)) + 
  geom_boxplot() +
  theme_bw()
  

site_info <- read.csv(file.path('data', 'site_info.csv'))
if (!"IGBP" %in% names(outcome)) {
  outcome <- outcome %>% left_join(site_info[, c('site_ID', 'IGBP', 'Climate_class')], by = 'site_ID')
}

IGBPs <- unique(outcome$IGBP)
Climate_classes <- unique(outcome$Climate_class)

df.TAS.IGBP <- data.frame(IGBP = character(), n = integer(), total0 = double(), direct0 = double(), app0 = double())
for (i in 1:length(IGBPs)) {
  df.TAS.IGBP[i, 1] <- IGBPs[i]
  df.TAS.IGBP[i, 2] <- sum(outcome$IGBP == IGBPs[i])
  if (df.TAS.IGBP[i, 2] >= 3) {
    df.TAS.IGBP[i, 3] <- t.test(outcome$TAS_tot[outcome$IGBP == IGBPs[i]], mu=0)$p.value
    df.TAS.IGBP[i, 4] <- t.test(outcome$TAS[outcome$IGBP == IGBPs[i]], mu=0)$p.value
    df.TAS.IGBP[i, 5] <- t.test(outcome$TAS_app[outcome$IGBP == IGBPs[i]], mu=0)$p.value
  }
}

t.test(outcome$TAS_tot[outcome$IGBP %in% c("OSH", "SAV", "WSA")], mu=0) # 0.01805551, p-value = 0.015
t.test(outcome$TAS_tot[outcome$IGBP %in% c("GRA", "ENF")], mu=0) # 

t.test(outcome$TAS_tot[outcome$IGBP %in% c("EBF", "DBF", "MF", "DNF", "CSH", "WET")], mu=0)

# broadleaf forests had trends for enhancing TAS. 
t.test(outcome$TAS[outcome$IGBP %in% c("DBF", "EBF", "MF", "DNF", "CSH", "WET")], mu=0) # 0.01805551, p-value = 0.015

t.test(outcome$TAS[outcome$IGBP %in% c("WET")], mu=0) 
t.test(outcome$TAS_app[outcome$IGBP %in% c("WET")], mu=0)  # -0.006572169, t = -1.7525, df = 13, p-value = 0.1032



df.TAS.climate <- data.frame(climate = character(), n = integer(), total0 = double(), direct0 = double(), app0 = double())
for (i in 1:length(Climate_classes)) {
  df.TAS.climate[i, 1] <- Climate_classes[i]
  df.TAS.climate[i, 2] <- sum(outcome$Climate_class == Climate_classes[i])
  if (df.TAS.climate[i, 2] >= 3) {
    df.TAS.climate[i, 3] <- t.test(outcome$TAS_tot[outcome$Climate_class == Climate_classes[i]], mu=0)$p.value
    df.TAS.climate[i, 4] <- t.test(outcome$TAS[outcome$Climate_class == Climate_classes[i]], mu=0)$p.value
    df.TAS.climate[i, 5] <- t.test(outcome$TAS_app[outcome$Climate_class == Climate_classes[i]], mu=0)$p.value
  }
}

t.test(outcome$TAS[outcome$Climate_class %in% c("Csa")], mu=0)  # 0.02203213, p-value = 0.09531
t.test(outcome$TAS_tot[outcome$Climate_class %in% c("Bsh", "Bsk", "Bwk")], mu=0)  # -0.08888042, p-value = 0.0001909


# look at total response
t.test(outcome$TAS_tot[outcome$Climate_class %in% c("Cfa", "Cfb", "Cfc")], mu=0)  # p-value = 0.09921
t.test(outcome$TAS_tot[outcome$Climate_class %in% c("Dfa", "Dfb")], mu=0)  # p-value = 0.06449
t.test(outcome$TAS_tot[outcome$Climate_class %in% c("Dfc", "Dfd", "Dwc")], mu=0)  # p-value = 0.2231
t.test(outcome$TAS_tot[outcome$Climate_class %in% c("Cfa", "Cfb", "Cfc", "Dfa", "Dfb", "Dfc", "Dfd", "Dwc", "Csa", "Csb", "ET", "Af", "Am")], mu=0) 
t.test(outcome$TAS_tot[outcome$Climate_class %in% c("Csa")], mu=0)  # p-value = 0.09921


var(outcome$TAS_tot)
var(outcome$TAS_app)
var(outcome$TAS)
cov(outcome$TAS_app, outcome$TAS)
(var(outcome$TAS_app) + cov(outcome$TAS_app, outcome$TAS)) / var(outcome$TAS_tot)
(var(outcome$TAS) + cov(outcome$TAS_app, outcome$TAS)) / var(outcome$TAS_tot)


summary(lm(data = outcome, TAS_tot ~ IGBP))  # OSH, SAV, WSA are significantly different. 
summary(lm(data = outcome, TAS_tot ~ Climate_class))  # Bsh, Bsk, Bwk, Csb, and ET are significantly different. 

summary(lm(data = outcome, TAS ~ IGBP))  # OSH, SAV, WSA are significantly different. p = 0.06791
summary(lm(data = outcome, TAS ~ Climate_class))  # Bsh, Bsk, Bwk, Csb, and ET are significantly different. p = 0.04067

summary(lm(data = outcome, TAS_app ~ IGBP))           #  no sig;   p-value: 0.07353
summary(lm(data = outcome, TAS_app ~ Climate_class))  #  yes sig;  p-value: 9.513e-09

#-------------------------------------------------------------------------------------------------------------------------------

outcome %>% tidyr::pivot_longer(col = c("TAS_tot", "TAS", "TAS_app"), names_to = 'TAS_type', values_to = "TAS_value") %>% 
  ggplot(aes(x = Climate_class, y = TAS_value, col = TAS_type)) + 
  geom_boxplot() +
  theme_bw()

outcome %>% tidyr::pivot_longer(col = c("TAS_tot", "TAS", "TAS_app"), names_to = 'TAS_type', values_to = "TAS_value") %>% 
  mutate(IGBP = factor(IGBP, levels = c('OSH', 'SAV', 'WSA', 'CSH', 'DBF', 'DNF', 'EBF', 'MF', 'ENF', 'GRA', 'WET'))) %>%
  ggplot(aes(x = IGBP, y = TAS_value, col = TAS_type)) + 
  geom_boxplot() +
  theme_bw()
