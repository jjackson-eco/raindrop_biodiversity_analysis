####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##           Incorporating weather data           ##
##                                                ##
##                 Mar 11th 2022                  ##
##                                                ##
####################################################

rm(list = ls())
options(width = 100)

library(tidyverse)
library(psych)
library(MASS)
library(brms)
library(patchwork)
library(viridis)

# colours
raindrop_colours <-
  tibble(treatment = c("Ambient", "Control", "Drought", "Irrigated"),
         num = rep(100,4),
         colour = c("#417A0B", "#BFBFBF", "#E69F00", "#6ECCFF"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load and wrangle data ####

#_________________________________________________________________
## 1a. Biodiversity data
load("../../RainDropRobotics/Data/raindrop_biodiversity_2016_2021.RData",
     verbose = TRUE)

# Converting to total biomass - Here we are considering TRUE 0 observations - need to chat to Andy about this
totbiomass <- biomass %>%
  group_by(year, harvest, block, treatment) %>%
  summarise(tot_biomass = sum(biomass_g*4),
            log_tot_biomass = log(tot_biomass + 1)) %>%
  ungroup() %>%
  mutate(year_f = as.factor(year),
         year_s = as.numeric(scale(year)))

## group biomass - Focussing on only graminoids, forbs and legumes
biomass_gr <- biomass %>%
  filter(group %in% c("Forbs", "Graminoids", "Legumes")) %>%
  mutate(biomass_log = log((biomass_g*4) + 1),
         year_f = as.factor(year),
         year_s = as.numeric(scale(year)))

# Percent cover - total biodiversity indices
percent_cover <- percent_cover %>%
  drop_na(percent_cover) %>%
  filter(species_level == "Yes")

diversity <- percent_cover %>%
  # first convert percentages to a proportion
  mutate(proportion = percent_cover/100) %>%
  group_by(year, month, block, treatment) %>%
  # diversity indices for each group
  summarise(richness = n(),
            simpsons = sum(proportion^2),
            shannon = -sum(proportion * log(proportion))) %>%
  ungroup() %>%
  mutate(year_f = as.factor(year),
         year_s = as.numeric(scale(year)),
         shannon = as.numeric(scale(shannon)),
         simpsons = as.numeric(scale(simpsons)))

#_________________________________________________________________
## 1b. Weather data

# Julian day info - checked
date_info <- tibble(date = seq.Date(as.Date("2016-01-01"),
                                    as.Date("2020-12-31"), by = 1)) %>%
  mutate(day = as.numeric(format.Date(date, format = "%j")),
         year = as.numeric(format.Date(date, format = "%Y")),
         year_day = paste0(day, "_", year))

# load weather data
ECN <- read_csv("../../RainDropRobotics/Data/ECN_MA_Automated_Weather_Station_data_2016_2020_compatible_with_ECN_Database_headers.csv") %>%
  dplyr::select(year = 1, day = 2, hour = 3, temp = 7, rainfall = 10,
                humidity = 17, wind = 8) %>%
  mutate(year_day = paste0(day, "_", year)) %>%
  filter(temp >= -274 & rainfall >= 0 &
           humidity >= 0 & wind >= 0 &
           is.na(year) == F & is.na(day) == F) %>%
  left_join(x = ., y = date_info[,c(1,4)], by = "year_day") %>%
  dplyr::select(year, day, date, hour, temp, rainfall, humidity, wind)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Exploring raw weather data ####

# jpeg(filename = "output/weather_covariance.jpeg",
#      width = 15, height = 15, units = "cm", res = 600)
# pairs.panels(ECN[,5:8], ellipses = FALSE, smooth = FALSE)
# dev.off()

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Generating weather summary variables ####

weather_summary <- ECN %>%
  # Adding in seasons
  mutate(season = case_when(
    date < "2016-03-21" ~ "winter_2016",
    date >= "2016-03-21" & date < "2016-06-21" ~ "spring_2016",
    date >= "2016-06-21" & date < "2016-09-23" ~ "summer_2016",
    date >= "2016-09-23" & date < "2016-12-21" ~ "autumn_2016",
    date >= "2016-12-21" & date < "2017-03-21" ~ "winter_2017",
    date >= "2017-03-21" & date < "2017-06-21" ~ "spring_2017",
    date >= "2017-06-21" & date < "2017-09-23" ~ "summer_2017",
    date >= "2017-09-23" & date < "2017-12-21" ~ "autumn_2017",
    date >= "2017-12-21" & date < "2018-03-21" ~ "winter_2018",
    date >= "2018-03-21" & date < "2018-06-21" ~ "spring_2018",
    date >= "2018-06-21" & date < "2018-09-23" ~ "summer_2018",
    date >= "2018-09-23" & date < "2018-12-21" ~ "autumn_2018",
    date >= "2018-12-21" & date < "2019-03-21" ~ "winter_2019",
    date >= "2019-03-21" & date < "2019-06-21" ~ "spring_2019",
    date >= "2019-06-21" & date < "2019-09-23" ~ "summer_2019",
    date >= "2019-09-23" & date < "2019-12-21" ~ "autumn_2019",
    date >= "2019-12-21" & date < "2020-03-21" ~ "winter_2020",
    date >= "2020-03-21" & date < "2020-06-21" ~ "spring_2020",
    date >= "2020-06-21" & date < "2020-09-23" ~ "summer_2020",
    date >= "2020-09-23" & date < "2020-12-21" ~ "autumn_2020",
    date >= "2020-12-21" ~ "winter_2021")) %>%
  group_by(season) %>%
  # Weather summary variables
  summarise(year = year[n()],
            # temperature
            max_temp = max(temp),
            min_temp = min(temp),
            mean_temp = mean(temp),
            cv_temp = sd(temp) / mean(temp),
            # precipitation
            total_precip = sum(rainfall),
            mean_precip = mean(rainfall),
            cv_precip = sd(rainfall) / mean(rainfall),
            # humidity
            mean_humidity = mean(humidity),
            # wind
            mean_wind = mean(wind))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Summary plots ####

# pairs panels plot
jpeg(filename = "output/weather_summary_covariance.jpeg",
     width = 20, height = 20, units = "cm",res = 1000)
pairs.panels(weather_summary[,-c(1:2)], ellipses = F, lm = TRUE)
dev.off()

# All temp vars >0.79 rho - just use av temp

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Biodiversity data ####

# load data
load("../../RainDropRobotics/Data/raindrop_biodiversity_2016_2021.RData",
     verbose = TRUE)

# Converting to total biomass - Here we are considering TRUE 0 observations - need to chat to Andy about this
totbiomass <- biomass %>%
  group_by(year, harvest, block, treatment) %>%
  summarise(tot_biomass = sum(biomass_g*4),
            log_tot_biomass = log(tot_biomass + 1)) %>%
  ungroup() %>%
  mutate(year_f = as.factor(year),
         year_s = as.numeric(scale(year)))

## group biomass - Focussing on only graminoids, forbs and legumes
biomass_gr <- biomass %>%
  filter(group %in% c("Forbs", "Graminoids", "Legumes")) %>%
  mutate(biomass_log = log((biomass_g*4) + 1),
         year_f = as.factor(year),
         year_s = as.numeric(scale(year)))

# Percent cover - correcting species names
species_names_corrected <- read_csv(file = "../../RainDropRobotics/Data/species_names_corrected.csv")

percent_cover <- percent_cover %>%
  left_join(x =., y = species_names_corrected,
            by = "species") %>%
  mutate(species = if_else(is.na(species_correct) == F,
                           species_correct, species)) %>%
  filter(!note %in% "not_species") %>%
  dplyr::select(-c(species_correct, note)) %>%
  drop_na(percent_cover) %>%
  filter(species_level == "Yes")

# Percent cover - total biodiversity indices
diversity <- percent_cover %>%
  # first convert percentages to a proportion
  mutate(proportion = percent_cover/100) %>%
  group_by(year, month, block, treatment) %>%
  # diversity indices for each group
  summarise(richness = n(),
            simpsons = sum(proportion^2),
            shannon = -sum(proportion * log(proportion))) %>%
  ungroup() %>%
  mutate(year_f = as.factor(year),
         year_s = as.numeric(scale(year)),
         shannon = as.numeric(scale(shannon)),
         simpsons = as.numeric(scale(simpsons)))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Merge weather with biodiversity data ####

## Core weather variables
# - Summer temperature
# - Spring temperature
# - Winter temperature
# - Summer precipitation
# - Spring precipitation
# - Winter precipitation

temp_av <- weather_summary %>%
  mutate(year = as.numeric(year),
         season = gsub("_20[0-9]{2}", "", season),
         temp_av = mean_temp) %>%
  filter(season != "autumn") %>%
  dplyr::select(year, season, temp_av) %>%
  slice(-16) %>%
  pivot_wider(id_cols = year,
              names_from = season,
              values_from = temp_av,
              names_prefix = "temp_")

precip_av <- weather_summary %>%
  mutate(year = as.numeric(year),
         season = gsub("_20[0-9]{2}", "", season),
         precip_av = mean_precip) %>%
  filter(season != "autumn") %>%
  dplyr::select(year, season, precip_av) %>%
  slice(-16) %>%
  pivot_wider(id_cols = year,
              names_from = season,
              values_from = precip_av,
              names_prefix = "precip_")

## Merge with biodiversity
div_simple <- diversity %>%
  filter(month == "June") %>%
  dplyr::select(year, block, treatment, richness, simpsons, shannon)

biodiversity_weather <- totbiomass %>%
  group_by(year, block, treatment) %>%
  summarise(tot_biomass = sum(tot_biomass),
            log_tot_biomass = log(tot_biomass + 1)) %>%
  ungroup() %>%
  left_join(x = ., y = div_simple, by = c("year", "block", "treatment")) %>%
  left_join(x = ., y = temp_av, by = "year") %>%
  left_join(x = ., y = precip_av, by = "year") %>%
  filter(year != 2021)


##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 7. Exploratory plots ####

#_________________________________________________________________
## 7a. Biomass

# Temp
biom_sumt <- ggplot(biodiversity_weather, aes(x = temp_summer,
                                 y = log_tot_biomass,
                                 colour = treatment,
                                 group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour) +
  labs(x = "Average summer temperature (\u00B0C)",
       y = "Above ground net primary productivity (ANPP)", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

biom_sprt <- ggplot(biodiversity_weather, aes(x = temp_spring,
                                              y = log_tot_biomass,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average spring temperature (\u00B0C)",
       y = "Above ground net primary productivity (ANPP)", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

biom_wint <- ggplot(biodiversity_weather, aes(x = temp_winter,
                                              y = log_tot_biomass,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average winter temperature (\u00B0C)",
       y = "Above ground net primary productivity (ANPP)", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

# Precip
biom_sump <- ggplot(biodiversity_weather, aes(x = precip_summer,
                                              y = log_tot_biomass,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour) +
  labs(x = "Average summer precipitation (mm)",
       y = "Above ground net primary productivity (ANPP)", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

biom_sprp <- ggplot(biodiversity_weather, aes(x = precip_spring,
                                              y = log_tot_biomass,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average spring precipitation (mm)",
       y = "Above ground net primary productivity (ANPP)", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

biom_winp <- ggplot(biodiversity_weather, aes(x = precip_winter,
                                              y = log_tot_biomass,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average winter precipitation (mm)",
       y = "Above ground net primary productivity (ANPP)", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

ggsave((biom_wint + biom_sprt + biom_sumt)/
         (biom_winp + biom_sprp + biom_sump),
       filename = "output/biomass_weather.jpeg",
       width = 39, height = 23, units = "cm", dpi = 1000)

#_________________________________________________________________
## 7b. Simpsons

# Temp
simp_sumt <- ggplot(biodiversity_weather, aes(x = temp_summer,
                                              y = simpsons,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour) +
  labs(x = "Average summer temperature (\u00B0C)",
       y = "Simpson's biodiversity index", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

simp_sprt <- ggplot(biodiversity_weather, aes(x = temp_spring,
                                              y = simpsons,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average spring temperature (\u00B0C)",
       y = "Simpson's biodiversity index", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

simp_wint <- ggplot(biodiversity_weather, aes(x = temp_winter,
                                              y = simpsons,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average winter temperature (\u00B0C)",
       y = "Simpson's biodiversity index", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

# Precip
simp_sump <- ggplot(biodiversity_weather, aes(x = precip_summer,
                                              y = simpsons,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour) +
  labs(x = "Average summer precipitation (mm)",
       y = "Simpson's biodiversity index", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

simp_sprp <- ggplot(biodiversity_weather, aes(x = precip_spring,
                                              y = simpsons,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average spring precipitation (mm)",
       y = "Simpson's biodiversity index", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

simp_winp <- ggplot(biodiversity_weather, aes(x = precip_winter,
                                              y = simpsons,
                                              colour = treatment,
                                              group = interaction(block, treatment))) +
  geom_line(show.legend = F) + geom_point(alpha = 0.8, size = 3) +
  scale_colour_manual(values = raindrop_colours$colour, guide = "none") +
  labs(x = "Average winter precipitation (mm)",
       y = "Simpson's biodiversity index", colour = "Treatment") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

ggsave((simp_wint + simp_sprt + simp_sumt)/
         (simp_winp + simp_sprp + simp_sump),
       filename = "output/simpsons_weather.jpeg",
       width = 39, height = 23, units = "cm", dpi = 1000)

#_________________________________________________________________
## 7c. Linear regressions (only 5 years so reduced df)

## Biomass
biom_sprt_mod <- lm(log_tot_biomass ~ temp_spring,
                    data = biodiversity_weather)
summary(biom_sprt_mod)

biom_sumt_mod <- lm(log_tot_biomass ~ temp_summer,
                    data = biodiversity_weather)
summary(biom_sumt_mod)

biom_sprp_mod <- lm(log_tot_biomass ~ precip_spring,
                    data = biodiversity_weather)
summary(biom_sprp_mod)

biom_sump_mod <- lm(log_tot_biomass ~ precip_summer,
                    data = biodiversity_weather)
summary(biom_sump_mod)


## Simpson
simp_sprt_mod <- lm(simpsons ~ temp_spring,
                    data = biodiversity_weather)
summary(simp_sprt_mod)

simp_sumt_mod <- lm(simpsons ~ temp_summer,
                    data = biodiversity_weather)
summary(simp_sumt_mod)

simp_sprp_mod <- lm(simpsons ~ precip_spring,
                    data = biodiversity_weather)
summary(simp_sprp_mod)




