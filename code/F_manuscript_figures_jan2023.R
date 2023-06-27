####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##               Manuscript figures               ##
##                                                ##
##                 Sep 11th 2022                  ##
##                                                ##
####################################################

rm(list = ls())
options(width = 100)

# packages
library(tidyverse)
library(brms)
library(tidybayes)
library(patchwork)
library(gghalves)
library(flextable)
library(vegan)

# colours
raindrop_colours <-
  tibble(treatment = c("Ambient", "Control", "Drought", "Irrigated"),
         num = rep(100,4),
         colour = c("#417A0B", "#BFBFBF", "#E69F00", "#6ECCFF"))


##Figure names shifted down by 1!!

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load and wrangle data ####

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

## Biomass variance
biomass_var <- totbiomass %>%
  mutate(plot = paste0(block, "_", treatment)) %>%
  group_by(plot, year) %>%
  summarise(treatment = treatment[1], block = block[1],
            tot_biomass = sum(tot_biomass)) %>%
  ungroup() %>%
  group_by(plot) %>%
  summarise(treatment = treatment[1], block = block[1],
            cv_biomass = -1*(sd(tot_biomass)/mean(tot_biomass))) %>%
  ungroup()

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
         shannon_sc = as.numeric(scale(shannon)),
         simpsons_sc = as.numeric(scale(simpsons)))

# summary stats
totbiomass %>%
  group_by(year, block, treatment) %>%
  summarise(tot = sum(tot_biomass)) %>%
  ungroup() %>%
  summarise(mn = mean(tot), md = median(tot), sd = sd(tot))

biomass_gr %>%
  mutate(totbio = sum(biomass_g)) %>%
  group_by(group) %>%
  summarise(perc = sum(biomass_g)/totbio[1])

percent_cover %>% summarise(n = n_distinct(species))

diversity %>% summarise(mn = mean(richness), md = median(richness),
                        sd = sd(richness))

totbiomass %>%
  group_by(treatment) %>%
  summarise(mn = mean(tot_biomass), md = median(tot_biomass), sd = sd(tot_biomass))

biomass_gr %>%
  group_by(group, treatment) %>%
  summarise(mn = mean(biomass_g), md = median(biomass_g), sd = sd(biomass_g))

biomass_var %>%
  group_by(treatment) %>%
  summarise(mn = -1*mean(cv_biomass))

diversity %>%
  group_by(year) %>%
  summarise(mn = mean(richness))

# top_species_perc %>% group_by(species, year) %>% summarise(mn = mean(proportion))

# group-level biomass percentages
biomass %>%
  mutate(tot_biomass = sum(biomass_g*4)) %>%
  group_by(group) %>%
  summarise(biomass = sum(biomass_g*4),
            biomass_percent = (biomass/tot_biomass[1])*100)

# scaling to 5 m2
rich_5m <- percent_cover %>%
  filter(month == "June") %>%
  mutate(id = paste0(block, "_", treatment)) %>%
  group_by(block, year) %>%
  summarise(richness = n_distinct(species))

summary(rich_5m$richness)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Figure 1 - Biomass reduction ####

load("data/totbiomass_models.RData", verbose = T)
load("data/biomass_models.RData", verbose = T)
load("data/biomass_var_models.RData", verbose = T)

## 2a. Figure 1a. Total biomass
tbsum <- totbiomass %>% group_by(treatment) %>%
  summarise(mn = mean(log_tot_biomass)) %>%
  mutate(trt = c(0.8,1.8, 2.8,3.8))

set.seed(666)
fig1a <- as.data.frame(brms::posterior_predict(totbiomass_treatment)) %>%
  mutate(sim = 1:8000) %>%
  pivot_longer(-sim) %>%
  bind_cols(., slice(totbiomass, rep(1:220, 8000))) %>%
  filter(harvest == "Mid") %>%
  ggplot(aes(x = treatment, y = value, fill = treatment)) +
  stat_halfeye(width = 0.2, alpha = 0.7) +
  geom_half_point(data = totbiomass,
                  aes(y = log_tot_biomass, x = treatment,
                      colour = treatment),
                  side = "l", ## draw jitter on the left
                  range_scale = 0.4, ## control range of jitter
                  alpha = 0.5, size = 3.4) +
  geom_point(data = tbsum, aes(y = mn, x = trt, colour = NULL),
                  alpha = 1, size = 4, colour = "black", shape = 17) +
  scale_fill_manual(values = raindrop_colours$colour,
                    guide = "none", aesthetics = c("fill", "colour")) +
  scale_x_discrete(labels = c("Ambient\nControl", "Procedural\nControl",
                              "Drought", "Irrigated")) +
  labs(tage = "a)", x = "Precipitation treatment",
       y = "Above-ground net primary productivity (ANPP)") +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
fig1a

## 2b. Figure S2. Group-level biomass - moving to supplementary
bgsum <- biomass_gr %>% group_by(group, treatment) %>%
  summarise(mn = mean(biomass_log)) %>%
  ungroup() %>%
  mutate(trt = rep(c(0.8,1.8, 2.8,3.8), 3))

set.seed(666)
figS2 <- as.data.frame(brms::posterior_predict(biomass_treatment_group_int,
                                               ndraws = 2000)) %>% # reducing to make data more manageable
  mutate(sim = 1:2000) %>%
  pivot_longer(-sim) %>%
  bind_cols(., slice(biomass_gr, rep(1:660, 2000))) %>%
  filter(harvest == "Mid") %>%
  ggplot(aes(x = treatment, y = value, fill = treatment, group = group)) +
  stat_halfeye(width = 0.2, alpha = 0.7) +
  geom_half_point(data = biomass_gr,
                  aes(y = biomass_log, x = treatment,
                      colour = treatment, group = NULL),
                  side = "l", ## draw jitter on the left
                  range_scale = 0.4, ## control range of jitter
                  alpha = 0.5, size = 3.2) +
  geom_point(data = bgsum, aes(y = mn, x = trt, colour = NULL),
             alpha = 1, size = 4, colour = "black", shape = 17) +
  scale_fill_manual(values = raindrop_colours$colour,
                    guide = "none", aesthetics = c("fill", "colour")) +
  scale_x_discrete(labels = c("Ambient\nControl", "Procedural\nControl",
                              "Drought", "Irrigated")) +
  labs(x = "Precipitation treatment",
       y = "Above-ground net primary productivity") +
  facet_wrap(~ group, nrow = 3) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11))
figS2

ggsave(figS2, filename = "output/SI/figureS2.jpeg",
              width = 14, height = 19, units = "cm", dpi = 1500)

## 2c. Figure 1b. Temporal stability
bvsum <- biomass_var %>% group_by(treatment) %>%
  summarise(mn = mean(cv_biomass)) %>%
  mutate(trt = c(0.8,1.8, 2.8,3.8))

set.seed(666)
fig1b <- as.data.frame(brms::posterior_predict(biomass_var_treat)) %>%
  mutate(sim = 1:8000) %>%
  pivot_longer(-sim) %>%
  bind_cols(., slice(biomass_var, rep(1:20, 8000))) %>%
  ggplot(aes(x = treatment, y = value, fill = treatment)) +
  stat_halfeye(width = 0.2, alpha = 0.7) +
  geom_half_point(data = biomass_var,
                  aes(y = cv_biomass, x = treatment,
                      colour = treatment),
                  side = "l", ## draw jitter on the left
                  range_scale = 0.2, ## control range of jitter
                  alpha = 0.5, size = 3.4) +
  geom_point(data = bvsum, aes(y = mn, x = trt, colour = NULL),
             alpha = 1, size = 4, colour = "black", shape = 17) +
  scale_fill_manual(values = raindrop_colours$colour,
                    guide = "none", aesthetics = c("fill", "colour")) +
  scale_x_discrete(labels = c("Ambient\nControl", "Procedural\nControl",
                              "Drought", "Irrigated")) +
  labs(tage = "b)", x = "Precipitation treatment",
       y = "Temporal stability in productivity") +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
fig1b

ggsave(fig1a + fig1b,
       filename = "output/manuscript_figures/figure1.jpeg",
       width = 43, height = 20, units = "cm", dpi = 1500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Figure 2 - Community resistance to Drought in RainDrop ####

load("data/NMDS_results.RData")

## fig 2a. Richness treatment
fig2a <- ggplot(diversity, aes(x = treatment, y = richness)) +
  geom_violin(aes(fill = treatment), colour = "white",
              alpha = 0.7, trim = F, show.legend = F) +
  geom_jitter(aes(colour = treatment),
              width = 0.1, alpha = 0.9, size = 2,
              show.legend = F) +
  stat_summary(geom = "point", fun = "mean", colour = "white",
               shape = 21, size = 4, fill = "white", show.legend = F) +
  scale_y_continuous(breaks = seq(0,50, by = 5)) +
  scale_x_discrete(labels = NULL) +
  scale_fill_manual(values = raindrop_colours$colour,
                    guide = "none", aesthetics = c("fill", "colour")) +
  labs(tag = "a)", x = NULL,
       y = "Species richness") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

## fig 2b. Shannon-Weiner treatment
fig2b <- ggplot(diversity, aes(x = treatment, y = shannon)) +
  geom_violin(aes(fill = treatment), colour = "white",
              alpha = 0.7, trim = F, show.legend = F) +
  geom_jitter(aes(colour = treatment),
              width = 0.1, alpha = 0.9, size = 2,
              show.legend = F) +
  stat_summary(geom = "point", fun = "mean", colour = "white",
               shape = 21, size = 4, fill = "white", show.legend = F) +
  scale_x_discrete(labels = NULL) +
  scale_fill_manual(values = raindrop_colours$colour,
                    guide = "none", aesthetics = c("fill", "colour")) +
  labs(tag = "b)", x = NULL,
       y = "Shannon-Weiner index (H)") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

## fig 2c. Simpson's treatment
fig2c <- ggplot(diversity, aes(x = treatment, y = simpsons)) +
  geom_violin(aes(fill = treatment), colour = "white",
              alpha = 0.7, trim = F, show.legend = F) +
  geom_jitter(aes(colour = treatment),
              width = 0.1, alpha = 0.9, size = 2,
              show.legend = F) +
  stat_summary(geom = "point", fun = "mean", colour = "white",
               shape = 21, size = 4, fill = "white", show.legend = F) +
  scale_x_discrete(labels = c("Ambient\nControl", "Procedural\nControl",
                              "Drought", "Irrigated")) +
  scale_fill_manual(values = raindrop_colours$colour,
                    guide = "none", aesthetics = c("fill", "colour")) +
  labs(tag = "c)", x = "Precipitation treatment",
       y = "Simpson's index (D)") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

# treatment NMDS plot
fig2d <- rd_scores %>%
  ggplot(aes(x = MDS1, y = MDS2, colour = treatment, fill = treatment)) +
  stat_ellipse(aes(colour = NULL), geom = "polygon", alpha = 0.4, level = 0.80, show.legend = F) +
  geom_point(alpha = 0.95, size = 7) +
  scale_colour_manual(values = raindrop_colours$colour, aesthetics = c("colour", "fill"),
                      labels = c("Ambient\nControl", "Procedural\nControl",
                                 "Drought", "Irrigated")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Precipitation\ntreatment",
       fill = "Precipitation\ntreatment", tag = "d)") +
  guides(colour = guide_legend(keyheight = 3)) +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank())

layout <- "
AADDD
BBDDD
CCDDD
"

ggsave(fig2a + fig2b + fig2c + fig2d + plot_layout(design = layout),
       filename = "output/manuscript_figures/figure2.jpeg",
       width = 38, height = 23, units = "cm", dpi = 1500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Figure 3 - Community composition robust but changing through time ####

load("data/NMDS_model_results.RData", verbose = TRUE)
load("data/richness_models.RData")

## fig 2d. Richness increase over time
preddat_richness <- expand_grid(year_s = seq(-1.15,1.6, length.out = 100),
                                treatment = unique(diversity$treatment),
                                month = "June",
                                block = unique(diversity$block),
                                year_f = unique(diversity$year_f),
                                value = NA)

richness_pred <-  brms::posterior_predict(richness_year_linear,
                                          newdata = preddat_richness)

# summary data for each year value
predatsum <- preddat_richness %>%
  group_by(year_s) %>%
  summarise(value = NA, upr = NA, lwr = NA) %>%
  ungroup() %>% as.data.frame()

for(i in 1:nrow(predatsum)){
  cpos = which(preddat_richness$year_s == predatsum[i,]$year_s)

  post_mean = mean(richness_pred[,cpos])
  cQuant = quantile(richness_pred[,cpos], c(0.1, 0.90)) # 80% levels

  predatsum[i,]$value <- post_mean
  predatsum[i,]$upr <- cQuant[2]
  predatsum[i,]$lwr <- cQuant[1]

}

fig3a <- ggplot(predatsum, aes(x = year_s, y = value)) +
  geom_jitter(data = diversity, aes(y = richness),
              width = 0.04, size = 3, alpha = 0.9,
              colour = "lightsteelblue4") +
  geom_smooth(stat = "identity", aes(ymax = upr, ymin = lwr), alpha = 0.1,
              colour = "black", fill = "lightsteelblue4") +
  scale_y_continuous(breaks = seq(0,50, by = 5)) +
  scale_x_continuous(breaks = c(-1.18, -0.63, -0.0788, 0.473, 1.02, 1.58),
                     labels = 2016:2021) +
  labs(tag = "a)", x = "Year", y = "Species richness") +
  theme_bw(base_size = 17) +
  theme(panel.grid = element_blank())

# year NMDS plot
fig3b <- rd_scores %>%
  ggplot(aes(x = MDS1, y = MDS2, colour = year_f, fill = year_f)) +
  stat_ellipse(aes(colour = NULL), geom = "polygon", alpha = 0.2, level = 0.80, show.legend = F) +
  geom_point(alpha = 0.9, size = 3) +
  scale_colour_viridis_d(begin = 0.1, end = 0.9,
                         aesthetics = c("colour", "fill")) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Year", fill = "Year",
       tag = "b)") +
  theme_bw(base_size = 17) +
  theme(panel.grid = element_blank())

## fig 3c. NMDS1 increase over time
preddat_nmds1 <- expand_grid(year_s = seq(-1.5,1.5, length.out = 100),
                             treatment = unique(rd_scores$treatment),
                             block = unique(rd_scores$block),
                             value = NA)

nmds1_pred <-  brms::posterior_predict(nmds1_year,
                                          newdata = preddat_nmds1)

# summary data for each year value
predatsum_nmds <- preddat_nmds1 %>%
  group_by(year_s) %>%
  summarise(value = NA, upr = NA, lwr = NA) %>%
  ungroup() %>% as.data.frame()

for(i in 1:nrow(predatsum_nmds)){
  cpos = which(preddat_nmds1$year_s == predatsum_nmds[i,]$year_s)

  post_mean = mean(nmds1_pred[,cpos])
  cQuant = quantile(nmds1_pred[,cpos], c(0.1, 0.90)) # 80% levels

  predatsum_nmds[i,]$value <- post_mean
  predatsum_nmds[i,]$upr <- cQuant[2]
  predatsum_nmds[i,]$lwr <- cQuant[1]

}

fig3c <- ggplot(predatsum_nmds, aes(x = year_s, y = value, colour = year_f)) +
  geom_jitter(data = rd_scores, aes(y = MDS1),
              width = 0.04, size = 3, alpha = 0.8) +
  geom_smooth(stat = "identity", aes(ymax = upr, ymin = lwr,
                                     colour = NULL), alpha = 0.3,
              colour = "black", fill = "lightsteelblue") +
  scale_x_continuous(breaks = c(-1.4577380, -0.8746428, -0.2915476,
                                0.2915476,  0.8746428, 1.4577380),
                     labels = 2016:2021) +
  scale_colour_viridis_d(begin = 0.1, end = 0.9, guide = "none") +
  labs(tag = "c)", x = "Year", y = "NMDS1") +
  theme_bw(base_size = 17) +
  theme(panel.grid = element_blank())

ggsave(fig3a + fig3b + fig3c,
       filename = "output/manuscript_figures/figure3.jpeg",
       width = 40, height = 12, units = "cm", dpi = 1500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Figure 4 - Species contributions to community dissimilarity over years ####

load("output/SIMPER_results.RData", verbose = TRUE)

pc_labs <- percent_cover %>%
  mutate(species = gsub("_", " ", species)) %>%
  group_by(species) %>%
  summarise(group = group[1]) %>% ungroup()

# contribution percentage to dissimilarity between years
yc <- year_cont %>%
  left_join(x  = ., y = pc_labs, by = "species")
yc[16,"group"] <- "Grass"

fig4a <- ggplot(yc, aes(y = reorder(species, -mean_contribution),
                                        size = n_comparisons, colour = group)) +
  geom_segment(aes(x = mean_contribution, xend = median_contribution,
                   yend = species, size = NULL), show.legend = F, size = 1.4) +
  geom_point(aes(x = mean_contribution, shape = "Mean")) +
  geom_point(aes(x = median_contribution, shape = "Median"), show.legend = F) +
  labs(x = "Average contribution to community\ndissimilarity between years (%)",
       y = NULL, shape = NULL, colour = NULL,
       size = "Number of\ncommunity\ncomparisons", tag = "a)") +
  scale_colour_viridis_d(begin = 0.2, end = 0.8, option = "C",
                         labels = c("Forb", "Graminoid", "Legume")) +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         colour = guide_legend(override.aes = list(size = 5))) +
  theme_bw(base_size = 21) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "italic"))

# most influential species time series
fig4b <- top_species_perc %>%
  mutate(species = factor(species,
                          levels = c("Arrhenatherum elatius",
                                     "Lotus corniculatus",
                                     "Brachypodium pinnatum"))) %>%
  ggplot(aes(x = year, y = proportion,
                               colour = treatment, group = plot)) +
  geom_line(size = 1, show.legend = F) +
  geom_point(size = 3) +
  facet_wrap(~ species, ncol = 1) +
  scale_colour_manual(values = raindrop_colours$colour,
                      labels = c("Ambient\nControl", "Procedural\nControl",
                                 "Drought", "Irrigated")) +
  labs(x = "Year", y = "Relative abundance, p", colour = "Treatment",
       tag = "b)") +
  theme_bw(base_size = 21) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "italic"),
        strip.background = element_blank())

ggsave(fig4a + fig4b + plot_layout(widths = c(5,4)),
       filename = "output/manuscript_figures/figure4_simple.jpeg",
       width = 43, height = 21, units = "cm", dpi = 1500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Accidental aRt ####

# accidental_Rt <- rd_scores %>%
#   ggplot(aes(x = MDS1, y = MDS2, colour = treatment, fill = treatment)) +
#   stat_ellipse(aes(colour = NULL), geom = "polygon", alpha = 0.3, level = 0.95, show.legend = F) +
#   stat_ellipse(aes(colour = NULL), geom = "polygon", alpha = 0.4, level = 0.75, show.legend = F) +
#   stat_ellipse(aes(colour = NULL), geom = "polygon", alpha = 0.5, level = 0.50, show.legend = F) +
#   stat_ellipse(aes(colour = NULL), geom = "polygon", alpha = 0.6, level = 0.25, show.legend = F) +
#   scale_colour_manual(values = raindrop_colours$colour, aesthetics = c("colour", "fill")) +
#   theme_void() +
#   theme(plot.background = element_rect(fill = "white"))
# 
# ggsave(accidental_Rt, filename = "../accidental_aRt.pdf",
#        width = 20, height = 20, units = "cm")
# 
# ggsave(accidental_Rt, filename = "../accidental_aRt.jpeg",
#        width = 10, height = 10, units = "cm", dpi = 1200)


