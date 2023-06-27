####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##   Total biomass models - Gamma distributions   ##
##                                                ##
##                Mar 7th 2022                    ##
##                                                ##
####################################################

## Building Bayesian heirarchical models of total biodiversity indices with respect to
## the rainfall treatments, including block and annual effects. Simple effects first.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(tidybayes)
library(brms)
library(patchwork)
library(viridis)
library(gghalves)
library(MASS)

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
         year_s = as.numeric(scale(year)),
         log_tot_biomass = if_else(log_tot_biomass == 0, 0.0001, log_tot_biomass))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Total biomass models ####

## Simpler model - year as random effect
set.seed(666)
totbiomass_base <- brm(
  log_tot_biomass ~ 1 + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = Gamma(link = "log"),
  prior = c(
    prior(normal(1.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(9), class = sd, group = "block"),
    prior(exponential(9), class = sd, group = "block:treatment"),
    prior(exponential(9), class = sd, group = "year_f"),
    prior(gamma(5,0.8), class = shape)),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Treatment effect
set.seed(666)
totbiomass_treatment <- brm(
  log_tot_biomass ~ 1 + treatment + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = Gamma(link = "log"),
  prior = c(
    prior(normal(1.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(9), class = sd, group = "block"),
    prior(exponential(9), class = sd, group = "block:treatment"),
    prior(exponential(9), class = sd, group = "year_f"),
    prior(gamma(5,0.8), class = shape)),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
totbiomass_year_linear <- brm(
  log_tot_biomass ~ 1 + treatment + harvest + year_s + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = Gamma(link = "log"),
  prior = c(
    prior(normal(1.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(9), class = sd, group = "block"),
    prior(exponential(9), class = sd, group = "block:treatment"),
    prior(exponential(9), class = sd, group = "year_f"),
    prior(gamma(5,0.8), class = shape)),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year interaction with harvest
set.seed(666)
totbiomass_year_harvest <- brm(
  log_tot_biomass ~ 1 + treatment + year_s*harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = Gamma(link = "log"),
  prior = c(
    prior(normal(1.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(9), class = sd, group = "block"),
    prior(exponential(9), class = sd, group = "block:treatment"),
    prior(exponential(9), class = sd, group = "year_f"),
    prior(gamma(5,0.8), class = shape)),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment interaction with year
set.seed(666)
totbiomass_year_treat <- brm(
  log_tot_biomass ~ 1 + treatment*year_s + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = Gamma(link = "log"),
  prior = c(
    prior(normal(1.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(9), class = sd, group = "block"),
    prior(exponential(9), class = sd, group = "block:treatment"),
    prior(exponential(9), class = sd, group = "year_f"),
    prior(gamma(5,0.8), class = shape)),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)


## Model comparisons
totbiomass_base <- add_criterion(totbiomass_base, criterion = c("loo","waic"))
totbiomass_treatment <- add_criterion(totbiomass_treatment, criterion = c("loo","waic"))
totbiomass_year_linear <- add_criterion(totbiomass_year_linear, criterion = c("loo","waic"))
totbiomass_year_harvest <- add_criterion(totbiomass_year_harvest, criterion = c("loo","waic"))
totbiomass_year_treat <- add_criterion(totbiomass_year_treat, criterion = c("loo","waic"))

mod_comp_totbiomass <- as.data.frame(loo_compare(totbiomass_base, totbiomass_treatment,
                                                 totbiomass_year_linear, totbiomass_year_treat,
                                                 totbiomass_year_harvest,
                                                 criterion = "loo"))
mod_comp_totbiomass

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Figure 1 - Biomass reduction ####

## 3a. Figure 1a. Total biomass
# colours
raindrop_colours <-
  tibble(treatment = c("Ambient", "Control", "Drought", "Irrigated"),
         num = rep(100,4),
         colour = c("#417A0B", "#BFBFBF", "#E69F00", "#6ECCFF"))

tbsum <- totbiomass %>% group_by(treatment) %>%
  summarise(mn = mean(log_tot_biomass)) %>%
  mutate(trt = c(0.8,1.8, 2.8,3.8))

set.seed(666)
figS2b <- as.data.frame(brms::posterior_predict(totbiomass_treatment)) %>%
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
  labs(tage = "b)", x = "Precipitation treatment",
       y = "Above-ground net primary productivity (ANPP)") +
  coord_cartesian(ylim = c(0,10)) +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())
figS2b

## Supplementary figure with distributions of raw data
figS2a <- ggplot(totbiomass, aes(x = log_tot_biomass, fill = treatment)) +
  geom_histogram(bins = 20, colour = "black", size = 0.2) +
  geom_vline(data = tbsum, aes(xintercept = mn)) +
  facet_wrap(~ treatment, ncol = 2) +
  scale_y_continuous(breaks = seq(0,10, by = 2)) +
  scale_fill_manual(values = raindrop_colours$colour, guide = "none") +
  labs(tag = "a)",
       x = "Above-ground net primary productivity (ANPP)",
       y = "Frequency") +
  theme_bw(base_size = 20) +
  theme(panel.grid = element_blank(), strip.background = element_blank())
figS2a

ggsave(figS2a + figS2b,
       filename = "output/manuscript_figures_jan2023/SI/figureS2_gamma.jpeg",
       width = 35, height = 20, units = "cm", dpi = 1500)

