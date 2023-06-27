####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##           Total biodiversity models            ##
##                                                ##
##                  June 2023                     ##
##                                                ##
####################################################

## Building Bayesian heirarchical models of total biodiversity indices with respect to
## the rainfall treatments, including block and annual effects. Simple effects first.

rm(list = ls())
options(width = 100)

library(tidyverse)
library(brms)
library(patchwork)
library(viridis)

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

## biomass variance for each quadrat
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
  mutate(plot = paste0(block,"_",treatment)) %>%
  group_by(plot) %>%
  # first convert percentages to a proportion of total i.e. relative abundance
  mutate(proportion = percent_cover/sum(percent_cover)) %>%
  ungroup() %>%
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

#%>%
  # pivot_longer(cols = c(richness, simpsons, shannon),
  #              names_to = "index") %>%
  # mutate(index_lab = case_when(
  #   index == "richness" ~ "Species richness",
  #   index == "simpsons" ~ "Simpson's index",
  #   index == "shannon" ~ "Shannon-Wiener index"
  # ))

# # full species crib
# species_crib <- percent_cover %>%
#   group_by(group, species) %>%
#   summarise(n_records = n(),
#             cover = mean(percent_cover)) %>%
#   arrange(-cover) %>%
#   mutate(species = gsub("_", " ", species))
#
# write_csv(species_crib, file = "data/species_crib_2016_2021.csv")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Prior predictive simulation ####

plot(density(rexp(1000, rate = 10)), main = "Rate = 10")
lines(density(rexp(1000, rate = 0.2)), col = "green")
lines(density(rexp(1000, rate = 1)), col = "red")
lines(density(rexp(1000, rate = 5)), col = "blue")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Total biomass models ####

## Simpler model - year as random effect
set.seed(666)
totbiomass_base <- brm(
  log_tot_biomass ~ 1 + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
  )

## Treatment effect
set.seed(666)
totbiomass_treatment <- brm(
  log_tot_biomass ~ 1 + treatment + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.99),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
totbiomass_year_linear <- brm(
  log_tot_biomass ~ 1 + treatment + harvest + year_s + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year interaction with harvest
set.seed(666)
totbiomass_year_harvest <- brm(
  log_tot_biomass ~ 1 + treatment + year_s*harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment interaction with year
set.seed(666)
totbiomass_year_treat <- brm(
  log_tot_biomass ~ 1 + treatment*year_s + harvest + (1|block/treatment) + (1|year_f),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.99),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year as an autocorrelated term
set.seed(666)
totbiomass_year_auto <- brm(
  log_tot_biomass ~ 1 + treatment + ar(gr = year_s, p = 1) + harvest + (1|block/treatment),
  data = totbiomass, family = gaussian(),
  prior = c(
    prior(normal(3.5, 0.5), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
totbiomass_base <- add_criterion(totbiomass_base, criterion = c("loo","waic"))
totbiomass_treatment <- add_criterion(totbiomass_treatment, criterion = c("loo","waic"))
totbiomass_year_linear <- add_criterion(totbiomass_year_linear, criterion = c("loo","waic"))
totbiomass_year_harvest <- add_criterion(totbiomass_year_harvest, criterion = c("loo","waic"))
totbiomass_year_treat <- add_criterion(totbiomass_year_treat, criterion = c("loo","waic"))
totbiomass_year_auto <- add_criterion(totbiomass_year_auto, criterion = c("loo","waic"))

mod_comp_totbiomass <- as.data.frame(loo_compare(totbiomass_base, totbiomass_treatment,
                                      totbiomass_year_linear, totbiomass_year_treat,
                                      totbiomass_year_harvest, totbiomass_year_auto,
                                      criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(totbiomass_year_linear$formula)[1],
                     as.character(totbiomass_treatment$formula)[1],
                     as.character(totbiomass_year_harvest$formula)[1],
                     as.character(totbiomass_year_treat$formula)[1],
                     as.character(totbiomass_base$formula)[1],
                     as.character(totbiomass_year_auto$formula)[1]))
mod_comp_totbiomass

save(totbiomass_treatment, mod_comp_totbiomass, file = "data/totbiomass_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Group-level biomass ####

## Simpler model - year as random effect
set.seed(666)
biomass_base <- brm(
  biomass_log ~ 1 + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment effect
set.seed(666)
biomass_treatment <- brm(
  biomass_log ~ 1 + treatment + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## group effect
set.seed(666)
biomass_group <- brm(
  biomass_log ~ 1 + group + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment and group effects
set.seed(666)
biomass_treatment_group <- brm(
  biomass_log ~ 1 + treatment + group + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment and group and interaction
set.seed(666)
biomass_treatment_group_int <- brm(
  biomass_log ~ 1 + treatment*group + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## treatment, group and year(linear) interactions
set.seed(666)
biomass_treatment_group_year <- brm(
  biomass_log ~ 1 + treatment*group*year_s + harvest + (1|block/treatment) + (1|year_f),
  data = biomass_gr, family = gaussian(),
  prior = c(
    prior(normal(1, 1), class =  Intercept),
    prior(normal(0, 1), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(4), class = sd, group = "block:treatment"),
    prior(exponential(4), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
biomass_base <- add_criterion(biomass_base, criterion = c("loo","waic"))
biomass_treatment <- add_criterion(biomass_treatment, criterion = c("loo","waic"))
biomass_group <- add_criterion(biomass_group, criterion = c("loo","waic"))
biomass_treatment_group <- add_criterion(biomass_treatment_group, criterion = c("loo","waic"))
biomass_treatment_group_int <- add_criterion(biomass_treatment_group_int, criterion = c("loo","waic"))
biomass_treatment_group_year <- add_criterion(biomass_treatment_group_year, criterion = c("loo","waic"))

mod_comp_biomass <- as.data.frame(loo_compare(biomass_base,
                          biomass_treatment,
                          biomass_group,
                          biomass_treatment_group,
                          biomass_treatment_group_int,
                          biomass_treatment_group_year,
                          criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(biomass_treatment_group_int$formula)[1],
                     as.character(biomass_treatment_group_year$formula)[1],
                     as.character(biomass_treatment_group$formula)[1],
                     as.character(biomass_group$formula)[1],
                     as.character(biomass_treatment$formula)[1],
                     as.character(biomass_base$formula)[1]))
mod_comp_biomass

save(biomass_treatment_group_int, mod_comp_biomass, file = "data/biomass_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Biomass variance - Temporal stability ####

## Only unit of replication here is block (no treatment replicated within block), so just block-level variance
set.seed(666)
biomass_var_base <- brm(
  cv_biomass ~ 1 + (1|block),
  data = biomass_var, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(exponential(3), class = sd, group = "block")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
biomass_var_treat <- brm(
  cv_biomass ~ 1 + treatment + (1|block),
  data = biomass_var, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.3), class = b), # all beta terms
    prior(exponential(3), class = sd, group = "block")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
biomass_var_base <- add_criterion(biomass_var_base, criterion = c("loo","waic"))
biomass_var_treat <- add_criterion(biomass_var_treat, criterion = c("loo","waic"))

mod_comp_biomass_var <- as.data.frame(loo_compare(biomass_var_base,
                                              biomass_var_treat,
                                              criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(biomass_var_treat$formula)[1],
                     as.character(biomass_var_base$formula)[1]))

mod_comp_biomass_var

save(biomass_var_treat, mod_comp_biomass_var, file = "data/biomass_var_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Shannon-Wiener index across treatments ####

## Base model - year as a random effect
set.seed(666)
shannon_base <- brm(
  shannon ~ 1 + month + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Simple treatment effect
set.seed(666)
shannon_treatment <- brm(
  shannon ~ 1 + month + treatment + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
shannon_year_linear <- brm(
  shannon ~ 1 + month + treatment + year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year treatment interaction
set.seed(666)
shannon_treatment_year_int <- brm(
  shannon ~ 1 + month + treatment*year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year as an autocorrelated term
set.seed(666)
shannon_treatment_year_auto <- brm(
  shannon ~ 1 + month + treatment + ar(gr = year_s, p = 1) + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
shannon_base <- add_criterion(shannon_base, criterion = c("loo","waic"))
shannon_treatment <- add_criterion(shannon_treatment, criterion = c("loo","waic"))
shannon_year_linear <- add_criterion(shannon_year_linear, criterion = c("loo","waic"))
shannon_treatment_year_int <- add_criterion(shannon_treatment_year_int, criterion = c("loo","waic"))
shannon_treatment_year_auto <- add_criterion(shannon_treatment_year_auto, criterion = c("loo","waic"))

mod_comp_shannon <- as.data.frame(loo_compare(shannon_base, shannon_treatment,
                                              shannon_year_linear,
                                              shannon_treatment_year_int,
                                              shannon_treatment_year_auto,
                                              criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(shannon_base$formula)[1],
                     as.character(shannon_treatment_year_int$formula)[1],
                     as.character(shannon_treatment_year_auto$formula)[1],
                     as.character(shannon_treatment$formula)[1],
                     as.character(shannon_year_linear$formula)[1]),
         formula = gsub("month", "harvest", formula))

mod_comp_shannon

save(shannon_treatment_year_int, mod_comp_shannon, file = "data/shannon_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 7. Simpson's index across treatments ####

## Base model - year as a random effect
set.seed(666)
simpsons_base <- brm(
  simpsons ~ 1 + month + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Simple treatment effect
set.seed(666)
simpsons_treatment <- brm(
  simpsons ~ 1 + month + treatment + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
simpsons_year_linear <- brm(
  simpsons ~ 1 + month + treatment + year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.95),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year treatment interaction
set.seed(666)
simpsons_treatment_year_int <- brm(
  simpsons ~ 1 + month + treatment*year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year as an autocorrelated term
set.seed(666)
simpsons_treatment_year_auto <- brm(
  simpsons ~ 1 + month + treatment + ar(gr = year_s, p = 1) + (1|block/treatment) + (1|year_f),
  data = diversity, family = gaussian(),
  prior = c(
    prior(normal(0, 0.5), class =  Intercept),
    prior(normal(0, 0.5), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment"),
    prior(exponential(5), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
simpsons_base <- add_criterion(simpsons_base, criterion = c("loo","waic"))
simpsons_treatment <- add_criterion(simpsons_treatment, criterion = c("loo","waic"))
simpsons_year_linear <- add_criterion(simpsons_year_linear, criterion = c("loo","waic"))
simpsons_treatment_year_int <- add_criterion(simpsons_treatment_year_int, criterion = c("loo","waic"))
simpsons_treatment_year_auto <- add_criterion(simpsons_treatment_year_auto, criterion = c("loo","waic"))

mod_comp_simpsons <- as.data.frame(loo_compare(simpsons_base, simpsons_treatment,
                                              simpsons_year_linear,
                                              simpsons_treatment_year_int,
                                              simpsons_treatment_year_auto,
                                              criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(simpsons_treatment_year_int$formula)[1],
                     as.character(simpsons_base$formula)[1],
                     as.character(simpsons_treatment_year_auto$formula)[1],
                     as.character(simpsons_treatment$formula)[1],
                     as.character(simpsons_year_linear$formula)[1]),
         formula = gsub("month", "harvest", formula))

mod_comp_simpsons

save(simpsons_treatment_year_int, mod_comp_simpsons, file = "data/simpsons_models.RData")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 8. Species richness across treatments ####

# looking at a better prior for the intercept
hist(log(diversity$richness))
hist(rnorm(1000, 3,0.25))

## Base model - year as a random effect
set.seed(666)
richness_base <- brm(
  richness ~ 1 + month + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Simple treatment effect
set.seed(666)
richness_treatment <- brm(
  richness ~ 1 + month + treatment + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year linear effect
set.seed(666)
richness_year_linear <- brm(
  richness ~ 1 + month + treatment + year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year treatment interaction
set.seed(666)
richness_treatment_year_int <- brm(
  richness ~ 1 + month + treatment*year_s + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## year as an autocorrelated term
set.seed(666)
richness_treatment_year_auto <- brm(
  richness ~ 1 + month + treatment + ar(gr = year_s, p = 1) + (1|block/treatment) + (1|year_f),
  data = diversity, family = poisson(link = "log"),
  prior = c(
    prior(normal(3, 0.25), class =  Intercept),
    prior(normal(0, 0.25), class = b), # all beta terms
    prior(exponential(6), class = sd, group = "block"),
    prior(exponential(8), class = sd, group = "block:treatment"),
    prior(exponential(8), class = sd, group = "year_f")),
  control = list(adapt_delta = 0.98),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model comparisons
richness_base <- add_criterion(richness_base, criterion = c("loo","waic"))
richness_treatment <- add_criterion(richness_treatment, criterion = c("loo","waic"))
richness_year_linear <- add_criterion(richness_year_linear, criterion = c("loo","waic"))
richness_treatment_year_int <- add_criterion(richness_treatment_year_int, criterion = c("loo","waic"))
richness_treatment_year_auto <- add_criterion(richness_treatment_year_auto, criterion = c("loo","waic"))

mod_comp_richness <- as.data.frame(loo_compare(richness_base, richness_treatment,
                                              richness_year_linear,
                                              richness_treatment_year_int,
                                              richness_treatment_year_auto,
                                              criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(richness_year_linear$formula)[1],
                     as.character(richness_base$formula)[1],
                     as.character(richness_treatment$formula)[1],
                     as.character(richness_treatment_year_auto$formula)[1],
                     as.character(richness_treatment_year_int$formula)[1]),
         formula = gsub("month", "harvest", formula))

mod_comp_richness

save(richness_year_linear, mod_comp_richness, file = "data/richness_models.RData")


