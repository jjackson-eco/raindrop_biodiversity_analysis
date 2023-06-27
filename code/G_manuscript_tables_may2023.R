####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##             Supplementary tables               ##
##                                                ##
##                May 3rd 2023                    ##
##                                                ##
####################################################

rm(list = ls())
options(width = 100)

## Change the .libPaths and R_LIBS_USER to the right thing if you're on a uni computer
if(Sys.info()["nodename"] == "BIO-W-LT-000083") {.libPaths("C:/Users/zool2541/R-4.1.1/library/")}

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
         colour = c("#61D94E", "#BFBFBF", "#EB8344", "#6ECCFF"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load ####

load("data/totbiomass_models.RData", verbose = T)
load("data/biomass_models.RData", verbose = T)
load("data/biomass_var_models.RData", verbose = T)
load("data/shannon_models.RData", verbose = T)
load("data/simpsons_models.RData", verbose = T)
load("data/richness_models.RData", verbose = T)
load("data/NMDS_model_comparisons.RData", verbose = T)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Create flextables  ####

## rounding function w/ arguments for numeric columns
roundfunc <- function(x){round(x, digits = 2)}

## table S1 - total biomass
mod_comp_totbiomass %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS1.png")

## table S2 - group-level biomass
mod_comp_biomass %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS2.png")

## table S3 - temporal stability in biomass
mod_comp_biomass_var %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS3.png")

## table S4 - Shannon index
mod_comp_shannon %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS4.png")

## table S5 - Simpsons index
mod_comp_simpsons %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS5.png")

## table S6 - Richness
mod_comp_richness %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS6.png")

## table S7 - NMDS1
mod_comp_nmds1 %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS7.png")

## table S8 - NMDS2
mod_comp_nmds2 %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS8.png")

## table S9 - NMDS3
mod_comp_nmds3 %>%
  mutate(across(1:8, roundfunc)) %>%
  dplyr::select(`R object` = model_name, `Model formula` = formula, `LOO elpd` = elpd_loo,
                `LOO elpd error` = se_elpd_loo, `elpd difference` = elpd_diff,
                `elpd error difference` = se_diff, `LOO information criterion` = looic) %>%
  flextable(cwidth = 1.5) %>%
  width(10, j = 2) %>%
  save_as_image("output/manuscript_figures_jan2023/SI/tableS9.png")

