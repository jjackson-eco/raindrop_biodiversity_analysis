####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##          Data cleaning and aquisition          ##
##                                                ##
##                 March 2022                     ##
##                                                ##
####################################################

## Loading, cleaning and merging broad biodiversity data from raindrop 2016-2020

rm(list = ls())
options(width = 100)

library(tidyverse)
library(patchwork)
library(viridis)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load data ####

## 1a. Above Ground Biomass
biomass <- read_csv("../../RainDropRobotics/Data/Raindrop_biomass_2016-20.csv")
glimpse(biomass)

## 1b. Cover
cover16_18 <- read_csv("../../RainDropRobotics/Data/Raindrop_(Oxford)_Cover_Data_2016-18.csv")
glimpse(cover16_18)
cover19_20 <- read_csv("../../RainDropRobotics/Data/RainDrop_cover_2019_2020.csv")
glimpse(cover19_20)

## 1c. 2021 data + some wrangling
load("../../DROUGHTNet/percent_cover_2021.RData")
percentcover_2021 <- percentcover_2021 %>%
  mutate(year = 2021, month = "June", species = gsub(" ", "_", species_name),
         species_level = if_else(species_level == 1, "Yes", "No"))

biomass_2021 <- read_csv("../../DROUGHTNet/droughtnet_biomass_july2021.csv") %>%
  mutate(group = case_when(
    `Functional Group` == "Grass" ~ "Graminoids",
    `Functional Group` == "Forb" ~ "Forbs",
    `Functional Group` == "Legume" ~ "Legumes",
    `Functional Group` == "Moss" ~ "Bryophytes"),
    Treatment = if_else(Treatment == "Procedural", "Control", Treatment),
    year = 2021, harvest = "Mid", biomass_g = `Dry Biomass (g)`,
    biomass_g = if_else(is.na(biomass_g) == TRUE, 0, biomass_g)) %>%
  rename_all(.funs = tolower) %>%
  dplyr::select(year, harvest, block, treatment, group, biomass_g)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Percentage Cover 2019/2020 ####

## 2a. unique names for plant types - 1 type per species (that has a type listed)
cover_planttypes <- cover16_18 %>%
  group_by(Species) %>%
  summarise(n_types = n_distinct(Plant_Type),
            n_na_types = length(which(is.na(unique(Plant_Type)))),
            Plant_Type = na.omit(object = unique(Plant_Type))[1]) %>%
  dplyr::select(Species, Plant_Type) %>%
  drop_na(Plant_Type) %>%
  bind_rows(.,
            tibble(Species = c("Odontites_verna", "Taraxacum_sect._vulgaria",
                               "Vicia_cracca", "Sonchus_oleraceus",
                               "Anacamptis_pyramidalis", "Orobanche_minor",
                               "Convolvolus_arvensis", "Hypericum_hirsutum",
                               "Senecio_jacobaea", "Phleum_bertinolii",
                               "Ranunculus_bulbosus"),
                   Plant_Type = c("Forb", "Forb", "Legume", "Forb",
                                  "Forb", "Forb", "Forb", "Forb", "Forb",
                                  "Grass", "Forb")))

# Check/ change those with more than 1 plant type - needs verification from AH
cover_planttypes[cover_planttypes$Species == "Trifolium_repens","Plant_Type"]
cover_planttypes[cover_planttypes$Species ==
                   "Tragopogon_pratensis","Plant_Type"] <- "Forb"

## 2b. tidying the cover 19-20 data
cover19_20_tidy <- cover19_20 %>%
  # date column conistent
  dplyr::rename(date = Year) %>%
  # remove oddities from species columns
  mutate(Species = gsub("[.]", "", x = gsub(" ", "_", Species))) %>%
  # long format
  pivot_longer(-c(date,Species), names_to = "Plot", values_to = "Percent_Cover") %>%
  # changing the plot column + adding treatment and block
  mutate(Plot = gsub("_Control", "_Ambient", Plot),
         Plot = gsub("_P. Control", "_Control", Plot),
         Plot = gsub("_Water", "_Irrigated", Plot),
         treatment = gsub("._", "", Plot),
         block = gsub("_.*", "", Plot)) %>%
  # adding plant type from previous cover data
  left_join(x = ., y = cover_planttypes, by = "Species") %>%
  # remove na cover values
  drop_na(Percent_Cover) %>%
  # year and month
  mutate(year = as.numeric(gsub(".*_", "", date)),
         month = gsub("_.*", "", date)) %>%
  # correct columns
  dplyr::select(year, month, block, treatment,
                species = Species, group = Plant_Type,
                percent_cover = Percent_Cover) %>%
  # adding whether info available at species level
  mutate(species_level =
           if_else(species %in% c("Prunus_sp","Quercus_sp",
                                  "Rosa_sp","Bareground","Bryophytes") == T,
         "No", "Yes")) %>%
  dplyr::select(c(1:4,8,5:7))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Percentage Cover 2016-2018 + Combine ####

cover16_18_tidy <- cover16_18 %>%
  # treatment changes
  mutate(Treatment = case_when(
    Treatment == "PC" ~ "Control",
    Treatment == "C" ~ "Ambient",
    Treatment == "D" ~ "Drought",
    Treatment == "I" ~ "Irrigated"
  )) %>%
  # year and month
  mutate(year = as.numeric(gsub(".*_", "", date)),
         month = gsub("_.*", "", date)) %>%
  # dealing with species column inconsistencies
  mutate(Species = gsub("[.]sp", "_sp", Species),
         Species = gsub("sp[.]", "sp", Species),
         Species = gsub("_species", "_sp", Species),
         Species = gsub("_$", "_sp", Species),
         Species = gsub(".praecox.subspecies.P", "_praecox_p", Species),
         Species = gsub("Tragopogon$", "Tragopogon_sp", Species)) %>%
  # info available at the species level
  mutate(species_level = if_else(grepl("_sp$", Species) == T, "No", "Yes")) %>%
  # adding plant type from cleaned version (just one type per species and no NA)
  dplyr::select(-Plant_Type) %>%
  left_join(x = ., y = cover_planttypes, by = "Species") %>%
  # correct columns
  dplyr::select(year, month, block = Block, treatment = Treatment,
                species_level, species = Species, group = Plant_Type,
                percent_cover = Percent_Cover)

## combine data -16-20
percent_cover <- bind_rows(cover16_18_tidy, cover19_20_tidy)

## adding 2021
percentcover_2021 <- percentcover_2021 %>%
  left_join(x = ., y = cover_planttypes, by = c("species" = "Species")) %>%
  mutate(species_level = if_else(species == "Quercus_sp", "No", species_level)) %>%
  dplyr::select(year, month, block, treatment, species_level,
                species, group = Plant_Type, percent_cover)

filter(percentcover_2021, is.na(group)) %>% as.data.frame() ## all non species-level

percent_cover <- bind_rows(percent_cover, percentcover_2021)
biomass <- bind_rows(biomass, biomass_2021)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Save ####

### Full workflow here but data itself available on request
save(biomass, percent_cover, file = "../../RainDropRobotics/Data/raindrop_biodiversity_2016_2021.RData")



