####################################################
##                                                ##
##             RainDrop biodiversity              ##
##                                                ##
##           Community structure analysis         ##
##                                                ##
##               June 26th 2023                   ##
##                                                ##
####################################################

rm(list = ls())
options(width = 100)

## Change the .libPaths and R_LIBS_USER to the right thing if you're on a uni computer
if(Sys.info()["nodename"] == "BIO-W-LT-000083") {.libPaths("C:/Users/zool2541/R-4.1.1/library/")}

library(tidyverse)
library(brms)
library(vegan)
library(patchwork)
library(viridis)
library(GGally)

# colours
raindrop_colours <-
  tibble(treatment = c("Ambient", "Control", "Drought", "Irrigated"),
         num = rep(100,4),
         colour = c("#417A0B", "#BFBFBF", "#E69F00", "#6ECCFF"))

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 1. Load and wrangle data ####

# load data
load("../../RainDropRobotics/Data/raindrop_biodiversity_2016_2021.RData",
     verbose = TRUE)

# Saving species names to go through and check
species_names_raw <- percent_cover %>%
  filter(species_level == "Yes") %>%
  distinct(species) %>%
  arrange()

write_csv(species_names_raw, file = "../../RainDropRobotics/Data/species_names_raw.csv")
species_names_corrected <- read_csv(file = "../../RainDropRobotics/Data/species_names_corrected.csv")

percent_cover <- percent_cover %>%
  left_join(x =., y = species_names_corrected,
            by = "species") %>%
  mutate(species = if_else(is.na(species_correct) == F,
                           species_correct, species)) %>%
  filter(!note %in% "not_species") %>%
  filter(month == "June") %>%
  dplyr::select(-c(species_correct, note))

# data frame of the community matrix
species_community_df <- percent_cover %>%
  mutate(sample = paste(month, year, block, treatment, sep = "_")) %>%
  filter(species_level == "Yes" &
           is.na(percent_cover) == FALSE) %>%
  dplyr::select(sample, species, percent_cover) %>%
  pivot_wider(id_cols = sample, names_from = species, values_from = percent_cover,
              values_fill = 0)

# convert to a matrix
species_community_mat <- as.matrix(species_community_df[,2:95])
rownames(species_community_mat) <- species_community_df$sample

# sample information data for analyses
sample_info <- species_community_df %>%
  select(sample) %>%
  mutate(month = unlist(lapply(str_split(sample, "_"), `[[`, 1)),
         year = as.numeric(unlist(lapply(str_split(sample, "_"), `[[`, 2))),
         block = unlist(lapply(str_split(sample, "_"), `[[`, 3)),
         treatment = unlist(lapply(str_split(sample, "_"), `[[`, 4)))

# Bray distance matrix from Vegan
sp_dist <- vegdist(species_community_mat)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 2. Non-metric multidimensional scaling (NMDS) ####

# Three dimensions needed. Bray distances
raindrop_NMDS <- metaMDS(species_community_mat, k = 3, trymax = 1000)

# scores and species points
rd_scores <- as.data.frame(raindrop_NMDS$points)
rd_species <- as.data.frame(raindrop_NMDS$species)

# converting to a nice data frame
rd_scores <- rd_scores %>%
  mutate(plot = rownames(.)) %>%
  mutate(month = as.character(map(strsplit(plot, "_"), 1)),
         year = as.numeric(map(strsplit(plot, "_"), 2)),
         block = as.character(map(strsplit(plot, "_"), 3)),
         treatment = as.character(map(strsplit(plot, "_"), 4))) %>%
  mutate(block_treatment = paste0(block, "_", treatment),
         year_s = as.numeric(scale(year)),
         year_f = as.factor(year))

# and for species
species_types <- percent_cover %>%
  group_by(species) %>%
  summarise(group = group[1])

rd_species <- rd_species %>%
  mutate(species = rownames(.)) %>%
  left_join(x = ., y = species_types, by = "species")

# treatment plot
rd_scores %>%
  ggplot(aes(x = MDS1, y = MDS2, colour = treatment)) +
  geom_point(alpha = 0.8, size = 4) +
  stat_ellipse(level = 0.8, size = 1, show.legend = F) +
  scale_colour_manual(values = raindrop_colours$colour) +
  labs(x = "NMDS1", y = "NMDS2", colour = "Rainfall\ntreatment") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

save(rd_scores, rd_species, file = "data/NMDS_results.RData")

# individual
nmds1_t <- rd_scores %>%
  ggplot(aes(x = treatment, y = MDS1, colour = treatment)) +
  geom_violin(aes(fill = treatment), colour = "white",
              alpha = 0.4, trim = F, show.legend = F) +
  geom_jitter(alpha = 0.9, size = 3, width = 0.1) +
  scale_colour_manual(values = raindrop_colours$colour,
                      aesthetics = c("colour", "fill"), guide = "none") +
  labs(x = "Rainfall treatment", y = "NMDS1") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

nmds2_t <- rd_scores %>%
  ggplot(aes(x = treatment, y = MDS2, colour = treatment)) +
  geom_violin(aes(fill = treatment), colour = "white",
              alpha = 0.4, trim = F, show.legend = F) +
  geom_jitter(alpha = 0.9, size = 3, width = 0.1) +
  scale_colour_manual(values = raindrop_colours$colour,
                      aesthetics = c("colour", "fill"), guide = "none") +
  labs(x = "Rainfall treatment", y = "NMDS2") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

nmds3_t <- rd_scores %>%
  ggplot(aes(x = treatment, y = MDS3, colour = treatment)) +
  geom_violin(aes(fill = treatment), colour = "white",
              alpha = 0.4, trim = F, show.legend = F) +
  geom_jitter(alpha = 0.9, size = 3, width = 0.1) +
  scale_colour_manual(values = raindrop_colours$colour,
                      aesthetics = c("colour", "fill"), guide = "none") +
  labs(x = "Rainfall treatment", y = "NMDS3") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank())

ggsave(nmds1_t / nmds2_t / nmds3_t, filename = "output/NMDS_treatments.jpeg",
       width = 20, height = 30, units = "cm", dpi = 1500)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 3. Temporal trends in NMDS scores ####

#_______________________________________________________________
## NMDS1

ggplot(rd_scores, aes(x = year_f, y = MDS1)) + geom_boxplot()

set.seed(666)
nmds1_base <- brm(
  MDS1 ~ 1 + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(exponential(5), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
nmds1_year <- brm(
  MDS1 ~ 1 + year_s + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(5), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
nmds1_year_rs <- brm(
  MDS1 ~ 1 + year_s + (1 + year_s|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(5), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
nmds1_year_blocktreatment <- brm(
  MDS1 ~ 1 + year_s*block_treatment + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(5), class = sd, group = "block"),
    prior(exponential(5), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model selection
nmds1_base <- add_criterion(nmds1_base, criterion = c("loo","waic"))
nmds1_year <- add_criterion(nmds1_year, criterion = c("loo","waic"))
nmds1_year_rs <- add_criterion(nmds1_year_rs, criterion = c("loo","waic"))
nmds1_year_blocktreatment <- add_criterion(nmds1_year_blocktreatment, criterion = c("loo","waic"))

mod_comp_nmds1 <- as.data.frame(loo_compare(nmds1_base, nmds1_year,
                                            nmds1_year_rs,
                                            nmds1_year_blocktreatment,
                                            criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(nmds1_year$formula)[1],
                     as.character(nmds1_year_rs$formula)[1],
                     as.character(nmds1_year_blocktreatment$formula)[1],
                     as.character(nmds1_base$formula)[1]),
         formula = gsub("month", "harvest", formula))
mod_comp_nmds1

save(mod_comp_nmds1, nmds1_year, file = "data/NMDS_model_results.RData")

#_______________________________________________________________
## NMDS2

ggplot(rd_scores, aes(x = year_f, y = MDS2)) + geom_boxplot()

set.seed(666)
nmds2_base <- brm(
  MDS2 ~ 1 + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(exponential(1), class = sd, group = "block"),
    prior(exponential(1), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
nmds2_year <- brm(
  MDS2 ~ 1 + year_s + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(1), class = sd, group = "block"),
    prior(exponential(1), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.97),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
nmds2_yearf <- brm(
  MDS2 ~ 1 + year_f + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(1), class = sd, group = "block"),
    prior(exponential(1), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)


set.seed(666)
nmds2_year_blocktreatment <- brm(
  MDS2 ~ 1 + year_s*block_treatment + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(1), class = sd, group = "block"),
    prior(exponential(1), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model selection
nmds2_base <- add_criterion(nmds2_base, criterion = c("loo","waic"))
nmds2_year <- add_criterion(nmds2_year, criterion = c("loo","waic"))
nmds2_yearf <- add_criterion(nmds2_yearf, criterion = c("loo","waic"))
nmds2_year_blocktreatment <- add_criterion(nmds2_year_blocktreatment, criterion = c("loo","waic"))

mod_comp_nmds2 <- as.data.frame(loo_compare(nmds2_base, nmds2_year,
                                            nmds2_yearf, nmds2_year_blocktreatment,
                                            criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(nmds2_base$formula)[1],
                     as.character(nmds2_yearf$formula)[1],
                     as.character(nmds2_year$formula)[1],
                     as.character(nmds2_year_blocktreatment$formula)[1]),
         formula = gsub("month", "harvest", formula))
mod_comp_nmds2

#_______________________________________________________________
## NMDS3

ggplot(rd_scores, aes(x = year_f, y = MDS3)) + geom_boxplot()

set.seed(666)
nmds3_base <- brm(
  MDS3 ~ 1 + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(exponential(2), class = sd, group = "block"),
    prior(exponential(2), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
nmds3_year <- brm(
  MDS3 ~ 1 + year_s + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(2), class = sd, group = "block"),
    prior(exponential(2), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

set.seed(666)
nmds3_yearf <- brm(
  MDS3 ~ 1 + year_f + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(2), class = sd, group = "block"),
    prior(exponential(2), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)


set.seed(666)
nmds3_year_blocktreatment <- brm(
  MDS3 ~ 1 + year_s*block_treatment + (1|block/treatment),
  data = rd_scores, family = gaussian(),
  prior = c(
    prior(normal(0, 0.1), class =  Intercept),
    prior(normal(0, 0.1), class = b), # all beta terms
    prior(exponential(2), class = sd, group = "block"),
    prior(exponential(2), class = sd, group = "block:treatment")),
  control = list(adapt_delta = 0.96),
  chains = 4, cores = 4, iter = 4000, warmup = 2000
)

## Model selection
nmds3_base <- add_criterion(nmds3_base, criterion = c("loo","waic"))
nmds3_year <- add_criterion(nmds3_year, criterion = c("loo","waic"))
nmds3_yearf <- add_criterion(nmds3_yearf, criterion = c("loo","waic"))
nmds3_year_blocktreatment <- add_criterion(nmds3_year_blocktreatment, criterion = c("loo","waic"))

mod_comp_nmds3 <- as.data.frame(loo_compare(nmds3_base, nmds3_year,
                                            nmds3_yearf, nmds3_year_blocktreatment,
                                            criterion = "loo")) %>%
  mutate(model_name = rownames(.),
         formula = c(as.character(nmds3_yearf$formula)[1],
                     as.character(nmds3_year$formula)[1],
                     as.character(nmds3_base$formula)[1],
                     as.character(nmds3_year_blocktreatment$formula)[1]),
         formula = gsub("month", "harvest", formula))
mod_comp_nmds3

# all model comparison tables
save(mod_comp_nmds1, mod_comp_nmds2, mod_comp_nmds3, file = "data/NMDS_model_comparisons.RData")

#_______________________________________________________________
## Species contributing most to NMDS1

rd_species %>%
  mutate(abs_MDS1 = abs(MDS1)) %>%
  arrange(desc(abs_MDS1))

percent_cover %>%
  filter(species == "Cirsium_vulgare")

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 4. Analysis of similarity ####

# Treatments
anosim_treatment <- with(sample_info, anosim(sp_dist, treatment))
summary(anosim_treatment)

# Year
anosim_year <- with(sample_info, anosim(sp_dist, year))
summary(anosim_year)

# block
anosim_block <- with(sample_info, anosim(sp_dist, block))
summary(anosim_block)

#_______________________________________________________________
## Plots
anosim_treatment_plot <-
  tibble(group = anosim_treatment$class.vec,
       dis_rank = anosim_treatment$dis.rank) %>%
  ggplot(aes(x = group, y = dis_rank, fill = group)) +
  geom_violin() +
  geom_jitter(width = 0.09, alpha = 0.1) +
  geom_boxplot(outlier.shape = NA, coef = 0,
               varwidth = 0.5) +
  scale_fill_manual(values = c("white", raindrop_colours$colour),
                    guide = "none") +
  labs(x = "Treatment", y = "Bray dissimilarity rank") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

anosim_year_plot <-
  tibble(group = anosim_year$class.vec,
         dis_rank = anosim_year$dis.rank) %>%
  ggplot(aes(x = group, y = dis_rank, fill = group)) +
  geom_violin() +
  geom_jitter(width = 0.02, alpha = 0.05) +
  geom_boxplot(outlier.shape = NA, coef = 0,
               varwidth = 0.5) +
  scale_fill_manual(values = c("white",
                               viridis(6, begin = 0.1, end = 0.9)),
                    guide = "none") +
  labs(x = "Year", y = "Bray dissimilarity rank") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

anosim_block_plot <-
  tibble(group = anosim_block$class.vec,
         dis_rank = anosim_block$dis.rank) %>%
  ggplot(aes(x = group, y = dis_rank, fill = group)) +
  geom_violin(show.legend = F) +
  geom_jitter(width = 0.02, alpha = 0.05, show.legend = F) +
  geom_boxplot(outlier.shape = NA, coef = 0,
               varwidth = 0.5, show.legend = F) +
  labs(x = "Block", y = "Bray dissimilarity rank") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank())

ggsave(anosim_treatment_plot/anosim_year_plot/anosim_block_plot,
       filename = "output/anosim_disrank.jpeg",
       width = 14, height = 25, units = "cm",
       dpi = 1000)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 5. Multivariate Analysis of Variance ####

adonis2(sp_dist ~ treatment*block, data = sample_info)
adonis2(sp_dist ~ block, data = sample_info)
adonis2(sp_dist ~ treatment, data = sample_info)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 6. Similarity percentage for species ####

simper_treatment <- with(sample_info, simper(species_community_mat, treatment))
simper_year <- with(sample_info, simper(species_community_mat, year))
simper_block <- with(sample_info, simper(species_community_mat, block))

#_______________________________________________________________
## Extracting most influential species

# Treatment
treatment_cont <- bind_rows(lapply(simper_treatment, function(x){

  av_con = x$average[order(x$average, decreasing = TRUE)]
  spp = names(av_con)

  return(tibble(species = spp[1:10], contribution = av_con[1:10]))})) %>%
  mutate(comparison = rep(names(simper_treatment), each = 10)) %>%
  group_by(species) %>%
  summarise(mean_contribution = mean(contribution)*100,
            median_contribution = median(contribution)*100,
            n_comparisons = n()) %>%
  ungroup() %>%
  mutate(species = gsub(x = species, "_", " "))

# year
year_cont <- bind_rows(lapply(simper_year, function(x){

  av_con = x$average[order(x$average, decreasing = TRUE)]
  spp = names(av_con)

  return(tibble(species = spp[1:10], contribution = av_con[1:10]))})) %>%
  mutate(comparison = rep(names(simper_year), each = 10)) %>%
  group_by(species) %>%
  summarise(mean_contribution = mean(contribution)*100,
            median_contribution = median(contribution)*100,
            n_comparisons = n()) %>%
  ungroup() %>%
  mutate(species = gsub(x = species, "_", " "))

# block
block_cont <- bind_rows(lapply(simper_block, function(x){

  av_con = x$average[order(x$average, decreasing = TRUE)]
  spp = names(av_con)

  return(tibble(species = spp[1:10], contribution = av_con[1:10]))})) %>%
  mutate(comparison = rep(names(simper_block), each = 10)) %>%
  filter(comparison %in% comparison[grep("C", comparison)]) %>%
  group_by(species) %>%
  summarise(mean_contribution = mean(contribution)*100,
            median_contribution = median(contribution)*100,
            n_comparisons = n()) %>%
  ungroup() %>%
  mutate(species = gsub(x = species, "_", " "))

#_______________________________________________________________
## Plots
treatment_cont_plot <- ggplot(treatment_cont, aes(y = reorder(species, -mean_contribution),
                           size = n_comparisons)) +
  geom_segment(aes(x = mean_contribution, xend = median_contribution,
                   yend = species, size = NULL), show.legend = F, size = 1.4) +
  geom_point(aes(x = mean_contribution, shape = "Mean")) +
  geom_point(aes(x = median_contribution, shape = "Median"), show.legend = F) +
  labs(x = "Average contribution to community\ndissimilarity between rainfall treatments (%)",
       y = NULL, shape = NULL, size = "Number of\ncommunity\ncomparisons", tag = "a)") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "italic"))

year_cont_plot <- ggplot(year_cont, aes(y = reorder(species, -mean_contribution),
                                                  size = n_comparisons)) +
  geom_segment(aes(x = mean_contribution, xend = median_contribution,
                   yend = species, size = NULL), show.legend = F, size = 1.4) +
  geom_point(aes(x = mean_contribution, shape = "Mean")) +
  geom_point(aes(x = median_contribution, shape = "Median"), show.legend = F) +
  labs(x = "Average contribution to community\ndissimilarity between years (%)",
       y = NULL, shape = NULL, size = "Number of\ncommunity\ncomparisons", tag = "b)") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "italic"))

block_cont_plot <- ggplot(block_cont, aes(y = reorder(species, -mean_contribution),
                                                  size = n_comparisons)) +
  geom_segment(aes(x = mean_contribution, xend = median_contribution,
                   yend = species, size = NULL), show.legend = F, size = 1.4) +
  geom_point(aes(x = mean_contribution, shape = "Mean")) +
  geom_point(aes(x = median_contribution, shape = "Median"),
             show.legend = F) +
  labs(x = "Average contribution to community\ndissimilarity for block C (%)",
       y = NULL, size = "Number of\ncommunity\ncomparisons", tag = "c)",
       shape = NULL) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "italic"))

ggsave(treatment_cont_plot + year_cont_plot + block_cont_plot,
       filename = "output/community_contributions_by_species.jpeg",
       width = 52, height = 18, units = "cm", dpi = 1000)

ggsave(year_cont_plot +
         labs(x = "Average contribution to community dissimilarity between years (%)",
              tag = NULL),
       filename = "output/community_contributions_years.jpeg",
       width = 26, height = 18, units = "cm", dpi = 1000)

#_______________________________________________________________
## Year top four species

# year top contributors
year_top_spp <- year_cont %>%
  filter(mean_contribution > 5) %>%
  mutate(species = gsub(" ", "_", species)) %>%
  pull(species)

# percent cover data per plot
top_species_perc <- percent_cover %>%
  filter(species %in% year_top_spp) %>%
  mutate(plot = paste0(block, "_", treatment),
         species = gsub("_", " ", species)) %>%
  group_by(plot) %>%
  # first convert percentages to a proportion of total i.e. relative abundance
  mutate(proportion = percent_cover/sum(percent_cover)) %>%
  ungroup()

save(simper_year, year_cont, top_species_perc, file = "output/SIMPER_results.RData")

top_species_plot <-
  ggplot(top_species_perc, aes(x = year, y = percent_cover,
                             colour = treatment, group = plot)) +
  geom_line(size = 1, show.legend = F) +
  geom_point(size = 3) +
  facet_wrap(~ species) +
  scale_y_continuous(limits = c(0,100)) +
  scale_colour_manual(values = raindrop_colours$colour) +
  labs(x = "Year", y = "Percentage cover", colour = "Treatment",
       tag = "b)") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "italic"),
        strip.background = element_blank())

ggsave(top_species_plot + labs(tag = NULL),
       filename = "output/top_species_percent_cover.jpeg",
       width = 22, height = 17, units = "cm", dpi = 1000)


ggsave((year_cont_plot + labs(tag = "a)")) + (top_species_plot),
       filename = "output/top_species.jpeg",
       width = 45, height = 17, units = "cm", dpi = 1000)

##____________________________________________________________________________________________________________________________________________________________________________________________________________
#### 7. Looking at most abundant species ####

species_rankings <- percent_cover %>%
  filter(species_level == "Yes" &
           is.na(percent_cover) == FALSE) %>%
  group_by(species) %>%
  summarise(n_obs = n(),
            max_cov = max(percent_cover),
            median_cov = median(percent_cover),
            min_cover = min(percent_cover))

# Are these associated with each other?
ggpairs(species_rankings[,-1])

# Get the most abundant species
sp_highest_coverage <- species_rankings %>%
  filter(n_obs >= 3) %>%
  arrange(-median_cov) %>%
  pull(species) %>%
  .[1:10]

# Get the most widespread species
sp_most_obs <- species_rankings %>%
  arrange(-n_obs) %>%
  pull(species) %>%
  .[1:10]


# cover data for most abundant species
most_abundant_species <- percent_cover %>%
  filter(species_level == "Yes" &
           is.na(percent_cover) == FALSE &
           species %in% sp_highest_coverage == TRUE) %>%
  mutate(plot = paste0(block,"_",treatment)) %>%
  group_by(plot) %>%
  # first convert percentages to a proportion of total i.e. relative abundance
  mutate(proportion = percent_cover/sum(percent_cover),
         species = gsub("_", " ", species)) %>%
  ungroup()

plot_most_abundant <- ggplot(most_abundant_species, aes(x = year, y = proportion,
                                  colour = treatment, group = plot)) +
  geom_line(size = 1) +
  facet_wrap(~ species, ncol = 5) +
  scale_colour_manual(values = raindrop_colours$colour,
                      labels = c("Ambient\nControl", "Procedural\nControl",
                                 "Drought", "Irrigated")) +
  labs(x = "Year", y = "Relative abudance, p", colour = "Treatment") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "italic"),
        strip.background = element_blank())

ggsave(plot_most_abundant, filename = "output/most_abundant_10_species.jpeg",
       width = 38, height = 14, units = "cm", dpi = 1200)



