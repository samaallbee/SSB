


#Title: Code supplement for Ecology submission
#Author: Samantha Allbee


#load libraries 
library(lmerTest)
library(emmeans)
library(readr)
library(factoextra)
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(goeveg)
library(plot3D)
library(vegan3d)
library(plotly)
library(cowplot)
library(gridExtra)
library(devtools)
library(abdiv)
library(multcomp)
library(tidyr)

#load data
emergence <- read_csv("seed bank/cleaned_final_MS.csv")

####FIGURE 2 // GLMMS ANALYSIS####

#richness across entire dataset:
species_richness <- emergence %>%
  group_by(microhabitat, site,sample) %>%
  summarise(species_richness = n_distinct(species)) %>%
  ungroup()
#model
glmm_model <- glmer(species_richness ~ microhabitat + offset(log(sample)) + (1 | site), 
                    data = species_richness, 
                    family = poisson)
summary(glmm_model)
tukey_results <- glht(glmm_model, linfct = mcp(microhabitat = "Tukey"))
summary(tukey_results)

#shannons across entire SB
new1 <- emergence %>%
  group_by(site, site_micro,microhabitat,dispersal_mode,species,growth_form) %>%
  dplyr::summarise(across(c(count),sum))
div <- new1 %>%
  group_by(microhabitat,site) %>%
  summarise_at(vars(count), c("shannon", "invsimpson", "simpson", "richness"))
div
#model
glmm_model <- lmer(shannon ~ microhabitat + (1 | site), 
                   data = div)
# Summary of the model
summary(glmm_model)


#woody
#richness
emergence_combined <- emergence %>%
  mutate(growth_form_combined = ifelse(growth_form %in% c("shrub", "tree"), "woody", growth_form))
richness_woody <- emergence_combined %>%
  filter(growth_form_combined == "woody") %>%
  group_by(microhabitat, site,sample) %>%
  summarise(richness_woody = n_distinct(species)) %>%
  ungroup() %>%
  complete(microhabitat, site, fill = list(richness_woody = 0))
#model
glmm_model <- glmer(richness_woody ~ microhabitat +  offset(log(sample)) +(1 | site), 
                    data = richness_woody, 
                    family = poisson)
summary(glmm_model)
# Perform Tukey's post hoc tests
tukey_results <- glht(glmm_model, linfct = mcp(microhabitat = "Tukey"))
# Summarize the results
summary(tukey_results)

#woody shannons
shannons_diversity <- emergence_combined %>%
  filter(growth_form_combined == "woody") %>%
  group_by(microhabitat, site) %>%
  summarise(shannons_diversity = diversity(count, index = "shannon")) %>%
  spread(site, shannons_diversity, fill = 0)
# Convert from wide to long format
shannons_diversity_long <- shannons_diversity %>%
  pivot_longer(cols = -microhabitat, 
               names_to = "site", 
               values_to = "shannons_diversity")
#model
glmm_model <- lmer(shannons_diversity ~ microhabitat + (1 | site), 
                   data = shannons_diversity_long)
# Summary of the model
summary(glmm_model)# Fit the GLMM
# Perform Tukey's post hoc tests
tukey_results <- glht(glmm_model, linfct = mcp(microhabitat = "Tukey"))
# Summarize the results
summary(tukey_results)


#herbs
#richness
richness_herb <- emergence %>%
  filter(growth_form == "herb") %>%
  group_by(microhabitat, site,sample) %>%
  summarise(richness_herb = n_distinct(species)) %>%
  ungroup() %>%
  complete(microhabitat, site, fill = list(richness_herb = 0))
#model
glmm_model <- glmer(richness_herb ~ microhabitat +  offset(log(sample)) +(1 | site), 
                    data = richness_herb, 
                    family = poisson)
# Summary of the model
summary(glmm_model)# Fit the GLMM
# Perform Tukey's post hoc tests
tukey_results <- glht(glmm_model, linfct = mcp(microhabitat = "Tukey"))
# Summarize the results
summary(tukey_results)

#herbs shannons
shannons_herb <- emergence %>%
  filter(growth_form == "herb") %>%
  group_by(microhabitat, site) %>%
  summarise(shannons_herb = diversity(count, index = "shannon")) %>%
  spread(site, shannons_herb, fill = 0)
# Convert from wide to long format
shannons_diversity_long <- shannons_herb %>%
  pivot_longer(cols = -microhabitat, 
               names_to = "site", 
               values_to = "shannons_diversity")
#model
glmm_model <- lmer(shannons_diversity ~ microhabitat + (1 | site), 
                   data = shannons_diversity_long)
summary(glmm_model)
# Perform Tukey's post hoc tests
tukey_results <- glht(glmm_model, linfct = mcp(microhabitat = "Tukey"))
# Summarize the results
summary(tukey_results)



####FIGURE 3 // NMDS #### 

#NMDS for entire SB community
data <- emergence %>%
  group_by(site, site_micro,microhabitat,species) %>%
  summarise(across(c(count), sum))
data <- data %>%
  pivot_wider(names_from = species, values_from = count)
# Replace NAs with 0s
data[is.na(data)] <- 0
Treatment <-  as.factor(data$microhabitat)
Treatment
data
species <- data [,4:165]
env <- data[,1:3]
ord <- metaMDS(species, distance = "bray", k = 2, try = 20, trymax = 70, autotransform = F, expand = F)
ord
#make  distance matrix
dat<- vegdist(species, "bray")
adonis2(dat ~ microhabitat , data = env,
        strata = env$site)


#woody NMDS
data1 <- read_csv("seed bank/traitordwoody_merged.csv")
Treatment <-  as.factor(data1$microhabitat)
Treatment
data1
species <- data1 [,9:27]
env <- data1[,1:8]
ord <- metaMDS(species, distance = "bray", k = 2, try = 20, trymax = 70, autotransform = F, expand = F)
ord
#make  distance matrix
dat<- vegdist(species, "bray")
adonis2(dat ~ microhabitat , data = env,
        strata = env$site)


#herbs NMDS
herbs <- emergence %>%
  filter(growth_form == "herb") %>%
  group_by(microhabitat, site, site_micro,species) %>%
  summarise(across(c(count), sum))
herbs <- herbs %>%
  pivot_wider(names_from = species, values_from = count)
# Replace NAs with 0s
herbs[is.na(herbs)] <- 0
Treatment <-  as.factor(herbs$microhabitat)
Treatment
herbs
species <- herbs [,4:57]
env <- herbs[,1:3]
ord <- metaMDS(species, distance = "bray", k = 3, try = 20, trymax = 70, autotransform = F, expand = F)
ord
#make  distance matrix
dat<- vegdist(species, "bray")
adonis2(dat ~ microhabitat , data = env,
        strata = env$site)


####FIGURE 4 // GLMMS####

#for all species
species_richness <- emergence %>%
  group_by(micro_tree, site,sample_tree) %>%
  summarise(species_richness = n_distinct(species)) %>%
  ungroup() %>%
  complete(micro_tree, site, fill = list(species_richness = 0))
#model
glmm_model <- glmer(species_richness ~ micro_tree + offset(log(sample_tree))+ (1 | site), 
                    data = species_richness, 
                    family = poisson)
summary(glmm_model)
# Perform Tukey's post hoc tests
tukey_results <- glht(glmm_model, linfct = mcp(micro_tree = "Tukey"))
# Summarize the results
summary(tukey_results)


#bird dispersed species
species_richness <- emergence %>%
  group_by(micro_tree, site, dispersal_mode,sample_tree) %>%
  filter(dispersal_mode == "endozoochory") %>%
  summarise(species_richness = n_distinct(species)) %>%
  ungroup() %>%
  complete(micro_tree, site, dispersal_mode, fill = list(species_richness = 0))
#model
glmm_model <- glmer(species_richness ~ micro_tree  + offset(log(sample_tree)) +(1 | site), 
                    data = species_richness, 
                    family = poisson)
summary(glmm_model)# Fit the GLMM
plot(glmm_model)
# Perform Tukey's post hoc tests
tukey_results <- glht(glmm_model, linfct = mcp(micro_tree = "Tukey"))
# Summarize the results
summary(tukey_results)


#abiotic dispersed seeds
emergence_combined <- emergence %>%
  mutate(disp_combined = ifelse(dispersal_mode %in% c("autochory", "ectozoochory", "abiotic"), "abiotic", dispersal_mode))
richness_abiotic <- emergence_combined %>%
  filter(disp_combined == "abiotic") %>%
  group_by(site, micro_tree,sample_tree) %>%
  summarise(richness_abiotic = n_distinct(species))%>%
  complete(micro_tree, site, fill = list(richness_abiotic = 0))
#model
glmm_model <- glmer(richness_abiotic ~ micro_tree + offset(log(sample_tree)) + (1 | site), 
                    data = richness_abiotic, 
                    family = poisson)
summary(glmm_model)# Fit the GLMM
plot(glmm_model)
# Perform Tukey's post hoc tests
tukey_results <- glht(glmm_model, linfct = mcp(micro_tree = "Tukey"))
# Summarize the results
summary(tukey_results)