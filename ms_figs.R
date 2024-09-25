
#Title: Figures for MS submission
#Author: Sam Allbee

library(readr)
library(dplyr)
library(lme4)
library(multcomp)
library(readr)
library(factoextra)
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(tidyverse)
library(vegan)
library(goeveg)
library(plot3D)
library(vegan3d)
library(cowplot)
library(gridExtra)
library(devtools)
library(abdiv)
library(scales)

col <- c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')

col <- c('#F6CF71', '#ce6693', '#87C55F')


col <- c('#F6CF71', '#bf812d', '#87C55F')

col <- c('#bf812d', '#dfc27d', '#01665e')

####FIG 1####

emergence <- read_csv("seed bank/cleaned_final_MS.csv")

#fig 1 total richness
species_richness <- emergence %>%
  group_by(microhabitat, site,sample) %>%
  summarise(species_richness = n_distinct(species)) %>%
  ungroup()
rich <- ggplot(species_richness, aes(x = reorder(microhabitat, species_richness), y = species_richness/sample, fill = microhabitat)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  ylab("Spp richness") +
  xlab("Habitat") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") +
  scale_x_discrete(labels = c("open pasture" = "OP", "isolated tree" = "IT", "forest fragment" = "FF")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

#fig 1 total shannons
new1 <- emergence %>%
  group_by(site, site_micro,microhabitat,dispersal_mode,species,growth_form) %>%
  dplyr::summarise(across(c(count),sum))
div <- new1 %>%
  group_by(microhabitat,site) %>%
  summarise_at(vars(count), c("shannon", "invsimpson", "simpson", "richness"))
div
#plot
shan <-ggplot(div, aes(x = reorder(microhabitat, shannon), y = shannon, fill = microhabitat)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  ylab("Diversity (H')") +
  xlab("Habitat") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"), # Change size and color as needed
        axis.title.y = element_text(size = 14, color = "black"), # Change size and color as needed
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") +
  scale_x_discrete(labels = c("open pasture" = "OP", "isolated tree" = "IT", "forest fragment" = "FF")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

####WOODY FIG 1####
#fig 1 woody richness
emergence_combined <- emergence %>%
  mutate(growth_form_combined = ifelse(growth_form %in% c("shrub", "tree"), "woody", growth_form))
richness_woody <- emergence_combined %>%
  filter(growth_form_combined == "woody") %>%
  group_by(microhabitat, site,sample) %>%
  summarise(richness_woody = n_distinct(species)) %>%
  ungroup() %>%
  complete(microhabitat, site, fill = list(richness_woody = 0))
#boxplot
wood_rich <- ggplot(richness_woody, aes(x = reorder(microhabitat, richness_woody), y = richness_woody/sample, fill = microhabitat)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  ylab("Woody spp richness") +
  xlab("Habitat") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none")+
  scale_x_discrete(labels = c("open pasture" = "OP", "isolated tree" = "IT", "forest fragment" = "FF")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

#fig 1 woody shannon

# Calculate Shannon's diversity index for woody stems across each microhabitat and the ten sites
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

# Print the resulting long-format dataset
print(shannons_diversity_long)

shan_wood <- ggplot(shannons_diversity_long, aes(x = reorder(microhabitat, shannons_diversity), y = shannons_diversity, fill = microhabitat)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  ylab("Woody diversity (H')") +
  xlab("Habitat") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") +
  scale_x_discrete(labels = c("open pasture" = "OP", "isolated tree" = "IT", "forest fragment" = "FF")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))


####HERBS####

richness_herb <- emergence %>%
  filter(growth_form == "herb") %>%
  group_by(microhabitat, site,sample) %>%
  summarise(richness_herb = n_distinct(species)) %>%
  ungroup() %>%
  complete(microhabitat, site, fill = list(richness_herb = 0))

#boxplot
herb <- ggplot(richness_herb, aes(x = reorder(microhabitat, richness_herb), y = richness_herb/sample, fill = microhabitat)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  ylab("Herb spp richness") +
  xlab("Habitat") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") +
  scale_x_discrete(labels = c("open pasture" = "OP", "isolated tree" = "IT", "forest fragment" = "FF")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))

#### herbs shannons #####

shannons_herb <- emergence %>%
  filter(growth_form == "herb") %>%
  group_by(microhabitat, site) %>%
  summarise(shannons_herb = diversity(count, index = "shannon")) %>%
  spread(site, shannons_herb, fill = 0)

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

library(vegan)

# Print the resulting long-format dataset
print(shannons_diversity_long)

herb_sh <- ggplot(shannons_diversity_long, aes(x = reorder(microhabitat, shannons_diversity), y = shannons_diversity, fill = microhabitat)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  ylab("Herb diversity (H')") +
  xlab("Habitat") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") +
  scale_x_discrete(labels = c("open pasture" = "OP", "isolated tree" = "IT", "forest fragment" = "FF")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
#####FIG 1 ALL PANEL#####
#this is the one
combined_plot <- grid.arrange(rich, shan, wood_rich, shan_wood, herb, herb_sh, ncol = 2)

####FIG THREE NMDS####

set.seed(1)

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
ord <- metaMDS(species, distance = "bray", k = 2, try = 20, trymax = 70, autotransform = F, expand = F)

ordination_scores <- ord$points
scrs <- cbind(as.data.frame(ordination_scores), Treatment = Treatment)

#NMDS full dataset
entire <- ggplot(scrs,
       aes(x = MDS1,
           y = MDS2,
           color = Treatment)) + 
  geom_point() + 
  stat_ellipse(size = 1.2,level = 0.85) +  # Adjust the level parameter to make the ellipses closer to the points
  theme_bw() +
  xlab("NMDS1") + 
  ylab("NMDS2") +
  scale_color_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  theme_bw() + #add borders
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 16,
                                  color = "black"),#face = "bold" 
        axis.text = element_text(size = 16,
                                 color = "black"),
        axis.text.y = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.position = "none", 
        aspect.ratio = 1)


#woody nmds
data1 <- read_csv("seed bank/traitordwoody_merged.csv")
Treatment <-  as.factor(data1$microhabitat)
Treatment
data1
species <- data1 [,9:27]
ord <- metaMDS(species, distance = "bray", k = 2, try = 20, trymax = 70, autotransform = F, expand = F)
ord
ordination_scores <- ord$points
scrs <- cbind(as.data.frame(ordination_scores), Treatment = Treatment)

library(scales)
woody <- ggplot(scrs,
                 aes(x = MDS1,
                     y = MDS2,
                     color = Treatment)) + 
  geom_point() + 
  stat_ellipse(size = 1.2,level = 0.80) +  # Adjust the level parameter to make the ellipses closer to the points
  theme_bw() +
  xlab("NMDS1") + 
  ylab("NMDS2") +
  scale_color_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  theme_bw() + #add borders
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 16,
                                  color = "black"),#face = "bold" 
        axis.text = element_text(size = 16,
                                 color = "black"),
        axis.text.y = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.position = "none", 
        aspect.ratio = 1)

ggplot(scrs,
       aes(x = MDS1,
           y = MDS2,
           color = Treatment)) + 
  geom_point() + 
  stat_ellipse(size = 1.2, level = 0.80) +  # Adjust the level parameter to make the ellipses closer to the points
  theme_bw() +
  xlab("NMDS1") + 
  ylab("NMDS2") +
  scale_color_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  theme_bw() + # Add borders
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 16, color = "black"), # face = "bold" 
        axis.text = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.position = "none", 
        aspect.ratio = 1) +
  scale_x_continuous(labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(labels = label_number(accuracy = 0.1))

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
ord <- metaMDS(species, distance = "bray", k = 2, try = 20, trymax = 70, autotransform = F, expand = F)
ord
ordination_scores <- ord$points
scrs <- cbind(as.data.frame(ordination_scores), Treatment = Treatment)

herb <- ggplot(scrs,
       aes(x = MDS1,
           y = MDS2,
           color = Treatment)) + 
  geom_point() + 
  stat_ellipse(size = 1.2,level = 0.8) +  # Adjust the level parameter to make the ellipses closer to the points
  theme_bw() +
  xlab("NMDS1") + 
  ylab("NMDS2") +
  scale_color_manual(values = c("open pasture" = "#F6CF71", "isolated tree" = "#ce6693", "forest fragment" = "#87C55F")) +
  theme_bw() + #add borders
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 16,
                                  color = "black"),#face = "bold" 
        axis.text = element_text(size = 16,
                                 color = "black"),
        axis.text.y = element_text(size = 16),
        axis.ticks = element_blank(),
        legend.position = "none", 
        aspect.ratio = 1) 


####FIG 4####

#everything 4 micros

desired_order <- c("open pasture", "non-fleshy fruit tree", "fleshy fruit tree", "forest fragment")

# Convert the "micro_tree" variable to a factor with the desired order
species_richness$micro_tree <- factor(species_richness$micro_tree, levels = desired_order)

species_richness <- emergence %>%
  group_by(micro_tree, site,sample_tree) %>%
  summarise(species_richness = n_distinct(species)) %>%
  ungroup() %>%
  complete(micro_tree, site, fill = list(species_richness = 0))

all <- ggplot(species_richness, aes(x = micro_tree, y = species_richness/sample_tree, fill = micro_tree)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "fleshy fruit tree" = "#ce6693", "non-fleshy fruit tree" = "#E69F00", "forest fragment" = "#87C55F")) +
  ylab("Spp richness") +
  xlab("Habitat") +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none") +
  scale_x_discrete(labels = c("open pasture" = "OP", "non-fleshy fruit tree" = "NFT", "fleshy fruit tree" = "FFT", "forest fragment" = "FF"))

#this code is to get abiotic dispersed
emergence_combined <- emergence %>%
  mutate(disp_combined = ifelse(dispersal_mode %in% c("autochory", "ectozoochory", "abiotic"), "abiotic", dispersal_mode))

richness_abiotic <- emergence_combined %>%
  filter(disp_combined == "abiotic") %>%
  group_by(site,micro_tree,sample_tree) %>%
  summarise(richness_abiotic = n_distinct(species))

desired_order <- c("open pasture", "non-fleshy fruit tree", "fleshy fruit tree", "forest fragment")

# Convert the "micro_tree" variable to a factor with the desired order
richness_abiotic$micro_tree <- factor(richness_abiotic$micro_tree, levels = desired_order)

# Plot the boxplots with the specified order
abiotic_fleshy <- ggplot(richness_abiotic, aes(x = micro_tree, y = richness_abiotic/sample_tree, fill = micro_tree)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "fleshy fruit tree" = "#ce6693", "non-fleshy fruit tree" = "#E69F00", "forest fragment" = "#87C55F")) +
  ylab("AD spp richness") +
  xlab("Habitat") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none")+
  scale_x_discrete(labels = c("open pasture" = "OP", "non-fleshy fruit tree" = "NFT", "fleshy fruit tree" = "FFT", "forest fragment" = "FF"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))


#bird dispersed

species_richness <- emergence %>%
  group_by(micro_tree, site, dispersal_mode,sample_tree) %>%
  filter(dispersal_mode == "endozoochory") %>%
  summarise(species_richness = n_distinct(species)) %>%
  ungroup() %>%
  complete(micro_tree, site, dispersal_mode, fill = list(species_richness = 0))

species_richness <- read_csv("seed bank/species_richness.csv")

bird_disp <- ggplot(species_richness, aes(x = reorder(micro_tree, species_richness, FUN = median), y = species_richness/sample_tree, fill = micro_tree)) +
  geom_boxplot() + 
  scale_fill_manual(values = c("open pasture" = "#F6CF71", "fleshy fruit tree" = "#ce6693", "non-fleshy fruit tree" = "#E69F00", "forest fragment" = "#87C55F")) +
  ylab("BD spp richness") +
  xlab("Habitat") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 0, size = 14, color = "black"),
        axis.title.x = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.1, "cm"),
        legend.position = "none")+
  scale_x_discrete(labels = c("open pasture" = "OP", "non-fleshy fruit tree" = "NFT", "fleshy fruit tree" = "FFT", "forest fragment" = "FF"))+
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1))


bird_disp

desired_order <- c("open pasture", "non-fleshy fruit tree", "fleshy fruit tree", "forest fragment")

# Convert the "micro_tree" variable to a factor with the desired order
species_richness$micro_tree <- factor(species_richness$micro_tree, levels = desired_order)

grid.arrange(all, arrangeGrob(bird_disp, abiotic_fleshy), ncol = 1)
