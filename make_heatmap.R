library(tidyverse)
library(ape)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)


tree <- read.tree("C:/Users/Nnamdi/Desktop/thesisMain/data/rawData/lgc_protein_rooted.tree")  # phylogenetic tree of lactobacilli species
treeTips = as.vector(tree$tip.label)

clusters = read.csv("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/exploratoryData/groups.csv")  # dataset generated with geneCluster.py
family = read.csv("C:/Users/Nnamdi/Desktop/thesisMain/data/generatedData/exploratoryData/taxofam.csv")

## heatmap of prophage clusters
head(clusters)
clusters %>%
  rename(accession = genomeID, lactobacilli_species = lactobacilli) %>%
  mutate(accession = factor(accession, levels = tree$tip.label)) %>%
  arrange(accession) %>%
  mutate(lactobacilli_species = factor(lactobacilli_species, levels = lactobacilli_species)) %>%
  gather(
    key = "phage_clusters", value = "copy_number", starts_with("clusters")
  ) %>%
  ggplot(aes(x = phage_clusters, y = lactobacilli_species, fill = copy_number)) +
  geom_tile(color = "black") +
  geom_rug(aes(y = lactobacilli_species, col = lifestyle),size = 1, sides = "lr",outside=T) +
  coord_cartesian(clip = "off") +
  scale_color_distiller(type = "div" , palette = "Accent",aesthetics ="fill",guide="legend") +
  theme_dark() +
  theme(
    axis.text.y = element_text(vjust = 0.7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("C:/Users/Nnamdi/Desktop/thesisMain/results/heatmapCluster.png", units = "cm", height = 40, width = 20)


## heatmap of phage family distribution

cleaned = family #remove the non-essential columns
# select columns to serve as ID variavle and columns to be merged into rows (count values will be obtained from these rows)
meltCleaned <- melt(cleaned,var1='lactobacilli',var2=c('DNA','Herelleviridae','Myoviridae','Siphoviridae','Podoviridae'))
colnames(meltCleaned) <- c('lactobacilli_species','phageFamily','value') # rename columns 
ggplot(meltCleaned, aes(phageFamily,lactobacilli_species)) +
  geom_tile(aes(fill = value), colour = "black") +
  coord_cartesian(clip = "off") +
  scale_color_distiller(type = "div" , palette = "Accent",aesthetics ="fill",guide="legend") +
  theme_dark() +
  theme(
    axis.text.y = element_text(vjust = 0.7),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave("C:/Users/Nnamdi/Desktop/thesisMain/results/heatmapFam.png", units = "cm", height = 40, width = 20)



