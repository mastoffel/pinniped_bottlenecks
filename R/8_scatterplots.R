# plots for the talk
# install_github("mastoffel/sealABC", dependecies = TRUE)
library(stringr)
library(readxl)
library(reshape2)
library(dplyr)
library(ggthemes)
library(tidyr)
library(stringr)
source("R/multiplot.R")
library(hierfstat)
library(ggthemes)
library(gridExtra)
library(viridis)
library(ggplot2)
library(ggthemr)
library(cowplot)
library(sealABC)
library(ggrepel)
# load phyologenetic sequence

# load all datasets
seals <- read_excel("data/processed/seal_data_complete.xlsx")
names(seals)[11] <- "IUCN_rating"
seals$weight_ratio <- seals$male_weight / seals$female_weight

# check which variables are non-numeric
non_numeric <- apply(seals, 2, is.numeric)
seal_names_real <- seals$seal_names_real

# set theme
ggthemr(palette = "fresh", layout = "scientific", spacing = 2.5, line_weight = 0.7, text_size = 18,
    type = "outer")

ggplot(seals, aes(log(harem_size), TPM90_ratio, label = seal_names_real)) +
    geom_smooth(method = "lm", alpha = 0.1, size = 0.5, colour = "grey") + #, se = FALSE
    geom_point(aes( color = BreedingType),size = 3, alpha = 0.8) +
    geom_text_repel(size = 3) +
    theme(legend.position = "top") +
    xlab("log(harem size)") +
    ylab("genetic bottleneck strength") +
    legend_right()


library(lme4)
mod <- lm(TPM95_ratio~weight_ratio + nloc + nind, data = seals)
summary(mod)

arrow <- arrow(length = unit(0.4, "cm"), type = "closed")


# ssd vs harem size
ggplot(seals, aes(mean_SSD, sqrt(harem_size))) +
    geom_point(size = 3, alpha = 0.5) +
    theme_minimal() 

summary(lm(mean_SSD ~ harem_size, data = seals))
summary(glm(TPM70_ratio ~ mean_SSD + nloc + nind, data = seals))


# several summary statistics
library(ggvis)
seals %>%  
ggvis(x = ~species, y = input_select(names(seals)[10:length(names(seals))], map = as.name)) %>% 
    layer_points(fill =~ BreedingType) %>% 
    add_axis("x", properties = axis_props(labels = list(angle = 35, 
        align = "left", baseline = "middle")), title_offset = 100) 
    #add_legend("size", properties = legend_props(legend = list(y = 200)))


# 

seals %>%  
    ggvis(x = input_select(names(seals)[10:length(names(seals))], map = as.name), y =~ TPM70_ratio) %>% 
    add_axis("x", properties = axis_props(labels = list(angle = 35, 
    align = "left", baseline = "middle")))  %>% 
    layer_points(fill =~ BreedingType)

