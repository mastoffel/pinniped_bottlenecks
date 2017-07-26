library(readr)
library(readxl)
library(dplyr)
library(ggthemr)
library(viridis)
seal_data <- read_xlsx("data/processed/all_data_seals_handbook.xlsx", skip = 1)

seal_data <- seal_data %>% mutate(seal_hunt = factor(bottleneck_hist))
levels(seal_data$seal_hunt)
hunts <- c("near extinction", "heavy hunting", 
    "hunted", "other severe anthropogenic influences","unknown",
    "never hunted substantially" )
seal_data$seal_hunt <- factor(seal_data$seal_hunt, levels = hunts)

ggthemr(palette = "fresh", layout = "clear", spacing = 4, text_size = 15)
ggplot(data = seal_data, aes(x = seal_hunt)) +
    geom_bar() +
    scale_x_discrete(labels = c("near \nextinction", "heavy \nhunting", 
        "hunted", "other severe \nanthropogenic \ninfluences","unknown",
        "never hunted \nsubstantially" )) +
    labs(x="reported hunting",y="species count")

