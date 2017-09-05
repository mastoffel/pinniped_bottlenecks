library(readr)
library(readxl)
library(dplyr)
library(ggthemr)
library(viridis)
library(extrafont)
source("Martin.R")
seal_data <- read_xlsx("data/processed/all_data_seals_handbook.xlsx", skip = 1)

seal_data <- seal_data %>% mutate(seal_hunt = factor(bottleneck_hist))
levels(seal_data$seal_hunt)
hunts <- c("near extinction", "heavy hunting", 
    "hunted", "other severe anthropogenic influences","unknown",
    "never hunted substantially" )
seal_data$seal_hunt <- factor(seal_data$seal_hunt, levels = hunts)



# ggthemr(palette = "fresh", layout = "clear", spacing = 4, text_size = 15)
p <- ggplot(data = seal_data, aes(x = seal_hunt, fill = seal_hunt)) +
    geom_dotplot(stackdir = "up", binwidth = 0.3, width = 2) +
    scale_x_discrete(labels = c("near \nextinction", "heavy \nhunting", 
        "hunted", "other severe \nanthropogenic \ninfluences","unknown",
        "never hunted \nsubstantially" )) +
    labs(x="reported hunting",y="species count") +
    scale_y_continuous("Number of Species", breaks = NULL, limits = c(0,0.5)) +
    theme_martin() +
    scale_fill_brewer(palette = "YlGnBu", type = "seq", direction = -1) +
    theme(
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0)),
        legend.position="none" ,
        panel.grid.major = element_blank()
    )
  
p
ggsave(p, filename = "figures/hunting.jpg", width = 6, height = 4)

