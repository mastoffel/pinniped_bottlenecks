library(waffle)
library(extrafont)

library(readr)
library(readxl)
library(dplyr)
library(ggthemr)
library(viridis)
library(extrafont)
source("Martin.R")
seal_data <- read_xlsx("data/processed/hunting.xlsx")

seal_data <- seal_data %>% mutate(seal_hunt = factor(bottleneck_hist))
levels(seal_data$seal_hunt)
hunts <- c("near extinction", "heavy hunting", 
    "hunted", "unknown",
    "never hunted substantially" )
seal_data$seal_hunt <- factor(seal_data$bottleneck_hist, levels = hunts)


?waffle

waffle_parts <- as.numeric(table(seal_data$seal_hunt))
names(waffle_parts) <- c("near extinction", "heavily hunted", 
    "hunted", "unknown",
    "never hunted substantially" )

waffle(parts = waffle_parts, rows = 3, colors = c("#bd0026", "#f03b20", 
    "#fd8d3c", "grey", "#33a02c"), legend_pos = "right", glyph_size = 10)
