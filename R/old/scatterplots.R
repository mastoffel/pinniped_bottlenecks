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
# load phyologenetic sequence

# load all datasets
seals <- read_excel("data/processed/seal_data_complete.xlsx")
names(seals)[11] <- "IUCN_rating"

# check which variables are non-numeric
non_numeric <- apply(seals, 2, is.numeric)

seal_names_real <- seals$seal_names_real
str(seals)

arrow <- arrow(length = unit(0.4, "cm"), type = "closed")


# ssd vs harem size
ggplot(seals, aes(mean_SSD, sqrt(harem_size))) +
    geom_point(size = 3, alpha = 0.5) +
    theme_minimal() + 
    theme(
        axis.line = element_line(arrow = arrow)
    )

summary(lm(mean_SSD ~ harem_size, data = seals))


summary(glm(TPM70_ratio ~ mean_SSD + nloc + nind, data = seals))

