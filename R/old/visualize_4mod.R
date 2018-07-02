# visualise output from 4 model analysis

library(readr)
library(tidyr)
library(abc)

source("R/martin.R")
mod4_ms <- read_delim(file = "data/processed/sims_1000k_4mods_model_selection_30.txt", delim = " ",
           col_names = c("species", "bot_iceage", "bot", "neut_iceage", "neut"), skip = 1)

mod4_ms_lf <- gather(mod4_ms, model, probability, bot_iceage:neut)

?melt

ggplot(mod4_ms_lf, aes(x = model, y = species, fill = probability)) +
    geom_tile(color = "grey", size = 0.1) +
        scale_fill_viridis(name = "Model probabilities", direction = -1) +
    theme_martin()


#confusion matrix

