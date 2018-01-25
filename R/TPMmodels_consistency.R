names(all_stats)
library(tidyr)
library(rptR)

all_stats <- read_csv("data/processed/all_stats_tree_29.csv")
# repeatability of het-exc across species
all_stats_long <- all_stats %>% 
    select(species, TPM70_ratio, TPM80_ratio, TPM90_ratio, SMM_ratio) %>% 
    gather(model, value, -species)

rep_mutmod <- rpt(formula = value ~ model + (1|species), data = all_stats_long, datatype = "Gaussian", 
                  grname = "species", nboot = 100, npermut = 100)


all_stats_long <- all_stats %>% 
    select(species, TPM70_ratio, TPM80_ratio, TPM90_ratio, SMM_ratio) %>% 
    gather(model, value, -species) %>% 
    group_by(model) %>% 
    summarise(mean_hetexc = mean(value))

