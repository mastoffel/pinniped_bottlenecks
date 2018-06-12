# prepare data frame for modeling

# prepare the data.frame for modeling

library(readr)
library(ape)
library(dplyr)
library(forcats)
# load (modified) phylogeney. 26 species from 10ktrees plus 3 subspecies of ringed seal
tree_final <- read.tree("data/raw/phylogeny/30_species_10ktrees_final.tre")

# all_stats_tree is from 10_visualise_phylogeny.R
all_stats <- read_csv("data/processed/all_stats_tree_30.csv") %>% 
    mutate(SSD = male_weight/female_weight) %>% 
    mutate(abc_out = ifelse(bot > 0.5, "bot", "neut")) %>% 
    mutate(BreedingType = factor(BreedingType, levels = c("ice", "land", "both"))) %>% 
    mutate(logAbundance = log(Abundance),
        logharem_size = log(harem_size),
        logmale_weight = log(male_weight),
        logbreed_season = log(breeding_season_length),
        loglactation_length = log(lactation_length),
        logSSD = log(SSD)) %>% 
    # order factors according to tree
    mutate(tip_label = fct_inorder(factor(tip_label)),
        species = fct_inorder(factor(species)),
        latin = fct_inorder(factor(latin)),
        common = fct_inorder(factor(common)),
        short = fct_inorder(factor(short)))
# count grey and harbour seal to land breeding
all_stats[all_stats$BreedingType == "both", "BreedingType"] <- "land" 
all_stats <- all_stats %>% mutate(BreedingType = as.factor(as.character(BreedingType))) %>% data.frame()

write_csv(all_stats, "data/processed/all_stats_30_modeling.csv")


# add on
library(dplyr)
library(tibble)
library(tidyr)
# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_30_modeling.csv"))

all_stats_sub <- all_stats %>% 
    dplyr::select(species, lactation_length, breeding_season_length, Generation_time, life_span_years)

all_stats <- all_stats %>% 
                replace_na(list(lactation_length = 56, life_span_years = 46, breeding_season_length = 60)) %>% 
                mutate(logbreed_season = log(breeding_season_length)) 
# Breeding season length mediterranean monk seal (mediterranean colonies)
all_stats[all_stats$species == "mediterranean_monk_seal", "breeding_season_length"] <- 60

# wtrite to file
write_csv(all_stats, "data/processed/all_stats_30_modeling.csv")


