# commonality analysis

library(sealABC)
library(readr)
library(readxl)
library(dplyr)
# data prep ----------------------------------------------------------------------------------------

# load descriptive data
load("data/processed/seal_ss_rarefaction16ind.RData")
# load all datasets
seals <- read_excel("data/processed/seal_data_complete.xlsx")
names(seals)[11] <- "IUCN_rating"
# load model probablities from ABC
model_probs <- read.table("data/processed/sims_5000k_large_bot_model_selection.txt")
# are all names overlapping? --> if 28 then yes..
sum(rownames(all_sumstats_full) %in% seals$species)
# reorder all_sumstats_full (output from rarefaction summary stats calculation with sealABC::mssumstats) 
# and seals (full descriptive dataset)
# reorder LIKE seals$species
ind_reorder <- NULL
for (i in seals$species){
    ind_reorder <- c(ind_reorder, which(rownames(all_sumstats_full) == i))
}
# reorder sumstats according to seal descriptives file
sumstats <- all_sumstats_full[ind_reorder, ]
modelprobs <- model_probs[ind_reorder, ] # works because modelprobs has same sequence as sumstats
# quick check
sum(rownames(sumstats) == rownames(modelprobs))
sum(rownames(sumstats) == seals$species)
# create new data.frame as a combination

all_stats <- cbind(seals[c(1:17, 36:56)], sumstats, modelprobs) %>% 
    mutate(abc = ifelse(bot > 0.5, "bot", "neut"))
all_stats_rownames <- rownames(all_stats)

# select important variables
all_stats <- as_tibble(all_stats)
all_stats


