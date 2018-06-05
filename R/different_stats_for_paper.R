# Different stats for paper
library(readr)
seals <- read_csv("data/processed/all_stats_30_modeling.csv")
names(seals)

# genetic summary
median(seals$nind)
median(seals$nloc)
# correlation
cor(seals$num_alleles_mean, seals$obs_het_mean)
