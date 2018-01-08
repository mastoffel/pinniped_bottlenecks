## This script produces the figure for the models of bottleneck signatures
# explained by life-history traits.

# phylogenetic comparative analysis
library(ggtree)
library(ape)
library(phytools)
library(dplyr)
library(readxl)
library(stringr)
library(viridis)
library(ggtree)
library(ggthemr)
library(reshape2)
library(cowplot)
library(ggthemes)
library(ggimage)
library(RColorBrewer)
library(scales)
library(forcats)
library(readr)
# for comparative analysis
library(caper)
library(yhat)
library(dplyr)
library(ggrepel)
library(GGally)
library(ggthemr)
source("R/martin.R")
library(MCMCglmm)
library(purrr)
library(readr)
library(xtable)
library(stargazer)
## what should this script do:


# load data and prepare mixed models

# load (modified) phylogeney. 26 species from 10ktrees plus 3 subspecies of ringed seal
tree_final <- read.tree("data/raw/phylogeny/29_species_10ktrees.tre")

# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_29_modeling.csv"))

# number of genotypes
sum(all_stats$nind * all_stats$nloc)



# Table 1: IUCN, life history data and sample size/locus information
all_stats_table <- all_stats %>% 
    dplyr::select(common, latin, IUCN_rating, Abundance, BreedingType, SSD, harem_size, nloc, nind) %>% 
    dplyr::mutate(Genotypes = nloc * nind) %>% 
    dplyr::rename(`Common name` = common,
        `Scientific` = latin,
        `IUCN status` = IUCN_rating,
        `Breeding habitat` = BreedingType,
        `Harem size` = harem_size,
        `Loci` = nloc,
        `Individuals` = nind) %>% 
    dplyr::arrange(-row_number())

stargazer(all_stats_table, summary = FALSE)