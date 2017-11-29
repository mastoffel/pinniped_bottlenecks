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
                  
ndecimal <- c(0, 0,0,0,0,0,2,0,0,0,0)

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

print(xtable(all_stats_table, digits = ndecimal,  align = c("l", "l",">{\\itshape}l", rep("l", 8))),  #hline.after = c(1), 
      include.rownames=FALSE, sanitize.colnames.function=bold, booktabs = TRUE, scalebox = 0.7, floating = TRUE) #floating.environment = "sidewaystable"

# Table 2: Genetic information and outcome of the ABC
transform_ss <- function(PE, CIlow, CIhigh) {
    out <- paste0(PE, " (", CIlow," ,", CIhigh, ")")
    out
}

all_stats_table2 <- all_stats %>% 
                    dplyr::select(common, latin,
                                  num_alleles_mean, num_alleles_mean.CIlow, num_alleles_mean.CIhigh,
                                  obs_het_mean, obs_het_mean.CIlow, obs_het_mean.CIhigh,
                                  exp_het_mean, exp_het_mean.CIlow, exp_het_mean.CIhigh,
                                  prop_low_afs_mean, prop_low_afs_mean.CIlow, prop_low_afs_mean.CIhigh,
                                  mean_allele_range, mean_allele_range.CIlow, mean_allele_range.CIhigh,
                                  mratio_mean, mratio_mean.CIlow, mratio_mean.CIhigh
                                  ) %>% 
                    mutate_if(is.numeric, funs(round(., 2))) 

# bin ci and means together
all_stats_table2_short <- data.frame(do.call(cbind, lapply(seq(from = 3, to = ncol(all_stats_table2 ), by = 3), 
                                function(x) transform_ss(all_stats_table2[[x]], all_stats_table2[[x + 1]], 
                                    all_stats_table2[[x + 2]])))) %>% 
                         dplyr::rename(`Allelic richness` = X1,
                                        `Obs. heterozygosity` = X2,
                                        `Exp. heterozygosity` = X3,
                                        `Prop. low frequency alleles` = X4,
                                        `Allelic range` = X5,
                                        `M-ratio` = X6) %>% 
                          bind_cols(all_stats_table2[1:2], .) %>% 
                             dplyr::arrange(-row_number())

 

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

print(xtable(all_stats_table2_short, align = c("l", "l",">{\\itshape}l", rep("l", 6))),  #hline.after = c(1), 
    include.rownames=FALSE, sanitize.colnames.function=bold, booktabs = TRUE, scalebox = 0.7, 
    floating.environment = "sidewaystable")



