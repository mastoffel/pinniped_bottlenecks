# create a table for the paper

# phylogenetic comparative analysis
library(ggtree)
library(ape)
library(phytools)
library(dplyr)
library(readxl)
library(stringr)
library(ggthemr)
library(reshape2)
library(scales)
library(forcats)
library(readr)
# for comparative analysis
library(caper)
library(yhat)
library(dplyr)
library(GGally)
library(reshape2)
library(magrittr)
library(tibble)
library(mcmcR2)
library(ggrepel)
#source("martin.R")
library(xtable)

#load modified phylogeny and all stats
tree_final <- read.tree("data/raw/phylogeny/28_species_10ktrees.tre")
#load modified phylogeny and all stats
#tree_final <- read.tree("data/raw/phylogeny/higdon_mod2_28.tre")
# produce short names for plotting
short <- c("W", "NFS", "SSL", "CSL", "GSL", "SASL", "AFS", "NZSL", "AntFS", "NZFS", "SAFS", "GFS", 
    "BS", "HoS", "GS", "HS", "ARS", "SRS", "BRS", "LRS", "MMS", "HMS", "NES", "SES", "CS", "RS", "LS", "WS")

# all_stats_tree is from 10_visualise_phylogeny.R
all_stats <- read_csv("data/processed/all_stats_tree.csv") %>% 
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





# Table 1: IUCN, life history data and sample size/locus information
all_stats_table <- all_stats %>% 
                    dplyr::select(common, latin, IUCN_rating, Abundance, BreedingType, SSD, harem_size, nloc, nind) %>% 
                    dplyr::rename(`Common name` = common,
                           `Scientific` = latin,
                           `IUCN status` = IUCN_rating,
                           `Breeding habitat` = BreedingType,
                           `Harem size` = harem_size,
                           `Loci` = nloc,
                           `Sample size` = nind) %>% 
                    dplyr::arrange(-row_number())
                  
ndecimal <- c(0, 0,0,0,0,0,2,0,0,0)

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

print(xtable(all_stats_table, digits = ndecimal,  align = c("l", "l",">{\\itshape}l", rep("l", 7))),  #hline.after = c(1), 
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



