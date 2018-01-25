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
all_stats_origin <- read_xlsx("data/raw/table_data.xlsx")
all_stats_table <- all_stats %>% 
                    left_join(all_stats_origin, by = c("species", "common", "tip_label", "latin")) %>% 
                    dplyr::mutate(Genotypes = nloc * nind) %>% 
                    dplyr::select(common, latin, IUCN_rating, Abundance, BreedingType, 
                                  SSD, harem_size, nloc, nind, Genotypes, origin, published_short) %>% 
                    dplyr::rename(`Common name` = common,
                           `Scientific` = latin,
                           `IUCN status` = IUCN_rating,
                           `Breeding habitat` = BreedingType,
                           `Harem size` = harem_size,
                           `Loci` = nloc,
                           `Individuals` = nind,
                           `Sampling location` = origin,
                           `Published` = published_short) %>% 
                    dplyr::arrange(-row_number())
library(magrittr)
library(huxtable)
ht <- as_hux(all_stats_table[1:5], add_colnames = TRUE) %>% 
    set_bold(1, everywhere, TRUE)       %>%
    set_bottom_border(1, everywhere, 1) %>%
  #  set_align(everywhere, 2, 'right')   %>%
    set_right_padding(10)               %>%
    set_left_padding(10)                %>%
    set_width(0.9)                     %>%
    set_number_format(2)
quick_pdf(ht, file = "hux_out.pdf")

ht
right_padding(ht) <- 10
left_padding(ht) <- 10
bold(ht)[1,] <- TRUE
bottom_border(ht)[1,] <- 1
number_format(ht) <- 3
width(ht) <- 0.1
wrap(ht) <- TRUE

quick_docx(ht, file = "hux_out.docx")


ndecimal <- c(0, 0,0,0,0,0,2,0,0,0,0,0,0)

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}

print(xtable(all_stats_table, digits = ndecimal,  align = c("l", "l",">{\\itshape}l", rep("l", 10))),  #hline.after = c(1), 
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



