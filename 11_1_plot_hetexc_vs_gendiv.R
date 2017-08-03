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
source("martin.R")

# phylogeny
tree_final <- read.tree("data/raw/phylogeny/higdon_mod2_28.tre")
# produce short names for plotting
short <- c("W", "NFS", "SSL", "CSL", "GSL", "SASL", "AFS", "NZSL", "AntFS", "NZFS", "SAFS", "GFS", 
    "BS", "HoS", "GS", "HS", "ARS", "SRS", "BRS", "LRS", "MMS", "HMS", "NES", "SES", "CS", "RS", "LS", "WS")

all_stats <- read_csv("data/processed/all_stats_tree.csv") %>% 
    mutate(SSD = male_weight/female_weight) %>% 
    mutate(abc_out = ifelse(bot > 0.5, "bot", "neut")) %>% 
    mutate(BreedingType = factor(BreedingType, levels = c("ice", "land", "both"))) %>% 
    mutate(logAbundance = log(Abundance),
        logharem_size = log(harem_size),
        logmale_weight = log(male_weight),
        logbreed_season = log(breeding_season_length),
        loglactation_length = log(lactation_length)) %>% 
    # order factors according to tree
    mutate(tip_label = fct_inorder(factor(tip_label)),
        species = fct_inorder(factor(species)),
        latin = fct_inorder(factor(latin)),
        common = fct_inorder(factor(common)),
        short = fct_inorder(factor(short)))
# count grey and harbour seal to land breeding
all_stats[all_stats$BreedingType == "both", "BreedingType"] <- "land" 
all_stats <- all_stats %>% mutate(BreedingType = as.factor(as.character(BreedingType))) %>% data.frame()



## modeling and plotting 

# plot (1) het-exc vs genetic diversity
stats_mod <- all_stats %>% 
    dplyr::select(TPM70_ratio, TPM90_ratio, num_alleles_mean, obs_het_mean, 
        mean_allele_range, prop_low_afs_mean, mratio_mean,
        nloc, nind, bot) %>% 
    data.frame()

p1 <- ggplot(aes(num_alleles_mean, TPM70_ratio), data = stats_mod) +
    geom_smooth(method = "lm", col = "lightgrey", fill = "grey",  alpha = 0.2) +
    geom_point(size = 3, alpha = 0.4) + # abc_out
    geom_point(size = 3, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    xlab("Allelic Richness") +
    ylab("Heterozygosity-excess")
p1

p2 <- ggplot(aes(obs_het_mean, TPM70_ratio), data = stats_mod) +
    geom_smooth(method = "lm", col = "lightgrey", fill = "grey",  alpha = 0.2) +
    geom_point(size = 3, alpha = 0.4) + # abc_out
    geom_point(size = 3, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    xlab("Observed heterozygosity") +
    ylab("Heterozygosity-excess")
p2

p3 <- ggplot(aes(prop_low_afs_mean, TPM70_ratio), data = stats_mod) +
    geom_smooth(method = "lm", col = "lightgrey", fill = "grey",  alpha = 0.2) +
    geom_point(size = 3, alpha = 0.4) + # abc_out
    geom_point(size = 3, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    xlab("Proportion of low frequency alleles") +
    ylab("Heterozygosity-excess") +
    theme(panel.grid.minor = element_blank())
p3

# load model output
mod_out <- read_delim("data/processed/models/mod_gen.txt", delim = " ")
mod_out <- mod_out[-1, 1:4]
names(mod_out) <- c("comps", "pe", "cilow", "cihigh")

p4 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line.x = element_line(color = '#333333'),
          axis.title.y = element_blank()) +
    scale_y_discrete(labels = c("Allelic richness", "Observed\nheterozygosity", 
        "Proportion of\nlow frequence alleles")) +
    xlab("Effect size") +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
    

p_final <- plot_grid(p3, p4, ncol = 2)


ggsave('figures/gen_div_simple.jpg',p_final,  width=8, height=3)



