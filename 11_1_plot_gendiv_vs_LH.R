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
library(extrafont)
library(extrafontdb)
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


mod_out <- read_delim("data/processed/models/mod_gen_vs_lh.txt", delim = " ")
mod_out <- mod_out[-1, 1:4]
names(mod_out) <- c("comps", "pe", "cilow", "cihigh")

# plot

stats_mod <- all_stats %>% 
    dplyr::select(num_alleles_mean, logbreed_season, SSD,  BreedingType, logAbundance) %>% 
    data.frame()

p1 <- ggplot(aes(logAbundance, num_alleles_mean), data = stats_mod) +
    geom_point(size = 4, alpha = 0.4, aes(color = BreedingType)) + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") +
    geom_line(stat = "smooth", method = "lm",  alpha = 0.3,  aes(color = BreedingType)) +
    geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.1, 
        aes(fill = BreedingType)) +
    scale_color_manual(values = c("black", "cornflowerblue")) +
    scale_fill_manual(values = c("black", "cornflowerblue")) +
    theme_martin() +
    theme(legend.position=c(0.2, 0.8)) +
    xlab("log(Abundance)") +
    ylab("Allelic richness")
p1

mod_out
p2 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank()) +
    scale_y_discrete(labels = c("Ice-breeding", "Generation time", 
        "log(Global\nabundance)", "Sexual\nweight dimorphism")) +
    xlab("Effect size") +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p2


p_final <- plot_grid(p1, p2, ncol = 2)

ggsave('figures/gen_div_vs_lh.jpg',p_final,  width=8, height=3)
