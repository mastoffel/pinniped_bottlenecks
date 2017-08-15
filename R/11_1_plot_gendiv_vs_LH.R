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

# plot (1) het-exc vs genetic diversity
stats_mod <- all_stats %>% 
    dplyr::select(TPM80_ratio, TPM70_ratio, TPM90_ratio, num_alleles_mean, obs_het_mean, 
        mean_allele_range, prop_low_afs_mean, mratio_mean, logAbundance, SSD, BreedingType,
        nloc, nind, bot, Abundance) %>% 
    data.frame()

# load model output
mod_beta <- read_delim("data/processed/models/mod_div_beta.txt", delim = " ")
mod_R2 <- read_delim("data/processed/models/mod_div_R2.txt", delim = " ")
mod_SC <- read_delim("data/processed/models/mod_div_SC.txt", delim = " ")


names(mod_out) <- c("comps", "pe", "cilow", "cihigh")

p1 <- ggplot(aes(Abundance, num_alleles_mean), data = stats_mod) +
    geom_point(size = 4, alpha = 0.4, aes(color = BreedingType)) + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") +
    geom_line(stat = "smooth", method = "lm",  alpha = 0.6,  aes(color = BreedingType)) +
    geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.08, 
        aes(fill = BreedingType)) +
    scale_color_manual(values = c("cornflowerblue", "black"), name = "Breeding Habitat") +
    scale_fill_manual(values = c("cornflowerblue", "black"), name = "Breeding Habitat") +
    theme_martin() +
    scale_x_log10(breaks = c(100, 1000, 10000, 100000, 1000000), 
        labels = c(expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6})) ) + 
    theme(legend.position=c(0.2, 0.8),
        plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm")) +
    xlab("Abundance") +
    ylab("Allelic richness") 
    #guides(fill = guide_legend(title = "Breeding Habitat"))
p1


# beta coefficients
mod_out <- mod_beta[-1, c("components", "post_median", "lower", "upper")]
names(mod_out) <- c("comps", "pe", "cilow", "cihigh")

p2 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(1,0,0.5,1.5), "cm")) +
    #scale_y_discrete(labels = c("Breeding\nhabitat\nice vs. land", "log(Abundance)", 
    #    "Sexual Size\nDimorphism")) +
    xlab(expression(paste("Effect size ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p2


mod_out_SC <- mod_SC[, c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC ) <- c("comps", "pe", "cilow", "cihigh")

# structure coefficients
p3 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_SC) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = c(0.5)),
        plot.margin = unit(c(1,0.2,0.2,0.1), "cm")) +
    scale_y_discrete(labels = c("Breeding\nhabitat\nice vs. land", "Abundance", 
        "Sexual Size\nDimorphism (SSD)")) +
    xlab(expression(paste("Structure coefficient", " r(", hat(Y),",x)") )) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p3

col_legend <- "#969696"

mod_out_R2 <- mod_R2[, c("combinations", "medianR2", "lower", "upper")]
names(mod_out_R2) <- c("comps", "pe", "cilow", "cihigh")
mod_out_R2$comps <- rev(fct_inorder(factor(mod_out_R2$comps)))
# structure coefficients
p4 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_R2 ) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = c(0.5)),
        plot.margin = unit(c(0.3, 0.2, 0.3, 0.2), "cm")) +
    scale_y_discrete(labels = rev(c("Full model", "Abundance", "Breeding Habitat", "SSD", "Abundance &\nBreeding Habitat",
                    "Abundance &\nSSD", "Breeding Habitat &\nSSD"))) +
    xlab(expression(paste(R^{2}))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
    annotate("segment", x = 0.9, xend = 0.9, y = 0.8, yend = 3.5, color = col_legend) +
    annotate("text", x = 0.98, xend = 0.95, y = 2, yend = 3, color = col_legend, label = c("common"), angle = 270) +
    annotate("segment", x = 0.9, xend = 0.9, y = 4, yend = 6, color = col_legend) +
    annotate("text", x = 0.98, xend = 0.95, y = 5, yend = 6.5, color = col_legend, label = c("unique"), angle = 270) +
    annotate("segment", x = 0.9, xend = 0.9, y = 6.3, yend = 7.5, color = col_legend) +
    annotate("text", x = 0.98, xend = 0.95, y = 7, yend = 8, color = col_legend, label = c("marginal"), angle = 270)
p4


p_top <- plot_grid(p1, p4, rel_widths = c(1.4,1), labels = c("A", "B"),label_fontfamily = "Lato",
    label_x = 0.1, label_y = 0.98)
p_top 
p_bot <- plot_grid(p2, p3, labels = c("C", "D"),label_fontfamily = "Lato",
    label_x = 0.2, label_y = 1, rel_widths = c(1, 1.2))
p_bot

p_final <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(1.4,1))
p_final

ggsave('figures/gen_div_vs_lh.jpg',p_final,  width=6.7, height=5.5)











# 

stats_mod <- all_stats %>% 
    dplyr::select(num_alleles_mean, logbreed_season, SSD,  BreedingType, logAbundance) %>% 
    data.frame()

p1 <- ggplot(aes(logAbundance, num_alleles_mean), data = stats_mod) +
    geom_point(size = 4, alpha = 0.4, aes(color = BreedingType)) + # abc_out size = 4,
    geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") + #size = 4,
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
