## plot het-exc vs ssd for land breeding seals

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

# produce short names for plotting
short <- c("W", "NFS", "SSL", "CSL", "GSL", "AFS","SASL", "NZSL", "AntFS", "NZFS", "SAFS", "GFS", 
    "BS", "HoS", "HS", "GS", "ARS", "BRS", "SRS", "LRS", "MMS", "HMS", "NES", "SES", "RS", "CS", "LS", "WS")

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
    dplyr::select(TPM80_ratio, TPM70_ratio, TPM90_ratio, short, abc_out, bot, SSD, BreedingType, common,
        nloc, nind, harem_size) %>% 
    data.frame()

# load model output
mod_beta <- read_delim("data/processed/models/mod_het_land_beta.txt", delim = " ")
mod_R2 <- read_delim("data/processed/models/mod_het_land_R2.txt", delim = " ")
mod_SC <- read_delim("data/processed/models/mod_het_land_SC.txt", delim = " ")


p1 <- ggplot(aes(SSD, TPM80_ratio, col = BreedingType, fill = BreedingType, alpha = BreedingType), data = stats_mod) +
    geom_line(stat = "smooth", method = "lm") + #, aes(fill = BreedingType)
    geom_ribbon(stat='smooth', method = "lm", alpha = 0.1, colour=NA) +
    geom_point(size = 4) + # abc_out
    geom_point(size = 4, shape = 21, col = "black") + #, aes(fill = BreedingType)
    #scale_color_viridis(option = "magma", direction = -1,
    #    name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    theme_martin() +
    xlab("SSD") +
    ylab("Heterozygosity-excess") +
    scale_x_log10(breaks = c(1,2,3,4,5,6,7,8)) +
    #ylab("Heterozygosity-excess") +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.2,0.2), "cm") ,
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
        legend.direction = "vertical",
        legend.position = c(0.2,0.85),
        legend.title=element_text(size=10),
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))
    ) +
    scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    scale_alpha_manual(values = c("0.1", "0.8")) +
    guides(color = guide_legend("Breeding habitat"),
            fill = guide_legend("Breeding habitat"),
            alpha = guide_legend("Breeding habitat")) + #, label.position = "bottom"
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) +
    geom_text_repel(label = short,size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
        segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
        segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines"))
p1

### old plotting
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
        plot.margin = unit(c(0,0.5,0,0.45), "cm")) +
    xlim(c(-0.03, 0.3)) +
    xlab(expression(paste("Effect   size ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p2

col_legend <- "#969696"
mod_out_R2 <- mod_R2[, c("combinations", "medianR2", "lower", "upper")][-4, ]
names(mod_out_R2) <- c("comps", "pe", "cilow", "cihigh")
mod_out_R2$comps <- factor(mod_out_R2$comps, levels = c("BreedingType", "SSD", "full model"))
mod_out_R2 <- mod_out_R2[-2, ]
# mod_out_R2$comps <- rev(fct_inorder(factor(mod_out_R2$comps)))
# structure coefficients
p3 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_R2 ) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3.5, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = c(0.5)),
        plot.margin = unit(c(0,0.5,0,0), "cm")) +
    xlim(-0.1,1) +
    scale_y_discrete(labels = c("Full\nmodel")) +
    xlab(expression(paste(R^{2}))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1) 
p3

forestplots <- plot_grid(p3, p2, ncol = 1, labels = c("B", "C"), label_fontfamily = "Lato")
forestplots


p_final <- plot_grid(p1, forestplots, rel_widths = c(2.3,1), labels = c("A", ""), label_fontfamily = "Lato")

ggplot2::ggsave('figures/het_vs_SSD_landbreed.jpg',p_final,  width=7, height=4)
# Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
# embed_fonts("figures/het_vs_SSD_landbreed.pdf")

