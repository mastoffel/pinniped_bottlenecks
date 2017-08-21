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

# produce short names for plotting
short <- c("W", "NFS", "SSL", "CSL", "GSL", "AFS","SASL", "NZSL", "AntFS", "NZFS", "SAFS", "GFS", 
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
    dplyr::select(TPM80_ratio, TPM70_ratio, TPM90_ratio, short, abc_out, bot, SSD, BreedingType, common,
        nloc, nind, harem_size) %>% 
    data.frame()

# load model output
mod_beta <- read_delim("data/processed/models/mod_het_beta.txt", delim = " ")
mod_R2 <- read_delim("data/processed/models/mod_het_R2.txt", delim = " ")
mod_SC <- read_delim("data/processed/models/mod_het_SC.txt", delim = " ")

p1 <- ggplot(aes(harem_size, TPM80_ratio, col = bot), data = stats_mod) +
    geom_smooth(method = "lm", col = "lightgrey",   alpha = 0.1) + #, aes(fill = BreedingType)
    geom_point(size = 3.5, alpha = 0.6) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, col = "black") + #, aes(fill = BreedingType)
    #scale_color_viridis(option = "magma", direction = -1,
    #    name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    theme_martin() +
    xlab("SSD") +
    ylab("") +
    scale_x_log10(breaks = c(1,5,10,20,30,40)) + #breaks = c(1,2,3,4,5,6,7,8)
    #ylab("Heterozygosity-excess") +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.1,0.1,0.2,0.2), "cm") ,
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
        legend.direction = "vertical",
        legend.position = c(0.85,0.25),
        legend.title=element_text(size=10),
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))
    ) +
    scale_color_distiller(palette = "RdBu",
        direction = -1,
        name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    
    guides(color = guide_colorbar(barwidth = 0.5, barheight = 5, 
        title.position = "left")) + #, label.position = "bottom"
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) +
    geom_text_repel(label = short,size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
        segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
        segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines"))
p1

p2 <-  ggplot(aes(BreedingType, TPM80_ratio), data = stats_mod) +
    geom_boxplot(alpha = 0.3, col = "black",  size = 0.2, width = 0.4, aes(fill = BreedingType)) + #
    geom_point(size = 3.5, alpha = 0.6, aes(color = BreedingType)) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab("Breeding Habitat") +
    ylab("Heterozygosity-excess") +
    guides(fill=FALSE, color = FALSE) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.1,0.3,0.3,0.2), "cm") ,
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    ) 
# geom_text_repel(label = short,size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
#        segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
#       segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines"))
p2

plot_grid(p1, p2)


# beta coefficients
mod_out <- mod_beta[-1, c("components", "post_median", "lower", "upper")]
names(mod_out) <- c("comps", "pe", "cilow", "cihigh")

p3 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out) + 
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
        plot.margin = unit(c(1,0.2,0.33,0.1), "cm"),
        axis.title.x=element_text(margin=margin(t=12))) +
    scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2)) +
    scale_y_discrete(labels = c("Breeding\nhabitat",
        "SSD")) +
    xlab(expression(paste("Effect size ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p3

# SC
mod_out_SC <- mod_SC[, c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC ) <- c("comps", "pe", "cilow", "cihigh")

# structure coefficients
p4 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_SC) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3.5, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(1,0.2,0.4,0.3), "cm"),
        axis.text.x=element_text(margin=margin(t=0.5))) +
    scale_y_discrete(labels = c("Breeding\nhabitat",
        "SSD")) +
    xlab(expression(paste("Structure coefficient", " r(", hat(Y),",x)") )) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p4

plot_grid(p3, p4)

# R2
col_legend <- "#969696"
mod_out_R2 <- mod_R2[, c("combinations", "medianR2", "lower", "upper")][-4, ]
names(mod_out_R2) <- c("comps", "pe", "cilow", "cihigh")
mod_out_R2$comps <- factor(mod_out_R2$comps, levels = c("BreedingType", "SSD", "full model"))
# mod_out_R2$comps <- rev(fct_inorder(factor(mod_out_R2$comps)))
# structure coefficients
p5 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_R2 ) + 
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
        plot.margin = unit(c(0.3, 1, 0.6, 0.4), "cm")) +
    scale_y_discrete(labels = c( "Breeding\nhabitat", "SSD", "Full model")) +
    xlab(expression(paste(R^{2}))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
    #annotate("segment", x = 0.9, xend = 0.9, y = 0.8, yend = 3.5, color = col_legend) +
    #annotate("text", x = 0.98, xend = 0.95, y = 2, yend = 3, color = col_legend, label = c("common"), angle = 270) +
    annotate("segment", x = 1.2, xend = 1.2, y = 0.8, yend = 2.2, color = col_legend) +
    annotate("text", x = 1.3, xend = 1.3, y = 1.5, yend = 2.2, color = col_legend, label = c("unique"), angle = 270) +
    annotate("segment", x = 1.2, xend = 1.2, y = 2.4, yend = 3.5, color = col_legend) +
    annotate("text", x = 1.3, xend = 1.3, y = 2.9, yend = 3.5, color = col_legend, label = c("marginal"), angle = 270)
p5


p_top <- plot_grid(p2, p1, rel_widths = c(1.8, 3), labels = c("A", "B"))
p_top

p_bot <- plot_grid(p3, p4, p5, nrow = 1, rel_widths = c(1.5,1,1.6),
    labels = c("C","D","E"))
p_bot

p_final <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(1.6,1))
p_final

ggsave('figures/hetexc_vs_ecol.jpg',p_final,  width=9, height=6.5)

