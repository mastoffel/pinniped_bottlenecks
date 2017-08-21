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
    dplyr::select(TPM80_ratio, TPM70_ratio, TPM90_ratio, num_alleles_mean, obs_het_mean, 
        mean_allele_range, prop_low_afs_mean, mratio_mean, BreedingType,
        nloc, nind, bot) %>% 
    data.frame()

# load model output
mod_beta <- read_delim("data/processed/models/mod_gen_beta.txt", delim = " ")
mod_R2 <- read_delim("data/processed/models/mod_gen_R2.txt", delim = " ")
mod_SC <- read_delim("data/processed/models/mod_gen_SC.txt", delim = " ")


p1 <- ggplot(aes(prop_low_afs_mean, TPM80_ratio), data = stats_mod) +
    geom_smooth(method = "lm", col = "lightgrey",   alpha = 0.2) +
    geom_point(size = 3, alpha = 0.4) + # abc_out
    geom_point(size = 3, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    xlab("Proportion of low frequency alleles") +
    ylab("Heterozygosity-excess") +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.1,0.3,0.3,0.2), "cm")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
        )
p1

col_legend <- "#969696"

mod_out_R2 <- mod_R2[, c("combinations", "medianR2", "lower", "upper")]
names(mod_out_R2) <- c("comps", "pe", "cilow", "cihigh")
mod_out_R2 <- mod_out_R2[nrow(mod_out_R2):1, ]
mod_out_R2$comps <- fct_inorder(factor(mod_out_R2$comps))
mod_out_R2 <- mod_out_R2[-1, ]
# structure coefficients
p2 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_R2 ) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = c(0.5), color = "#333333", family = "Lato Light"),
        plot.margin = unit(c(0.3, 0.2, 0.3, 0.2), "cm")) +
    scale_y_discrete(labels = rev(c("Full model", "AR", "Het", "LFA","AR & Het", "AR & LFA", "Het & LFA"))) +
    xlab(expression(paste(R^{2}))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
    annotate("segment", x = 0.9, xend = 0.9, y = 0.8, yend = 3.2, color = col_legend) +
    annotate("text", x = 0.98, xend = 0.95, y = 2, yend = 3, color = col_legend, label = c("common"), angle = 270) +
    annotate("segment", x = 0.9, xend = 0.9, y = 3.6, yend = 6.4, color = col_legend) +
    annotate("text", x = 0.98, xend = 0.95, y = 5, yend = 6.5, color = col_legend, label = c("unique"), angle = 270) +
    annotate("segment", x = 0.9, xend = 0.9, y = 6.65, yend = 7.5, color = col_legend) +
    annotate("text", x = 0.98, xend = 0.95, y = 7.06, yend = 8, color = col_legend, label = c("marginal"), angle = 270)
p2

plot_grid(p1, p2, rel_widths = c(1.4,1))


# beta coefficients
mod_out <- mod_beta[-1, c("components", "post_median", "lower", "upper")]
names(mod_out) <- c("comps", "pe", "cilow", "cihigh")

p3 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out) + 
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
        plot.margin = unit(c(1,0,0.5,1.5), "cm"),
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))) +
    #scale_y_discrete(labels = c("Breeding\nhabitat\nice vs. land", "log(Abundance)", 
    #    "Sexual Size\nDimorphism")) +
    xlab(expression(paste("Effect size ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p3



mod_out_SC <- mod_SC[, c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC ) <- c("comps", "pe", "cilow", "cihigh")

# structure coefficients
p4 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_SC) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
       # axis.text.y = element_text(hjust = c(0.5), margin = margin(t = 0, r = 50, b = 0, l = 0)),
        axis.text.y = element_blank(),
        plot.margin = unit(c(1,0.2,0.65,0), "cm")) +
    scale_y_discrete(labels = c("Het", "AR", "LFA")) +
    xlab(expression(paste("Structure coefficient", " r(", hat(Y),",x)") )) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p4

# label plot
df <- data.frame(label = c("LFA", "AR", "Het"), value = c(1,1,1))
p5 <- ggplot(df, aes(value, label, label = label)) +
    geom_text(family = "Lato Light", color = '#333333') +
    theme_nothing() +
    theme(plot.margin = unit(c(1,-1,2,-1), "cm"))
  

p_top <- plot_grid(p1, p2, rel_widths = c(1.6,1), labels = c("A", "B"),label_fontfamily = "Lato",
    label_x = 0.1, label_y = 0.98)
p_top 

p_bot <- plot_grid(p3,p5,p4, labels = c("C","", ""),label_fontfamily = "Lato",
    label_x = 0.13, label_y = 1, rel_widths = c(1.3, 0.6, 1), ncol = 3)
p_bot

p_bot <- ggdraw(p_bot) + draw_label("D", colour = "black", x = 0.68, y = 0.97, fontfamily = "Lato", fontface = "bold")

p_final <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(1.4,1))
p_final

ggsave('figures/gen_div_vs_hetexc.jpg',p_final,  width=6.7, height=5.5)








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
          axis.title.y = element_blank()) +
    scale_y_discrete(labels = c("Allelic\nrichness", "Observed\nheterozygosity", 
        "Proportion of\nlow frequency\nalleles")) +
    xlab(expression(paste("Effect size ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p4

# structure coefficients
mod_out_SC <- mod_SC[c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC) <- c("comps", "pe", "cilow", "cihigh")

p5 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_SC) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank()) +
    scale_y_discrete(labels = c("Allelic\nrichness", "Observed\nheterozygosity", 
        "Proportion of\nlow frequency\nalleles")) +
    xlab(expression(paste("Structure coefficient", " r(", hat(Y),",x)") )) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p5    

mod_out_R2 <- mod_R2[c("combinations", "medianR2", "lower", "upper")]
names(mod_out_R2) <- c("comps", "pe", "cilow", "cihigh")

# commonality coefficients
# p6 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_R2) + 
#     # geom_point(size = 3, color = "grey69") + # abc_out
#     geom_errorbarh(alpha=0.4, color="black",height = 0) +
#     geom_point(size = 3, shape = 21, col = "black", fill = "grey69") +
#     # geom_errorbarh(alpha=0.4, color="black",height = 0) +
#     theme_martin() +
#     theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.line.x = element_line(color = '#333333'),
#         axis.title.y = element_blank()) +
#    # scale_y_discrete(labels = rev(mod_out_R2$comps)) +
#     xlab("R2") +
#     geom_vline(xintercept = 0, color = "black", alpha = 0.1)
# p6 

p_effect <- plot_grid(p4, p5, ncol = 1, labels = c("B", "C"), label_fontfamily = "Lato",
    label_x = 0.9, label_y = 0.95)
p_effect

p_final <- plot_grid(p3, p_effect, ncol = 2, rel_widths = c(2.4,1.6), labels = c("A", ""), label_fontfamily = "Lato",
    label_x = 0.9, label_y = 0.964)
p_final


ggsave('figures/gen_div_simple.jpg',p_final,  width=6, height=3)









# different
p3 <- ggplot(aes(prop_low_afs_mean, TPM70_ratio), data = stats_mod) +
    geom_smooth(method = "lm", col = "lightgrey", fill = "grey",  alpha = 0.2) +
    geom_point(size = 3, alpha = 0.4) + # abc_out
    geom_point(size = 3, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    xlab("Proportion of low frequency alleles") +
    ylab("Heterozygosity-excess") +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(1,1,0.5,1), "cm")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    )
p3



mod_out <- mod_beta[-1, c("components", "post_median", "lower", "upper")]
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
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0.15,0.2,1,1), "cm")) +
    scale_y_discrete(labels = c("Allelic\nrichness", "Observed\nheterozygosity", 
        "Proportion of\nlow frequency\nalleles")) +
    xlab(expression(paste("Effect size ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p4

# structure coefficients
mod_out_SC <- mod_SC[c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC) <- c("comps", "pe", "cilow", "cihigh")

p5 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_SC) + 
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
        plot.margin = unit(c(0.15,1,0.9,0.2), "cm")) +
    scale_y_discrete(labels = c("Allelic\nrichness", "Observed\nheterozygosity", 
        "Proportion of\nlow frequency\nalleles")) +
    xlab(expression(paste("Structure coefficient", " r(", hat(Y),",x)") )) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p5   


p_effect <- plot_grid(p4, p5, ncol = 2, rel_widths = c(1, 1.2), labels = c("B", "C"), label_fontfamily = "Lato",
    label_x = 0.83, label_y = 0.98)
p_effect
p_final <- plot_grid(p3,p_effect, ncol = 1, rel_heights = c(1.3,1), labels = c("A", ""), label_fontfamily = "Lato",
    label_x = 0.9, label_y = 0.964)
p_final

ggsave('figures/gen_div_simple_vt.jpg',p_final,  width=5, height=5.5)
