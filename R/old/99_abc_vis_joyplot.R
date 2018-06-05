# plot posterior distributions 

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggthemr)
library(cowplot)
library('ggthemes')
library(scales)
library(viridis)
library(ggjoy)
source("martin.R")
# load abc posterior data
load("data/processed/abc_estimates/abc_10000k_complete_30.RData")
# load parameter distributions
# abc_params <- fread("data/processed/abc_estimates/sims_1500k_params.txt")

# unnest the abc posterior values for plotting
abc <- unnest(abc_complete)
# filter(abc, pars == "nbot")

species_names <- c(
    "antarctic_fur_seal" = "Antarctic Fur Seal",
    "guadalupe_fur_seal" = "Guadalupe Fur Seal",
    "california_sea_lion" = "California Sea Lion",
    "galapagos_fur_seal" = "Galapagos Fur Seal", 
    "grey_seal_orkneys" = "Grey Seal",
    "hawaiian_monk_seal" = "Hawaiian Monk Seal",
    "lagoda_ringed_seal" = "Ladoga Ringed Seal",
    "mediterranean_monk_seal" = "Mediterranean Monk Seal",
    "nes" = "Northern Elephant Seal",
    "saimaa_ringed_seal" = "Saimaa Ringed Seal",
    "south_american_fur_seal" = "South American Fur Seal"
)

abc$species <- factor(abc$species, levels = rev(c("saimaa_ringed_seal", "mediterranean_monk_seal","hawaiian_monk_seal", "nes",
    "guadalupe_fur_seal", "galapagos_fur_seal",
    "lagoda_ringed_seal", "antarctic_fur_seal", "grey_seal_orkneys", "california_sea_lion",  "south_american_fur_seal")))

abc$species <- factor(abc$species, levels = c("saimaa_ringed_seal", "mediterranean_monk_seal","hawaiian_monk_seal", "nes", "galapagos_fur_seal",
    "lagoda_ringed_seal", "antarctic_fur_seal", "grey_seal_orkneys", "california_sea_lion",  "south_american_fur_seal"))

empty_names <- c(
    "antarctic_fur_seal" = " ",
    "california_sea_lion" = " ",
    "galapagos_fur_seal" = " ", 
    "grey_seal_orkneys" = " ",
    "hawaiian_monk_seal" = " ",
    "lagoda_ringed_seal" = " ",
    "mediterranean_monk_seal" = " ",
    "nes" = " ",
    "saimaa_ringed_seal" = " ",
    "south_american_fur_seal" = " "
)
source("R/martin.R")


p1 <- abc %>% filter(pars == "nbot") %>% 
        filter(species == "hawaiian_monk_seal" | species == "nes" | species == "saimaa_ringed_seal" | species == "mediterranean_monk_seal") %>% 
        ggplot(aes(x = unadj_vals, y = species, fill = species,  height = ..density..)) + 
        geom_joy(stat = "density",scale = 3, alpha = 1,  adjust = 2) +
        theme_martin() +
        scale_fill_cyclical(values = c("#deebf7", "#9ecae1")) +
        scale_fill_cyclical(values = c("black", "white")) +
        xlim(0, 300) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()
              ) +
        scale_y_discrete(expand = c(0.01, 0), labels = c("Northern Elephant Seal", "Hawaiian Monk Seal", "Mediterranean Monk Seal",
                                                         "Saimaa Ringed Seal"))+
        xlab("")

p2 <- abc %>% filter(pars == "nbot") %>% 
    filter(!(species == "hawaiian_monk_seal" | species == "nes" | species == "saimaa_ringed_seal" | species == "mediterranean_monk_seal")) %>% 
    ggplot(aes(x = unadj_vals, y = species, fill = species,  height = ..density..)) + 
    geom_joy(scale = 2, alpha = 1, stat = "density", adjust = 2) +
    theme_martin() +
    scale_fill_cyclical(values = c("lightgrey", "darkgrey")) +
    xlim(0, 1000) +
    scale_y_discrete(expand = c(0.01, 0), labels = c("South American Fur Seal", "California Sea Lion", "Grey Seal", "Antarctic Fur Seal",
        "Ladoga Ringed Seal","Galapagos Fur Seal" ))+
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    xlab("Bottleneck Ne")

    
plot_grid(p1, p2, ncol = 1, rel_heights = c(1.7, 2))


abc %>% filter(pars == "nbot") %>% 
    #filter(!(species == "hawaiian_monk_seal" | species == "saimaa_ringed_seal" | species == "mediterranean_monk_seal")) %>% 
    ggplot(aes(y = unadj_vals, y = species, fill = species,  height = ..density..)) + 
    geom_joy(scale = 4, alpha = 0.7, stat = "density", adjust = 2) +
    theme_martin() +
    scale_fill_cyclical(values = c("lightgrey", "darkgrey")) +
    xlim(0, 800) +
    scale_y_discrete(expand = c(0.01, 0))+
    xlab("")



## really cool way of plotting with beeswarm plot and boxplot
library(ggbeeswarm)
estimate_mode <- function(s) {
    d <- density(s, adjust = 2)
    d$x[which.max(d$y)]
}

p <- abc %>% filter(pars == "nbot") %>% 
    group_by(species) %>% 
    #sample_n(4000) %>% 
    #filter(!(species == "hawaiian_monk_seal" | species == "saimaa_ringed_seal" | species == "mediterranean_monk_seal")) %>% 
    ggplot(aes(y = adj_vals, x = species)) + 
    geom_quasirandom(alpha = 0.01, size = 2, color = "#053061", width = 0.45, bandwidth = 1.5) +
   # geom_beeswarm(priority='density', alpha = 0.5, cex = 0.2, color = "grey") +
    # geom_jitter(size = 0.5, alpha = 0.1, width = 0.2, color = "grey") +
    geom_boxplot(width = 0.4, outlier.shape = NA, color = "white", alpha = 0.5, size = 0.4) +
    stat_summary(fun.y = "estimate_mode", colour = "black", geom = "point", size = 2, shape = 21, fill = "grey") +
    #stat_summary(fun.y = "mean", colour = "blue", geom = "point") +
    theme_martin() +
   # scale_fill_cyclical(values = c("lightgrey", "darkgrey")) +
    #ylim(0, 900) +
    scale_x_discrete(labels = species_names) + 
    scale_y_continuous(breaks = c(seq(from = 100, to = 900, by = 200)), limits = c(0,900)) +
    #scale_y_discrete(expand = c(0.01, 0))+
    xlab("") +
    ylab(expression(Bottleneck~N[e])) +
    coord_flip() +
    theme(plot.margin = unit(c(0.5,1.5,0.5,0), "cm"))

p
ggsave(filename = "figures/abc_posteriors.jpg", plot = p, width = 3.5, height = 6.5)

abc %>% filter(pars == "nbot") %>% 
    group_by("species") %>% 
    #sample_n(1000) %>% 
    #filter(!(species == "hawaiian_monk_seal" | species == "saimaa_ringed_seal" | species == "mediterranean_monk_seal")) %>% 
    ggplot(aes(x= adj_vals)) + 
    #geom_density(adjust = 1.5) +
    geom_histogram(bins = 50)+
    facet_wrap(~species, scales = "free")
    #geom_quasirandom(alpha = 0.05, size = 2, color = "#4575b4", width = 0.45, bandwidth = 1) +
    # geom_beeswarm(alpha = 0.1, size = 0.1, color = "grey") +
    # geom_jitter(size = 0.5, alpha = 0.1, width = 0.2, color = "grey") +
    geom_boxplot(width = 0.4, outlier.shape = NA, color = "black", alpha = 0.5, size = 0.3) +
    stat_summary(fun.y = "estimate_mode", colour = "black", geom = "point", size = 2, shape = 21, fill = "grey") +
    #stat_summary(fun.y = "mean", colour = "blue", geom = "point") +
    theme_martin() +
    # scale_fill_cyclical(values = c("lightgrey", "darkgrey")) +
    #ylim(0, 900) +
    scale_x_discrete(labels = species_names) + 
    scale_y_continuous(breaks = c(seq(from = 0, to = 900, by = 100)), limits = c(0,900)) +
    #scale_y_discrete(expand = c(0.01, 0))+
    xlab("") +
    ylab(expression(Bottleneck~N[e])) +
    coord_flip() +
    theme(plot.margin = unit(c(0.5,1.5,0.5,0), "cm"))


