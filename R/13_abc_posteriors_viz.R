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
library(ggridges)
source("R/martin.R")
library(patchwork)   
library(ggforce)
library(extrafont)
library(extrafontdb)
## bottleneck posteriors -----

# load abc posterior data
load("data/processed/abc_estimates/abc_10000k_bot_complete.RData")
# load parameter distributions
# abc_params <- fread("data/processed/abc_estimates/sims_1500k_params.txt")

# unnest the abc posterior values for plotting
abc_bot <- unnest(abc_complete)
# filter(abc, pars == "nbot")

species_names_bot <- c(
    "antarctic_fur_seal" = "Antarctic Fur Seal",
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

species_names_bot_twolines <- c(
    "antarctic_fur_seal" = "Antarctic\nFur Seal",
    "california_sea_lion" = "California\nSea Lion",
    "galapagos_fur_seal" = "Galapagos\nFur Seal", 
    "grey_seal_orkneys" = "Grey Seal",
    "hawaiian_monk_seal" = "Hawaiian\nMonk Seal",
    "lagoda_ringed_seal" = "Ladoga\nRinged Seal",
    "mediterranean_monk_seal" = "Mediterranean\nMonk Seal",
    "nes" = "Northern\nElephant Seal",
    "saimaa_ringed_seal" = "Saimaa\nRinged Seal",
    "south_american_fur_seal" = "South American\nFur Seal"
)

# abc$species <- factor(abc$species, levels = rev(c("saimaa_ringed_seal", "mediterranean_monk_seal","hawaiian_monk_seal", "nes", "galapagos_fur_seal",
#    "lagoda_ringed_seal", "antarctic_fur_seal", "grey_seal_orkneys", "california_sea_lion",  "south_american_fur_seal")))

abc_bot$species <- factor(abc_bot$species, levels = c("saimaa_ringed_seal", "mediterranean_monk_seal","hawaiian_monk_seal", "nes", "galapagos_fur_seal",
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



## really cool way of plotting with beeswarm plot and boxplot
library(ggbeeswarm)
estimate_mode <- function(s) {
    d <- density(s, adjust = 2.5)
    d$x[which.max(d$y)]
}

#library(png)
#library(grid)
#img <- readPNG("other_stuff/AFS.png")
#g <- rasterGrob(img, interpolate=TRUE)

p <- abc_bot %>% filter(pars == "nbot") %>% 
    group_by(species) %>% 
    #sample_n(4000) %>% 
    #filter(!(species == "hawaiian_monk_seal" | species == "saimaa_ringed_seal" | species == "mediterranean_monk_seal")) %>% 
    ggplot(aes(y = adj_vals, x = species)) + 
    # geom_quasirandom(alpha = 0.01, size = 2, color = "#053061", width = 0.45, bandwidth = 1.5) +
    geom_quasirandom(alpha = 0.07, size = 1, color = "#053061", width = 0.47, bandwidth = 2.5) +
   # geom_beeswarm(priority='density', alpha = 0.5, cex = 0.2, color = "grey") +
    # geom_jitter(size = 0.5, alpha = 0.1, width = 0.2, color = "grey") +
    geom_boxplot(width = 0.4, outlier.shape = NA, color = "white", alpha = 0.5, size = 0.2) +
    stat_summary(fun.y = "estimate_mode", colour = "black", geom = "point", size = 2, shape = 21, fill = "grey") +
    #stat_summary(fun.y = "mean", colour = "blue", geom = "point") +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
   # scale_fill_cyclical(values = c("lightgrey", "darkgrey")) +
    #ylim(0, 900) +
    scale_x_discrete(labels = species_names_bot_twolines) + 
    scale_y_continuous(breaks = c(seq(from = 0, to = 800, by = 200)), limits = c(0,900)) +
    #scale_y_discrete(expand = c(0.01, 0))+
    xlab("") +
    ylab(expression(Bottleneck~N[e])) +
    coord_flip() +
    theme(plot.margin = unit(c(0.5,1.5,0.5,0), "cm"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 11),
        axis.title.y = element_text(size = 13),
        axis.text.y = element_text(size = 11))

p

# plot as pdf
ggplot2::ggsave(filename = "other_stuff/figures/abc_posteriors_vert.pdf", p,
    width = 3.25, height = 10, device = "pdf")

Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts("other_stuff/figures/abc_posteriors.pdf")

# plot as jpg
ggplot2::ggsave(filename = "other_stuff/figures/abc_posteriors_vert.jpg", p,
    width = 4, height = 6)



# SUPPLEMENTARY plots: Mut Rate

p_mut_bot <- abc_bot %>% filter(pars == "mut_rate") %>%
    ggplot(aes(adj_vals, y = species, fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, scale = 3, alpha = 0.7) +
    scale_fill_cyclical(values = c("#4040B0", "#9090F0"), guide = "legend") +
    theme_martin(legend.position='none', base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    xlab(expression(mutation~rate~mu~(x~10^-4))) +
    scale_y_discrete(labels = species_names_bot) +
    scale_x_continuous(limits = c(-0.00005, 0.0006), breaks = c(0, 0.0001, 0.0002, 
                                    0.0003, 0.0004, 0.0005), labels = c(0,1,2,3,4,5)) +
    ylab("") +
    ggtitle("Bottleneck model")
p_mut_bot

# ggsave(filename = "other_stuff/figures/figures_final/Sup_abc_posteriors_mutrate.jpg", plot = p_mut, width = 4.5, height = 6.5)


# neutral model posteriors

# load abc posterior data
load("data/processed/abc_estimates/abc_10000k_neut_complete.RData")
# load parameter distributions
# abc_params <- fread("data/processed/abc_estimates/sims_1500k_params.txt")

# unnest the abc posterior values for plotting
abc_neut <- unnest(abc_complete)
# filter(abc, pars == "nbot")

species_names_neut <- c(
    "arctic_ringed_seal" = "Ringed Seal",
    "atlantic_walrus"  = "Walrus",
    "australian_fur_seal" = "Australian Fur Seal",
    "baltic_ringed_seal" = "Baltic Ringed Seal",
    "bearded_seal" = "Bearded Seal",
    "crabeater_seal" = "Crabeater Seal",
    "galapagos_sea_lion" = "Galapagos Sea Lion",
    "harbour_seal_waddensee" = "Harbour Seal",
    "hooded_seal" = "Hooded Seal",
    "leopard_seal" = "Leopard Seal",
    "new_zealand_fur_seal" = "New Zealand Fur Seal",
    "new_zealand_sea_lion" = "New Zealand Sea Lion",
    "northern_fur_seal" = "Northern Fur Seal",
    "ross_seal" = "Ross Seal",
    "ses" = "Southern Elephant Seal",
    "south_american_sea_lion" = "South American Sea Lion",
    "stellers_sea_lion" =  "Steller Sea Lion",
    "weddell_seal" = "Weddel Seal"
)

#abc$species <- factor(abc$species, levels = c("saimaa_ringed_seal", "mediterranean_monk_seal","hawaiian_monk_seal", "nes", "galapagos_fur_seal",
#    "lagoda_ringed_seal", "antarctic_fur_seal", "grey_seal_orkneys", "california_sea_lion",  "south_american_fur_seal"))

p_mut_neut <- abc_neut %>% filter(pars == "mut_rate") %>%
    ggplot(aes(adj_vals, y = species, fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, scale = 3, alpha = 0.8) +
    scale_fill_cyclical(values = c("#4040B0", "#9090F0"), guide = "legend") +
    theme_martin(legend.position='none', base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    xlab(expression(mutation~rate~mu~(x~10^-4))) +
    scale_y_discrete(labels = species_names_neut) +
    scale_x_continuous(limits = c(-0.00005, 0.0006), breaks = c(0, 0.0001, 0.0002, 0.0003, 
                        0.0004, 0.0005), labels = c(0,1,2,3,4,5)) +
    ylab("") +
    ggtitle("Neutral model")

 

## final figure
p_final <- p_mut_bot + p_mut_neut
p_final
ggsave(filename = "other_stuff/figures/figures_final/Sup_abc_posteriors_mutrate.jpg", 
    plot = p_final, width = 7.5, height = 5.5)


ggsave(filename = "other_stuff/figures/figures_final/Sup_abc_posteriors_mutrate.pdf", 
    plot = p_final, width = 7.5, height = 5.5)

Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts("other_stuff/figures/figures_final/Sup_abc_posteriors_mutrate.pdf")






### for GSM

p_gsm_neut <- abc_neut %>% filter(pars == "gsm_param") %>%
    ggplot(aes(adj_vals, y = species, fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, scale = 3, alpha = 0.8) +
    scale_fill_cyclical(values = c("#4040B0", "#9090F0"), guide = "legend") +
    theme_martin(legend.position='none', base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    xlab(expression(Multistep~mutations~"("~GSM[par]~")")) +
    scale_y_discrete(labels = species_names_neut) +
    scale_x_continuous(limits = c(-0.05, 0.4), breaks = c(0, 0.1, 0.2, 0.3,0.4)) +
    ylab("") 

ggsave(filename = "other_stuff/figures/figures_final/Sup_abc_posteriors_gsm.jpg", 
    plot = p_gsm_neut, width = 4, height = 6)

ggsave(filename = "other_stuff/figures/figures_final/Sup_abc_posteriors_gsm.pdf", 
    plot = p_gsm_neut, width = 4, height = 6)

Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts("other_stuff/figures/figures_final/Sup_abc_posteriors_gsm.pdf")

























