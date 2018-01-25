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
    "south_american_fur_seal" = "South\nAmerican\nFur Seal"
)

# abc$species <- factor(abc$species, levels = rev(c("saimaa_ringed_seal", "mediterranean_monk_seal","hawaiian_monk_seal", "nes", "galapagos_fur_seal",
#    "lagoda_ringed_seal", "antarctic_fur_seal", "grey_seal_orkneys", "california_sea_lion",  "south_american_fur_seal")))

abc_bot$species <- factor(abc_bot$species, levels = rev(c("saimaa_ringed_seal", "mediterranean_monk_seal","hawaiian_monk_seal", "nes", "galapagos_fur_seal",
    "lagoda_ringed_seal", "antarctic_fur_seal", "grey_seal_orkneys", "california_sea_lion",  "south_american_fur_seal")))



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
    theme_martin() +
   # scale_fill_cyclical(values = c("lightgrey", "darkgrey")) +
    #ylim(0, 900) +
    scale_x_discrete(labels = species_names_bot_twolines) + 
    scale_y_continuous(breaks = c(seq(from = 0, to = 800, by = 200)), limits = c(0,900)) +
    #scale_y_discrete(expand = c(0.01, 0))+
    xlab("") +
    ylab(expression(Bottleneck~N[e])) +
   # coord_flip() +
    theme(plot.margin = unit(c(0.5,1.5,0.5,0), "cm"),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))

p
ggsave(filename = "other_stuff/figures/abc_posteriors_try.jpg", plot = p, 
       width = 10, height = 3.5)


# SUPPLEMENTARY plots: Mut Rate

p_mut_bot <- abc_bot %>% filter(pars == "mut_rate") %>%
    ggplot(aes(adj_vals, y = species, fill = species)) +
    geom_density_ridges(rel_min_height = 0.01, scale = 3, alpha = 0.7) +
    scale_fill_cyclical(values = c("#4040B0", "#9090F0"), guide = "legend") +
    theme_martin(legend.position='none') +
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
    theme_martin(legend.position='none') +
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





# ridgeline plot new try
p_bot <- abc_bot %>% 
    filter(pars == "nbot") %>% 
    filter(species %in% c("galapagos_fur_seal", "lagoda_ringed_seal", "antarctic_fur_seal", 
        "grey_seal_orkneys", "california_sea_lion", "south_american_fur_seal")) %>%
    ggplot(aes(x = adj_vals, y = species, fill = ..x..)) +
    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01, gradient_lwd = 1.) +
    scale_fill_viridis(alpha = 0.5) +
    theme_martin(legend.position='none') +
    scale_y_discrete(labels = species_names_bot) + 
    scale_x_continuous(breaks = c(seq(from = 100, to = 900, by = 200)), limits = c(0,900)) +
    #scale_y_discrete(expand = c(0.01, 0))+
    xlab("") +
    ylab(expression(Bottleneck~N[e])) +
    theme(plot.margin = unit(c(0.5,1.5,0.5,0), "cm"))
    #geom_density_ridges(stat = "density", rel_min_height = 0.01, scale = 1.2, alpha = 0.7,
    #    panel_scaling = TRUE) +
    
    #scale_fill_cyclical(values = c("#4040B0", "#9090F0"), guide = "legend") +
 
p_bot
    
    
    xlab(expression(mutation~rate~mu~(x~10^-4))) +
    scale_y_discrete(labels = species_names_bot) +
    scale_x_continuous(limits = c(-0.00005, 0.0006), breaks = c(0, 0.0001, 0.0002, 
        0.0003, 0.0004, 0.0005), labels = c(0,1,2,3,4,5)) +
    ylab("") +
    ggtitle("Bottleneck model")
p_mut_bot






















hist(rlnorm(100, 10.5, 1))
# mode and HDI
estimate_mode <- function(s) {
    d <- density(s)
    d$x[which.max(d$y)]
}

library(coda)

# summary function

abc_summed <- abc %>% 
    group_by(species, pars) %>% 
    summarise(mode_par = estimate_mode(adj_vals),
        HPD_low_95 = HPDinterval(mcmc(adj_vals), 0.95)[1],
        HPD_high_95 = HPDinterval(mcmc(adj_vals), 0.95)[2],
        HPD_low_50 = HPDinterval(mcmc(adj_vals), 0.50)[1],
        HPD_high_50 = HPDinterval(mcmc(adj_vals), 0.50)[2])

# ggthemr("fresh", text_size = 12,spacing = 4, layout = "clean")

create_dens_plot <- function(species, parameter){
    # subset summary data
    df_summed <- filter(abc_summed , (pars == !!parameter) & (species == !!species))
    # subset posterior data
    temp_data <- filter(abc, (pars == !!parameter) & (species == !!species))
    # create density object and extract maximum density
    dens_obj <- print(ggplot(temp_data, aes(x=adj_vals)) + geom_density(adjust = 2))
    max_dens <- max(dens_obj$data[[1]]$density)
    # plotting
    p <- ggplot(data =temp_data, aes(adj_vals)) +
        geom_line(stat = "density", adjust = 2, color = "#969696") +
        geom_errorbarh(data = df_summed, aes(x = mode_par, xmin = HPD_low_95, xmax = HPD_high_95, y = max_dens/2),  
            inherit.aes = FALSE, height = 0, size = 0.4, color = "#6B1D81FF") +
        geom_errorbarh(data = df_summed, aes(x = mode_par, xmin = HPD_low_50, xmax = HPD_high_50, y = max_dens/2),  
            inherit.aes = FALSE, height = 0, size = 1, color = "#6B1D81FF") +
        geom_point(data = df_summed, aes(x  =mode_par, y = max_dens/2), size = 2, color = "#FCFDBFFF") + #"#FCFDBFFF"
        geom_point(data = df_summed, aes(x  =mode_par, y = max_dens/2), size = 2, color = "black", shape = 1) +
        # scale_y_continuous(breaks = c(max_dens/2), labels = df_summed$species) +
        theme_tufte(ticks = TRUE, base_family = "Helvetica", base_size = 16) +
        theme(axis.title.x=element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            plot.margin=unit(c(8, 5 , 5, 5), "points"),
            axis.line.x = element_line(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(size = 9, colour = "#525252")) 
    
    
    if(parameter == "nbot"){
        p <- p +  
            scale_x_continuous(limits = c(-20, 1000), breaks = c(0,100,300,500,700,900)) +
            ggtitle(species_names[as.character(df_summed$species)])
    } else if (parameter == "mut_rate"){
        p <- p + scale_x_continuous(limits = c(-0.0001, 6e-04), breaks = c( 0.0001, 0.0003, 0.0005), 
            labels =  c(expression(1%*%10^{-4}), expression(3%*%10^{-4}), expression(5%*%10^{-4}))) +
            ggtitle(" ")
    } else if (parameter == "gsm_param"){
        p <- p +  ggtitle(" ")
    }
    
    if(species == levels(abc$species)[length(levels(abc$species))] |
            species == levels(abc$species)[length(levels(abc$species))/2]){
        p <- p + theme(axis.text.x = element_text(size = 8, angle = 50),
            axis.line.x = element_line(),
            plot.margin=unit(c(8, 5, 0, 5), "points"))
    }
    p 
}
# bottleneck plots
all_plots_nbot <- lapply(levels(abc$species), create_dens_plot, "nbot")














## ridgeline plot

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