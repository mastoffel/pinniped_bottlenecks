# plot posterior distributions 

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggthemr)
library(cowplot)
library('ggthemes')

ggthemr("fresh", text_size = 30,spacing = 4, layout = "clean")
# ggthemr("fresh", text_size = 12, spacing = 2, layout = "clean")
# ggthemr_reset()

# load abc posterior data
load("data/processed/abc_estimates/abc_5000k_complete.RData")
# load parameter distributions
# abc_params <- fread("data/processed/abc_estimates/sims_1500k_params.txt")

# unnest the abc posterior values for plotting
abc <- unnest(abc_complete)
# filter(abc, pars == "nbot")

species_names <- c(
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
# plot nbot
p1 <- ggplot(data = filter(abc, pars == "nbot"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(aes(y=..scaled..), stat="density", adjust = 2, size = 1, col = "#99000d") +
    # geom_line(aes(x=unadj_vals, y=..scaled..), stat="density", adjust = 2, size = 0.5, col = "grey", linetype = 2) +
    # geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 0.5, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, nrow = 10, labeller = as_labeller(empty_names)) + # scales = "free_y" #as_labeller(species_names)
    theme_tufte(ticks = TRUE, base_family = "Helvetica", base_size = 16) +
    # theme_classic(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           strip.background = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           strip.text = element_text() ,
           axis.line.y = element_blank(),
           plot.margin=unit(c(5.5, 15, 15, 15), "points")
           )+
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    scale_x_continuous(limits = c(-20, 1000)) +
    #xlab(expression(bold('Bottleneck N'[e]))) +
    xlab("Bottleneck Ne") +
   # xlab(paste0(c("bottleneck", bquote(atop(N[e]))))) +
    #xlab(paste0("bottleneck ", expression(N[e])), "\n  ") +
    ylab("scaled posterior density")


p2 <- ggplot(data = filter(abc, pars == "mut_rate"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(aes(y=..scaled..), stat="density", adjust = 2, size = 1, col = "#3690c0") +
    #geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    #geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 1, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 10, labeller = as_labeller(species_names)) +
    theme_tufte(ticks = TRUE, base_family = "Helvetica", base_size = 16) +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(5.5, 15 , 15, 15), "points")
        ) +
    scale_x_continuous(limits = c(-0.0001, 6e-04), breaks = c( 0.0001, 0.0003, 0.0005), 
        labels =  c(expression(1%*%10^{-4}), expression(3%*%10^{-4}), expression(5%*%10^{-4}))) +
    xlab("mutation rate") 
    # ylab("posterior density")

p3 <- ggplot(data = filter(abc, pars == "gsm_param"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(aes(y=..scaled..), stat="density", adjust = 2, size = 1, col = "#6a51a3") +
    #geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    #geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 1, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 10, labeller = as_labeller(empty_names)) +
    theme_tufte(ticks = TRUE, base_family = "Helvetica", base_size = 16) +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin=unit(c(5.5, 15 , 15, 15), "points")
        ) +
    # scale_x_continuous(limits = c(0, 600)) +
    xlab("% multistep mutations") 
    # ylab("posterior density")

p_final <- plot_grid(p1, p2, p3, ncol = 3)

save_plot("abc_plot.pdf", p_final,
    ncol = 2, # we're saving a grid plot of 2 columns
    nrow = 1, # and 2 rows
    # each individual subplot should have an aspect ratio of 1.3
    # base_aspect_ratio = 0.9,
    base_height = 12,
    base_width = 4.5
)



## try 2: mean and CI

# mode and HDI
estimate_mode <- function(s) {
    d <- density(s)
    d$x[which.max(d$y)]
}

library(coda)
abc_summed <- abc %>% 
    group_by(species, pars) %>% 
    summarise(mode_par = estimate_mode(adj_vals),
        HPD_lod = HPDinterval(mcmc(adj_vals), 0.95)[1],
        HPD_high = HPDinterval(mcmc(adj_vals), 0.95)[2])

data_summary <- function(adj_vals) {
    mode_par = estimate_mode(adj_vals)
    HPD_lod = HPDinterval(mcmc(adj_vals), 0.95)[1]
    HPD_high = HPDinterval(mcmc(adj_vals), 0.95)[2]
    return(c(x = mode_par, xmin = HPD_lod, xmax = HPD_high))
}
ggthemr("fresh", text_size = 10,spacing = 1, layout = "clean")


p1 <- ggplot(data = filter(abc, pars == "nbot"), aes(x = adj_vals, y = species)) +
    # geom_density(adjust = 1.5) +
    geom_line(aes(y=..scaled..), stat="density", adjust = 2, size = 1, col = "#99000d") +
    facet_wrap(~ species, nrow = 10) +
    stat_summary(fun.data=data_summary)

afs <- filter(abc, (pars == "nbot") & (species == "antarctic_fur_seal"))

ggplot(afs, aes(adj_vals)) + geom_violin()

p1 <- ggplot(data = filter(abc, pars == "nbot"), aes(x = species, y = adj_vals)) +
    geom_violin(trim = FALSE) +
    coord_flip()
    facet_wrap(~ species, nrow = 10)
    geom_line(aes(y=..scaled..), stat="density", adjust = 2, size = 1, col = "#99000d") +
    facet_wrap(~ species, nrow = 10) +
    stat_summary(fun.data=data_summary)
    + stat_summary(fun.data=data_summary) +
    
    geom_errorbarh(data = filter(abc_summed, pars == "nbot"), aes(x = mean_par, xmax = lower_ci_par, xmin = upper_ci_par)) +
# geom_line(aes(x=unadj_vals, y=..scaled..), stat="density", adjust = 2, size = 0.5, col = "grey", linetype = 2) +
# geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 0.5, alpha = 0.5) + 
# geom_histogram() + 
facet_wrap(~ species, nrow = 10, labeller = as_labeller(empty_names)) + # scales = "free_y" #as_labeller(species_names)
    theme_tufte(ticks = TRUE, base_family = "Helvetica", base_size = 16) +
    # theme_classic(base_family = "Helvetica") +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text() ,
        axis.line.y = element_blank(),
        plot.margin=unit(c(5.5, 15, 15, 15), "points")
    )+
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", ".5", "1")) +
    scale_x_continuous(limits = c(-20, 1000)) +
    #xlab(expression(bold('Bottleneck N'[e]))) +
    xlab("Bottleneck Ne") +
    # xlab(paste0(c("bottleneck", bquote(atop(N[e]))))) +
    #xlab(paste0("bottleneck ", expression(N[e])), "\n  ") +
    ylab("scaled posterior density")


