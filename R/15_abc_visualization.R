# plot posterior distributions 

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(ggthemr)
library(cowplot)
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
# plot nbot
p1 <- ggplot(data = filter(abc, pars == "nbot"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(stat="density", adjust = 2, size = 1, col = "cornflowerblue") +
    geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 0.5, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 10, labeller = as_labeller(species_names)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
    scale_x_continuous(limits = c(-20, 1000)) +
    xlab("bottleneck Ne") +
    ylab("posterior density")


p2 <- ggplot(data = filter(abc, pars == "mut_rate"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(stat="density", adjust = 2, size = 1, col = "cornflowerblue") +
    geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 1, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 10, labeller = as_labeller(species_names)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
    scale_x_continuous(limits = c(-0.0001, 6e-04)) +
    xlab("mutation rate") +
    ylab("posterior density")

p3 <- ggplot(data = filter(abc, pars == "gsm_param"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(stat="density", adjust = 2, size = 1, col = "cornflowerblue") +
    geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 1, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 10, labeller = as_labeller(species_names)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
    # scale_x_continuous(limits = c(0, 600)) +
    xlab("prop. multistep mutations") +
    ylab("posterior density")

plot_grid(p1, p2, p3, ncol = 3)


ggplot(data = filter(abc, pars == "tbotend"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(stat="density", adjust = 2, size = 1, col = "cornflowerblue") +
    geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 1, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 5, labeller = as_labeller(species_names)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
    # scale_x_continuous(limits = c(0, 600)) +
    xlab("bottleneck end (in gen. ago)") +
    ylab("posterior density")

ggplot(data = filter(abc, pars == "tbotstart"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(stat="density", adjust = 2, size = 1, col = "cornflowerblue") +
    geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    geom_hline(yintercept = 0, linetype="twodash", col = "black", size = 1, alpha = 0.5) + 
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 5, labeller = as_labeller(species_names)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
    # scale_x_continuous(limits = c(0, 600)) +
    xlab("bottleneck start (in gen. ago)") +
    ylab("posterior density")

pop_size_prior <- data.frame(pop_size_prior = round(rlnorm(10000, 10.5, 1)))
ggplot(data = filter(abc, pars == "nhist"), aes(x = adj_vals)) +
    # geom_density(adjust = 1.5) +
    geom_line(stat="density", adjust = 2, size = 1, col = "cornflowerblue") +
    geom_line(aes(x=unadj_vals), stat="density", adjust = 2, size = 0.5, col = "goldenrod") +
    geom_line(data = pop_size_prior, aes(x=pop_size_prior), stat="density", adjust = 1, size = 0.5) +
    # geom_histogram() + 
    facet_wrap(~ species, scales = "free_y", nrow = 5, labeller = as_labeller(species_names)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
    scale_x_continuous(limits = c(0, 150000)) +
    xlab("historical population size") +
    ylab("posterior density")
