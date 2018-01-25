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
show_col(magma(20))

# ggthemr("fresh", text_size = 30,spacing = 4, layout = "clean")
# ggthemr("fresh", text_size = 12, spacing = 2, layout = "clean")
# ggthemr_reset()

# load abc posterior data
load("data/processed/abc_estimates/abc_10000k_bot_complete.RData")
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

head(abc)

# abc %>% group_by(species) %>% summarise(mean_val = mean(adj_val), median_val = median(adj_val)) 

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

#c ggthemr("fresh", text_size = 12,spacing = 4, layout = "clean")

create_dens_plot <- function(species, parameter){
    # subset summary data
    df_summed <- filter(abc_summed , (pars == !!parameter) & (species == !!species))
    # subset posterior data
    temp_data <- filter(abc, (pars == !!parameter) & (species == !!species))
    # create density object and extract maximum density
    p_temp <- ggplot(temp_data, aes(x=adj_vals)) + geom_density(adjust = 2)
    dens_obj <- ggplot_build(p_temp)
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

#### all this code is for plotting all three parameter distributions
# p_nbot <- do.call(plot_grid, c(all_plots_nbot, ncol = 1, list(rel_heights = c(rep(1,9), 1.2))))
# # mutation plots
# all_plots_mut <- lapply(levels(abc$species), create_dens_plot, "mut_rate")
# p_mut <- do.call(plot_grid, c(all_plots_mut, ncol = 1, list(rel_heights = c(rep(1,9), 1.2))))
# # gsm plots
# all_plots_gsm <- lapply(levels(abc$species), create_dens_plot, "gsm_param")
# p_gsm <- do.call(plot_grid, c(all_plots_gsm, ncol = 1, list(rel_heights = c(rep(1,9), 1.2))))
# 
# # 
# all_plots <- c(all_plots_nbot, all_plots_mut, all_plots_gsm)
# 
# # resort
# inds <- unlist(lapply(1:10, function(x) out <- c(x, x+10, x+20)))
# p_all <- do.call(plot_grid, c(all_plots[inds], ncol = 3, list(rel_heights = c(rep(1,9), 1.3))))
# p_all

# maybe just the bottlenecks
inds <- c(1,6,2,7,3,8,4,9,5,10)
p_nbot <- do.call(plot_grid, c(all_plots_nbot[inds], ncol = 2, list(rel_heights = c(rep(1,9), 1.25))))
p_nbot

p <- ggdraw(add_sub(p_nbot, "Bottleneck Ne", size = 13))

#p + draw_label("posterior density", x = 0, y = 0.5,
#    vjust = 0, hjust = 0, size = 10, fontface = 'bold', angle = 90)

p_test <- ggdraw(add_sub(p, "posterior density",vpadding=grid::unit(-4, "lines"),
          x = 0, y = 14, angle = 90, size = 10, vjust = -0.2
    )) #vpadding=grid::unit(-5, "lines")
p_test

save_plot("abc_plot_2.pdf", p_test,
    ncol = 2, # we're saving a grid plot of 2 columns
    nrow = 1, # and 2 rows
    # each individual subplot should have an aspect ratio of 1.3
    # base_aspect_ratio = 0.9,
    base_height = 6,
    base_width = 2
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
p1
