# plots for the talk
# install_github("mastoffel/sealABC", dependecies = TRUE)
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(tidyr)
library(stringr)
source("R/multiplot.R")
library(hierfstat)
library(ggthemes)
library(gridExtra)
library(viridis)
# load phyologenetic sequence

# load all datasets
seals <- read_excel("data/processed/seal_data_complete.xlsx")
names(seals)[11] <- "IUCN_rating"

# check which variables are non-numeric
non_numeric <- apply(seals, 2, is.numeric)

seal_names_real <- seals$seal_names_real
# idea: make multiple heatmaps


## bottleneck

## ratio
bottleneck_tests <- seals %>% select_("species",  "IAM_ratio", "TPM70_ratio", "TPM90_ratio", "TPM95_ratio", "SMM_ratio")  
bot <- melt(bottleneck_tests, id.vars = "species")
bot$species <- factor(bot$species, levels = bot$species[28:1])

p_bot_ratio <- ggplot(bot, aes(x= variable, y = species, fill = value)) + 
    #facet_grid(.~dataset) + 
    geom_tile(color = "white", size = 0.1) +
    labs(x = "microsat mutation model", y = "") +
    #scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
     #   name = "prop. \nhet-exc")  +
    scale_fill_viridis(option="magma", direction = -1) +
    scale_x_discrete(labels=c("IAM", "TPM70", "TPM90", "TPM95", "SMM")) +
    scale_y_discrete(labels = rev(seal_names_real)) +
    # scale_fill_gradient(name = "p-value", label = comma,  breaks = c(0.05, 0.5)) +
    theme_tufte(base_family="Helvetica") +
    # coord_equal() +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9), 
        axis.text.x = element_text(angle = 50, hjust = 1),
        legend.position="non",
        plot.margin=unit(c(10,0,10,0),"points")) +
    coord_fixed(ratio = 0.7)




gen_div <- seals %>% select_("species", "allelic_richness", "exp_het_mean")
gen_div[2:ncol(gen_div)] <- lapply(gen_div[2:ncol(gen_div)], function(x){
    out <- base::scale(x)
})

# melt
bot <- melt(gen_div, id.vars = "species")
bot$species <- factor(bot$species, levels = bot$species[28:1])

p_diversity <- ggplot(bot, aes(x= variable, y = species, fill = value)) + 
    #facet_grid(.~dataset) + 
    geom_tile(color = "white", size = 0.1) +
    labs(x = "genetic diversity", y = "") +
    #scale_fill_gradientn(colours=c("#081d58", "#253494", "#1d91c0", "#7fcdbb","#c7e9b4", "#ffffd9"),
    #    name = "standardised \nvalue")  +
    scale_fill_viridis(option="magma") +
    scale_x_discrete(labels=c( "allelic\nrichness ", "NEIDI")) +
    # scale_y_discrete(labels = rev(seal_names_real)) +
    scale_y_discrete(labels = element_blank()) + 
    # scale_fill_gradient(name = "p-value", label = comma,  breaks = c(0.05, 0.5)) +
    theme_tufte(base_family="Helvetica") +
   #  coord_equal() +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9), 
        axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text.y = element_blank(),
        #axis.title.y = element_blank(),
        legend.position="non",
        plot.margin=unit(c(10,0,10,0),"points")) +
    coord_fixed(ratio = 0.7)


# multiplot(p_bot_ratio, p_diversity, cols = 2)

# 
# library(gridExtra)
# cclist <- list(p_bot_ratio, p_diversity)
# cclist$ncol <- 2
# do.call(grid.arrange, cclist)




options(scipen = 999)
seals_g2 <- seals
# seals_g2$species <- factor(seals_g2$species, levels = seals_g2$species[28:1])
seals_g2$seal_names_real <- factor(seals_g2$seal_names_real, levels = seals_g2$seal_names_real[28:1])

seals_g2$IUCN_rating[is.na(seals_g2$IUCN_rating)] <- "data deficient"
seals_g2$IUCN_rating <- factor(seals_g2$IUCN_rating, levels = c("least concern", "vulnerable", "near threatened", "endangered",
                                                                "data deficient"))
#levels(seals$IUCN_rating)
g2_plot <- ggplot(seals_g2, aes(x=seal_names_real, y=g2, color = IUCN_rating)) + 
    geom_errorbar(aes(ymin=CIlow, ymax=CIup), width=.01) +
    geom_point(aes(size = Abundance), alpha = 0.5) +
    # scale_size_continuous(breaks=c(1000, 10000, 100000, 1000000, 5000000),
    #     guide = guide_legend(title = "worldwide abundance"),
    #     range = c(1,7)) +
    scale_size_continuous(guide = guide_legend(title = "Worldwide Abundance"),
        range = c(0.1,5), trans = "log", breaks=c(500, 1000, 10000, 100000, 1000000)) +
    scale_color_manual(values = c("#66c2a5", "#fdae61", "#f46d43", "#d53e4f", "#000000"),
        guide = guide_legend(title = "IUCN red list rating")) +
    # scale_color_viridis(discrete=TRUE, option = "C") +
    #geom_text(angle = 0, vjust = 1.4, hjust = -0.5, size = 3) +
    #scale_x_discrete(name = aes(real_names_species)) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0))) +
    #scale_x_discrete(labels = real_names_species) +
    geom_hline(yintercept = 0) +
    labs( y = "identity disequilibrium") +
    # labs(y = "identity disequilibrium") +
    scale_x_discrete(labels = element_blank()) + 
    theme(legend.position = c(0.7, 0.76),
        axis.ticks = element_blank(),
        axis.title.y=element_blank(),
        plot.margin=unit(c(10,0,10,0),"points"))  +
    coord_flip() 
    
      

p <- grid.arrange(p_bot_ratio, p_diversity, g2_plot, widths=c(0.4, 0.15, 0.45))

ggsave(p,
    filename = "full_plot_with_full_datasets_bot.jpg",
    width = 11,
    height = 10, units = "in",
    dpi = 300)





