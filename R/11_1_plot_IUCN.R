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
source("R/martin.R")

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


# select variables for modeling -----------------------------------------------------------
stats_mod <- all_stats %>% 
    dplyr::select(common, TPM70_ratio,TPM80_ratio, TPM90_ratio, num_alleles_mean, harem_size, SSD, BreedingType, 
        male_weight, breeding_season_length, lactation_length, life_span_years, 
        Abundance, Generation_time, Longevity, tip_label, mratio_mean,
        obs_het_mean, mean_allele_range, IUCN_rating, bot) %>% 
    data.frame()

# plot
stats_mod_IUCN <- stats_mod %>% 
    mutate(IUCN_binary = case_when(IUCN_rating == "vulnerable" ~ "concern",
        IUCN_rating == "near threatened" ~ "concern",
        IUCN_rating == "endangered" ~ "concern",
        IUCN_rating == "least concern" ~ "least concern"))
## modeling and plotting 

p1 <-  ggplot(data = stats_mod_IUCN, aes(IUCN_binary, num_alleles_mean)) +
    geom_boxplot(alpha = 0.3, col = "black",  size = 0.2, width = 0.4) + #
    geom_point(size = 3.5, alpha = 0.6, aes(color = BreedingType)) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab(" ") +
    ylab("Allelic richness") +
    # guides(fill=FALSE, color = FALSE) +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(2,0.3,0.3,0.2), "cm") ,
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm")),
        legend.position = "none"

        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    ) 
p1

p2 <-  ggplot(data = stats_mod_IUCN, aes(IUCN_binary, TPM80_ratio)) +
    geom_boxplot(alpha = 0.3, col = "black",  size = 0.2, width = 0.4) + #
    geom_point(size = 3.5, alpha = 0.6, aes(color = BreedingType)) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab("IUCN status") +
    ylab("Heterozygosity-excess") +
    # guides(fill=FALSE, color = FALSE) +
    ylim(c(0, 1)) +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5,0.3,0.3,0.2), "cm") ,
        axis.title.x= element_text(margin=margin(t=0.5, unit = "cm")),
        legend.position="top"
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    ) 
p2

p3 <-  ggplot(data = stats_mod_IUCN, aes(IUCN_binary, bot)) +
    geom_boxplot(alpha = 0.3, col = "black",  size = 0.2, width = 0.4) + #
    geom_point(size = 3.5, alpha = 0.6, aes(color = BreedingType)) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab(" ") +
    ylab("ABC bottleneck model \nprobability") +
    # guides(fill=FALSE, color = FALSE) +
    ylim(c(0, 1)) +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(2,0.3,0.3,0.2), "cm") ,
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm")),
        legend.position="none"
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    ) 
p3


p_final <- plot_grid(p1, p2, p3, cols = 3, labels = c("A", "B", "C"))
p_final
ggsave(p_final, filename = "figures/figures_final/IUCN.jpg", width = 10, height = 3.5)
