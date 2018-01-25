# phylogenetic comparative analysis
library(ggtree)
library(ape)
library(phytools)
library(dplyr)
library(readxl)
library(stringr)
library(ggthemr)
library(reshape2)
library(scales)
library(forcats)
library(readr)
# for comparative analysis
library(caper)
library(yhat)
library(dplyr)
library(GGally)
library(reshape2)
library(magrittr)
library(tibble)
library(mcmcR2)
library(ggrepel)
source("R/martin.R")



#load modified phylogeny and all stats
tree_final <- read.tree("data/raw/phylogeny/28_species_10ktrees.tre")
#load modified phylogeny and all stats
#tree_final <- read.tree("data/raw/phylogeny/higdon_mod2_28.tre")
# produce short names for plotting
short <- c("W", "NFS", "SSL", "CSL", "GSL", "SASL", "AFS", "NZSL", "AntFS", "NZFS", "SAFS", "GFS", 
    "BS", "HoS", "GS", "HS", "ARS", "SRS", "BRS", "LRS", "MMS", "HMS", "NES", "SES", "CS", "RS", "LS", "WS")

# all_stats_tree is from 10_visualise_phylogeny.R
all_stats <- read_csv("data/processed/all_stats_tree.csv") %>% 
    mutate(SSD = male_weight/female_weight) %>% 
    mutate(abc_out = ifelse(bot > 0.5, "bot", "neut")) %>% 
    mutate(BreedingType = factor(BreedingType, levels = c("ice", "land", "both"))) %>% 
    mutate(logAbundance = log(Abundance),
        logharem_size = log(harem_size),
        logmale_weight = log(male_weight),
        logbreed_season = log(breeding_season_length),
        loglactation_length = log(lactation_length),
        logSSD = log(SSD)) %>% 
    # order factors according to tree
    mutate(tip_label = fct_inorder(factor(tip_label)),
        species = fct_inorder(factor(species)),
        latin = fct_inorder(factor(latin)),
        common = fct_inorder(factor(common)),
        short = fct_inorder(factor(short)))
# count grey and harbour seal to land breeding
all_stats[all_stats$BreedingType == "both", "BreedingType"] <- "land" 
all_stats <- all_stats %>% mutate(BreedingType = as.factor(as.character(BreedingType))) %>% data.frame()

library(tidyr)

# select variables for modeling -----------------------------------------------------------
stats_mod <- all_stats %>% 
    dplyr::select(common, TPM70_ratio,TPM80_ratio, TPM90_ratio, num_alleles_mean, harem_size, SSD, BreedingType, 
        male_weight, breeding_season_length, lactation_length, life_span_years, 
        Abundance, Generation_time, Longevity, tip_label, mratio_mean,
        obs_het_mean, mean_allele_range, IUCN_rating, prop_low_afs_mean,
        nloc, nind, bot, logAbundance, logbreed_season, logSSD) %>% 
    data.frame()


# phylogenetic mixed model -------------------------------------------------------------------------
# phylogenetic mixed model preparation
library("kinship2")
library(MCMCglmm)

# modeling determinants of genetic diversity
# tree_final$tip.label <- as.character(all_stats$species)
plot(tree_final)

# construct inverse phylo matrix and priors
inv_phylo <- inverseA(tree_final, nodes="TIPS",scale=FALSE)$Ainv #,scale=TRUE
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))


stats_mod_gen <- 
    stats_mod %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        TPM80_ratio  = (TPM80_ratio - mean(TPM80_ratio ) / (2*sd(TPM80_ratio))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD)))) %>% 
    as.data.frame()



# model for plotting
mod_div <- MCMCglmm(num_alleles_mean ~ logAbundance + TPM80_ratio + BreedingType + SSD, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=110000,burnin=10000,thin=100)
summary(mod_div)

stats_mod$logAbundance

# prediction 

pred_df <- data.frame(num_alleles_mean = 0, logAbundance = rep(seq(from = 5, to = 15.5, by = 0.5), each = 2), TPM80_ratio = mean(stats_mod$TPM80_ratio),
    BreedingType = c("land", "ice"), SSD = mean(stats_mod$SSD), tip_label = stats_mod$tip_label[1])

mod_preds <- data.frame(predict(mod_div, pred_df, interval = "confidence"))

mod_preds$logAbundance <- rep(seq(from = 5, to = 15.5, by = 0.5), each = 2)
mod_preds$BreedingType <- c("land", "ice")
#names(mod_preds)[3] <- "Abundance"


# calculate mean if land and ice
mod_preds2 <- mod_preds
fit <- sapply(seq(from = 1, to = nrow(mod_preds2), by = 2), function(x) out <- mean(c(mod_preds[x,"fit"], mod_preds[x+1,"fit"])))
mod_preds_one <- data.frame(fit, logAbundance = seq(from = 5, to = 15.5, by = 0.5))

p_div1 <- ggplot(aes(logAbundance, num_alleles_mean), data = stats_mod) +
    geom_point(size = 4, alpha = 0.4) + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") +
    geom_line(data = mod_preds_one, aes(y = fit), size = 1, alpha = 0.5) +
    #geom_line(stat = "smooth", method = "lm",  alpha = 0.6,  aes(color = BreedingType)) +
    #geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.08, 
    #aes(fill = BreedingType)) +
    #scale_color_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    # scale_fill_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    theme_martin() +
    scale_y_continuous(breaks = c(seq(from = 2, to = 10, by =2))) +
    scale_x_continuous(trans = "log", breaks = c(log(100), log(1000), log(10000), log(100000), log(1000000), log(10000000)), 
        labels = c(expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}),  expression(10^{7})),
        limits = c(4.5, 16.2)) + 
    theme(legend.position=c(0.3, 0.8),
        plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm")) +
    xlab("Global abundance") +
    ylab("Genetic diversity \n(Alleles per 10 individuals)") 

p_div1
ggsave(p_div1, filename = "figures/DZG_talk/gen_div1.jpg", width = 4.3, height = 3.5)


# per land and ice breeding
pred_df <- data.frame(num_alleles_mean = 0, logAbundance = rep(seq(from = 5, to = 15.5, by = 0.5), each = 2), TPM80_ratio = mean(stats_mod$TPM80_ratio),
    BreedingType = c("land", "ice"), SSD = mean(stats_mod$SSD), tip_label = stats_mod$tip_label[1])

mod_preds <- data.frame(predict(mod_div, pred_df, interval = "confidence"))

mod_preds$logAbundance <- rep(seq(from = 5, to = 15.5, by = 0.5), each = 2)
mod_preds$BreedingType <- c("land", "ice")


p_div <- ggplot(aes(logAbundance, num_alleles_mean), data = stats_mod) +
    geom_point(size = 4, alpha = 0.8, aes(color = BreedingType)) + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") +
    geom_line(data = mod_preds, aes(y = fit, color = BreedingType), size = 1, alpha = 0.5) +
    #geom_line(stat = "smooth", method = "lm",  alpha = 0.6,  aes(color = BreedingType)) +
    #geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.08, 
        #aes(fill = BreedingType)) +
    scale_color_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    # scale_fill_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    theme_martin() +
    scale_x_continuous(trans = "log", breaks = c(log(100), log(1000), log(10000), log(100000), log(1000000), log(10000000)), 
        labels = c(expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}),  expression(10^{7})),
        limits = c(4.5, 16.2)) + 
    scale_y_continuous(breaks = c(seq(from = 2, to = 10, by =2))) +
    theme(legend.position=c(0.3, 0.8),
        plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm")) +
    xlab("Global abundance") +
    ylab("Genetic diversity \n(Alleles per 10 individuals)") 

p_div
ggsave(p_div, filename = "figures/DZG_talk/gen_div.jpg", width = 4.3, height = 3.5)




p_div2 <- ggplot(aes(TPM80_ratio, num_alleles_mean), data = stats_mod) +
    geom_point(size = 5, alpha = 0.4) + # abc_out
    geom_point(size = 5, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    scale_y_continuous(breaks = c(seq(from = 2, to = 10, by =2))) +
    scale_x_continuous(breaks= c(seq(from = 0.1, to = 1, by = 0.2))) + 
    theme(legend.position=c(0.3, 0.8),
        plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm")) +
    xlab("Heterozygosity-excess") +
    ylab("Genetic diversity \n(Alleles per 10 individuals)") 
ggsave(p_div2, filename = "figures/DZG_talk/gen_div_hetexc.jpg", width = 4.3, height = 3.5)



library(dplyr)
# model2 for plotting // heterozygosity-excess
# stats_mod_land <- filter(stats_mod, BreedingType == "land")
# 
# mod_het <- MCMCglmm(TPM80_ratio ~ SSD, #
#     random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
#     family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
#     data=stats_mod_land,nitt=110000,burnin=10000,thin=100)
# summary(mod_het)
# 
# 
# pred_df <- data.frame(TPM80_ratio = 0, SSD = seq(from = 0.5, to = 8, by = 0.5), tip_label = stats_mod$tip_label[1])
# 
# mod_preds <- data.frame(predict(mod_het, pred_df, interval = "confidence"))
# mod_preds$SSD <- seq(from = 0.5, to = 8, by = 0.5)
# 
# ggplot(aes(SSD, TPM80_ratio), data = stats_mod_land) +
#     geom_line(data = mod_preds, aes(y = fit), size = 1, alpha = 0.8, color = "#d8b365") +
#     geom_point(size = 4, alpha = 0.4, color = "#d8b365") + # abc_out
#     geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") +
#     # scale_fill_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
#     theme_martin() +
#     theme(plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm")) +
#     xlab("Sexual Weight Dimorphism \n (A proxy for Breeding System)") +
#     ylab("Heterozygosity-excess \n (A bottleneck signature)") 



# full model but without land/ice
mod_het <- MCMCglmm(TPM80_ratio ~ logSSD, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=110000,burnin=10000,thin=100)
summary(mod_het)

pred_df <- data.frame(TPM80_ratio = 0, logSSD = seq(from = -0.5, to = 2.1, by = 0.1), tip_label = stats_mod$tip_label[1])
mod_preds <- data.frame(predict(mod_het, pred_df, interval = "confidence"))

mod_preds$logSSD <- seq(from = -0.5, to = 2.1, by = 0.1)

p_het <- ggplot(aes(logSSD, TPM80_ratio), data = stats_mod) +
    geom_line(data = mod_preds, aes(y = fit), size = 1, alpha = 0.8, color = "grey") +
    geom_point(size = 5, alpha = 0.8, color = "grey") + # abc_out
    geom_point(size = 5, alpha = 0.8, shape = 21, col = "black") +
    # scale_fill_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    theme_martin() +
    # labs(title = "Bottleneck signatures vary with Mating System") +
    theme(plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm"),
        legend.position = c(0.2, 0.85),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text=element_text(size=14),
        legend.title = element_text(size = 14)) +
    scale_x_continuous(breaks = log(c(1,2,3,4,5,6,7)), labels = c(1,2,3,4,5,6,7)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) +
    scale_color_manual("Breeding Habitat", values = c("cornflowerblue", "#d8b365")) +
    xlab("Sexual Weight Dimorphism\n(A proxy for Mating System)") +
    ylab("Heterozygosity-excess") 

p_het
ggsave(p_het,filename= "figures/DZG_talk/het_exc.jpg", width = 5.5, height = 4)
 

p_het_col <- ggplot(aes(logSSD, TPM80_ratio), data = stats_mod) +
    geom_line(data = mod_preds, aes(y = fit), size = 1, alpha = 0.8, color = "grey") +
    geom_point(size = 5, alpha = 0.8, aes(color = BreedingType)) + # abc_out
    geom_point(size = 5, alpha = 0.8, shape = 21, col = "black") +
    # scale_fill_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    theme_martin() +
   # labs(title = "Bottleneck signatures vary with Mating System") +
    theme(plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm"),
          legend.position = c(0.2, 0.85),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.text=element_text(size=14),
        legend.title = element_text(size = 14)) +
    scale_x_continuous(breaks = log(c(1,2,3,4,5,6,7)), labels = c(1,2,3,4,5,6,7)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) +
    scale_color_manual("Breeding Habitat", values = c("cornflowerblue", "#d8b365")) +
    xlab("Sexual Weight Dimorphism\n(A proxy for Mating System)") +
    ylab("Heterozygosity-excess") 

p_het_col
ggsave(p_het_col,filename= "figures/DZG_talk/het_exc_col.jpg", width = 5.5, height = 4)





p_het2 <- ggplot(aes(BreedingType, TPM80_ratio), data = stats_mod) +
    geom_boxplot(alpha = 0.3, col = "black",  size = 0.2, width = 0.4, aes(fill = BreedingType)) + #
    geom_point(size = 4, alpha = 0.8, aes(color = BreedingType)) + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") +
    theme_martin() +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab("Breeding Habitat") +
    ylab("Heterozygosity-excess") +
    guides(fill=FALSE, color = FALSE) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.1,0.3,0.3,0.2), "cm") ,
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"), size = 16),
        axis.title.y=element_text(margin=margin(r=0.5, unit = "cm"), size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

p_het2 
ggsave(p_het2,filename= "figures/DZG_talk/het_exc2.jpg", width = 3.5, height = 4)
