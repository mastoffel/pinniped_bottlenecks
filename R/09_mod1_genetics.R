## This script both models and produces the publication ready plot for hypothesis 1:
# Genetic diversity and bottleneck signatures (het-excess) are correlated.

# packages
library(tibble)
library(dplyr)
library(ggtree)
library(ape)
library(phytools)
library(readxl)
library(stringr)
library(viridis)
library(ggtree)
library(ggthemr)
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
source("R/martin.R")
library(tidyr)
library("kinship2")
library(MCMCglmm)

## what should this script do:

# modeling
modeling <- FALSE
save_models <- FALSE

# plotting
plotting <- TRUE
save_plots <- FALSE


# load data and prepare mixed models

# load (modified) phylogeney. 26 species from 10ktrees plus 3 subspecies of ringed seal
tree_final <- read.tree("data/raw/phylogeny/29_species_10ktrees.tre")

# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_29_modeling.csv"))

# phylogenetic mixed model preparation

# construct inverse phylo matrix and priors
inv_phylo <- inverseA(tree_final, nodes="TIPS",scale=FALSE)$Ainv #,scale=TRUE
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))






## model 1: genetic diversity vs. het excess -------------------------------------------------------

# make sure here that the df isn't a tibble
# (1) check the relationship between Het excess and the other genetic variables
# stats_mod %>%
#     dplyr::select(c("obs_het_mean","TPM90_ratio", "TPM70_ratio", "num_alleles_mean",
#         "mean_allele_range", "prop_low_afs_mean"
#                 )) %>% #"Generation_time", "Abundance", "logAbundance", "SSD", "BreedingType"
#     ggduo()


### modeling variable at the beginning
if(modeling){
    
# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)
stats_mod_gen <- all_stats %>% 
        mutate(num_alleles_mean = ((num_alleles_mean - mean(num_alleles_mean)) / (2*sd(num_alleles_mean))), 
            obs_het_mean = (obs_het_mean - mean(obs_het_mean) / (2*sd(obs_het_mean))),
            prop_low_afs_mean = (prop_low_afs_mean - mean(prop_low_afs_mean) / (2*sd(prop_low_afs_mean)))) %>% 
        as.data.frame()

# model specification
run_mod <- function(iter){
    MCMCglmm(TPM80_ratio ~ num_alleles_mean + obs_het_mean + prop_low_afs_mean, #
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_gen,nitt=1100000,burnin=100000,thin=1000)
}

library(purrr)
library(readr)
# check if model is saved
model_file_name <- "gen_hetexc_ar_het_afs.RData"

if (!file.exists(paste0("output/mcmcmodels/", model_file_name))){
    # run models
    set.seed(1234)
    models <- purrr::map(1:3, run_mod)
    saveRDS(models, file = paste0("output/mcmcmodels/", model_file_name))
}

models <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name))

# model checks according to holgers paper
plot(mcmc.list(models[[1]]$Sol, models[[2]]$Sol, models[[3]]$Sol))
gelman.diag(mcmc.list(models[[1]]$Sol, models[[2]]$Sol, models[[3]]$Sol))

# one model
mod_gen <- models[[1]]

# summary
summary(mod_gen)

# visually inspecting chain convergence
plot(mod_gen$Sol)
plot(mod_gen$VCV)
autocorr(mod_gen$Sol)
autocorr(mod_gen$VCV)

# variation explained by phylogeny
var_phy <- mod_gen$VCV[, "tip_label"] / (mod_gen$VCV[, "tip_label"] + mod_gen$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

# R2s

# check if model R2s are saved
model_file_name_R2 <- "gen_hetexc_ar_het_afs_R2.RData"

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_gen <- mcmcR2::partR2(mod_gen, partvars = c("num_alleles_mean", "obs_het_mean", "prop_low_afs_mean"),
        data = stats_mod_gen, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS(R2_gen, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_gen <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))

if (save_models){
    R2_gen$R2 %>% write_delim("data/processed/models/mod_gen_R2.txt")
    R2_gen$SC %>% write_delim("data/processed/models/mod_gen_SC.txt")
    # structure coefficients
    
    mod_sum <- summary(mod_gen)
    # save summary to file
    sum_mod_gen <- mod_gen %>% 
        summary() %$%
        solutions %>% 
        as.data.frame() %>% 
        rownames_to_column("components") %>% 
        mutate(post_median = apply(mod_gen$Sol, 2, median)) %>% 
        mutate(post_mode = posterior.mode(mod_gen$Sol)) %>% 
        .[c(1,2,7,8,3:6)] %>% 
        rename(post_mean= post.mean,
            lower =  "l-95% CI",
            upper = "u-95% CI") %>% 
        write_delim("data/processed/models/mod_gen_beta.txt")
}

}








# plots for paper ----------------------------------------------------------------------------------

# plot (1) het-exc vs genetic diversity

# "num_alleles_mean", "obs_het_mean", "prop_low_afs_mean"

# some controls
point_size <- 3.5
point_alpha <- 0.3

# (I) Heterozygosity-excess ------------------------------------------------------------------------
# (1)vs. allelic richness
mod_plot_AR <- MCMCglmm(num_alleles_mean ~ TPM80_ratio , #num_alleles_mean  + 
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=all_stats,nitt=110000,burnin=10000,thin=100)

mean(all_stats$TPM80_ratio)
pred_df_AR <- data.frame(num_alleles_mean =0,
     TPM80_ratio = seq(from = 0.2, to = 1, by = 0.05), 
    tip_label = all_stats$tip_label[1])

mod_preds_AR <- data.frame(predict(mod_plot_AR, pred_df_AR, interval = "confidence")) %>% 
             mutate(TPM80_ratio = seq(from = 0.2, to = 1, by = 0.05))


p1 <- ggplot(aes(TPM80_ratio, num_alleles_mean), data = all_stats) +
    geom_line(data = mod_preds_AR, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    xlab("Heterozygosity-excess") +
    ylab("Allelic richness") +
    theme_martin() +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
        )
p1


# (2) vs. heterozygosity

mod_plot_het <- MCMCglmm(obs_het_mean ~ TPM80_ratio , #num_alleles_mean  + 
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=all_stats,nitt=110000,burnin=10000,thin=100)

pred_df_het <- data.frame(obs_het_mean =0,
    TPM80_ratio = seq(from = 0.2, to = 1, by = 0.05), 
    tip_label = all_stats$tip_label[1])

mod_preds_het <- data.frame(predict(mod_plot_het, pred_df_het, interval = "confidence")) %>% 
    mutate(TPM80_ratio = seq(from = 0.2, to = 1, by = 0.05))


p2 <- ggplot(aes(TPM80_ratio, obs_het_mean), data = all_stats) +
    geom_line(data = mod_preds_het, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    xlab("Heterozygosity-excess") +
    ylab("Observed heterozygosity") +
    theme_martin() +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    )
p2


# (3) vs. low frequ alleles

mod_plot_afs <- MCMCglmm(prop_low_afs_mean ~ TPM80_ratio , #num_alleles_mean  + 
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=all_stats,nitt=110000,burnin=10000,thin=100)

pred_df_afs <- data.frame(prop_low_afs_mean =0,
    TPM80_ratio = seq(from = 0.2, to = 1, by = 0.05), 
    tip_label = all_stats$tip_label[1])

mod_preds_afs <- data.frame(predict(mod_plot_afs, pred_df_afs, interval = "confidence")) %>% 
    mutate(TPM80_ratio = seq(from = 0.2, to = 1, by = 0.05))


p3 <- ggplot(aes(TPM80_ratio, prop_low_afs_mean), data = all_stats) +
    geom_line(data = mod_preds_afs, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    xlab("Heterozygosity-excess") +
    ylab("% low frequency alleles") +
    theme_martin() +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    )
p3





# ABC bottleneck probability -----------------------------------------------------------------------

mod_plot_AR_abc <- MCMCglmm(num_alleles_mean ~ bot, #num_alleles_mean  + 
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=all_stats,nitt=110000,burnin=10000,thin=100)


pred_df_AR_abc <- data.frame(num_alleles_mean =0,
    bot = seq(from = 0, to = 1, by = 0.05), 
    tip_label = all_stats$tip_label[1])

mod_preds_AR_abc <- data.frame(predict(mod_plot_AR_abc, pred_df_AR_abc, interval = "confidence")) %>% 
    mutate(bot= seq(from = 0, to = 1, by = 0.05))


p4 <- ggplot(aes(bot, num_alleles_mean), data = all_stats) +
    geom_line(data = mod_preds_AR_abc, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    xlab("Bottleneck model probability") +
    ylab("Allelic richness") +
    theme_martin() +
    theme(#panel.grid.major = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    )
p4









# plots for Marine Mammal Conference

p1_smm <- ggplot(aes(TPM80_ratio, num_alleles_mean), data = all_stats) +
    #geom_line(data = mod_preds_AR, aes(y = fit), size = 1, alpha = 0.3) +
    geom_point(size = point_size, alpha = 0.3) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    xlab("Heterozygosity excess") +
    ylab("Allelic richness") +
    theme_martin(grid = TRUE) +
    #ggrepel::geom_text_repel(data = subset(all_stats, num_alleles_mean < 4),
    #        aes(label = common)) +
    #ggrepel::geom_text_repel(data = subset(all_stats, num_alleles_mean > 7),
    #    aes(label = common)) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1)) +
    theme(#panel.grid.major = element_blank(),
        #plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm") ,
        axis.line.x = element_line(color="#cccccc", size=0.3),
        axis.ticks.x = element_line(color="#cccccc", size=0.3),
        axis.text.x = element_text(margin = margin(t = 5))
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    )
p1_smm

ggplot2::ggsave(filename = "figures/SMM/p1_smm.jpg", width = 4, height = 3)




