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
source("martin.R")
library(extrafont)
library(extrafontdb)
# phylogeny
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


## model 1: het excess vs. LH -------------------------------------------------------

# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)

stats_mod_div <- 
    all_stats %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logAbundance = ((logAbundance - mean(logAbundance)) / (2*sd(logAbundance))),
        logSSD = (logSSD - mean(logSSD) / (2*sd(logSSD))),
        logbreed_season = ((logbreed_season - mean(logbreed_season, na.rm = TRUE)) / (2*sd(logbreed_season, na.rm = TRUE))))

stats_mod_div <- stats_mod_div%>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))



run_mod <- function(iter){
    MCMCglmm(num_alleles_mean ~ logAbundance + BreedingType + SSD, # , #+ Abundance BreedingType
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_div,nitt=1100000,burnin=100000,thin=1000)
}


library(purrr)
library(readr)
# check if model is saved
model_file_name <- "gendiv_lh.RData"

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
mod_gendiv_lh <- models[[1]]

# summary
summary(mod_gendiv_lh)

# visually inspecting chain convergence
plot(mod_gendiv_lh $Sol)
plot(mod_gendiv_lh $VCV)
autocorr(mod_gendiv_lh $Sol)
autocorr(mod_gendiv_lh $VCV)

# variation explained by phylogeny
var_phy <- mod_gendiv_lh$VCV[, "tip_label"] / (mod_gendiv_lh$VCV[, "tip_label"] + mod_gendiv_lh$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

# commonality analyses and R2
model_file_name_R2 <- "mod_gendiv_lh_R2.RData"

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_gendiv <- mcmcR2::partR2(mod_gendiv_lh, partvars = c("logAbundance", "BreedingType", "SSD"),
        data = stats_mod_div, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS( R2_gendiv, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_gendiv <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))
R2_gendiv
# out <- mcmcR2::R2mcmc(mod_hetexc)
# out$partR2
R2_gendiv$R2 %>% write_delim("data/processed/models/mod_gendiv_lh_R2.txt")
R2_gendiv$SC %>% write_delim("data/processed/models/mod_gendiv_lh_SC.txt")

# save summary to file
mod_gendiv_lh %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_gendiv_lh$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_gendiv_lh$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_gendiv_lh_beta.txt")






# model for plotting
mod_div <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name))[[1]]

mod_div <- MCMCglmm(num_alleles_mean ~ logAbundance +  BreedingType + SSD, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=all_stats,nitt=110000,burnin=10000,thin=100)
summary(mod_div)

# prediction 
pred_df <- data.frame(num_alleles_mean = 0, 
    logAbundance = rep(seq(from = 5, to = 15.5, by = 0.5), each = 2),
    BreedingType = c("land", "ice"), 
    SSD = mean(stats_mod_div$SSD), 
    tip_label = stats_mod_div$tip_label[1])

mod_preds <- data.frame(predict(mod_div, pred_df, interval = "confidence"))

mod_preds$logAbundance <- rep(seq(from = 5, to = 15.5, by = 0.5), each = 2)
mod_preds$BreedingType <- c("land", "ice")
#names(mod_preds)[3] <- "Abundance"

# calculate mean if land and ice
mod_preds2 <- mod_preds
fit <- sapply(seq(from = 1, to = nrow(mod_preds2), by = 2), function(x) out <- mean(c(mod_preds[x,"fit"], mod_preds[x+1,"fit"])))
mod_preds_one <- data.frame(fit, logAbundance = seq(from = 5, to = 15.5, by = 0.5))

p_div1 <- ggplot(aes(logAbundance, num_alleles_mean), data = all_stats) +
    geom_point(size = 4, alpha = 0.4, col = "black") + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col = "black") +
    geom_line(data = mod_preds_one, aes(y = fit), size = 1, alpha = 0.5, col = "black") +
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
        plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    xlab("Global abundance") +
    ylab("Genetic diversity \n(Alleles per 10 individuals)") 

p_div1
ggsave(p_div1, filename = "figures/SMM/gen_div1.jpg", width = 4.3, height = 3.5)


# per land and ice breeding
pred_df <- data.frame(num_alleles_mean = 0, logAbundance = rep(seq(from = 5, to = 15.5, by = 0.5), each = 2), 
    BreedingType = c("land", "ice"), SSD = mean(all_stats$SSD), tip_label = all_stats$tip_label[1])

mod_preds <- data.frame(predict(mod_div, pred_df, interval = "confidence"))

mod_preds$logAbundance <- rep(seq(from = 5, to = 15.5, by = 0.5), each = 2)
mod_preds$BreedingType <- c("land", "ice")


p_div <- ggplot(aes(logAbundance, num_alleles_mean), data = all_stats) +
    geom_point(size = 4, alpha = 0.8, aes(color = BreedingType)) + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col ="grey") +
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
        plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    xlab("Global abundance") +
    ylab("Genetic diversity \n(Alleles per 10 individuals)") 

p_div
ggsave(p_div, filename = "figures/SMM/gen_div.jpg", width = 4.3, height = 3.5)


p <- ggplot(aes(TPM80_ratio, num_alleles_mean), data = all_stats) +
    geom_point(size = 4, alpha = 0.5, col = "black") + # abc_out
    geom_point(size = 4, alpha = 0.8, shape = 21, col ="grey") +
 
    theme_martin() +
    scale_x_continuous(breaks = c(seq(from = 0.1, to = 1, by = 0.2))) + 
    scale_y_continuous(breaks = c(seq(from = 2, to = 10, by =2))) +
    theme(legend.position=c(0.3, 0.8),
        plot.margin = unit(c(1, 0.2, 0.3, 0.2), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    xlab("Heterozygosity-excess\n(A bottleneck signature)") +
    ylab("Genetic diversity \n(Alleles per 10 individuals)") 

p 
ggsave(p, filename = "figures/SMM/gen_div_hetexc.jpg", width = 4.3, height = 3.5)







