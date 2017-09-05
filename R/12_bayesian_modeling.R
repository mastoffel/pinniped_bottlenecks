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
#source("martin.R")



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
# pairs plot
all_stats %>% # interesting: breed_season, logAbundance, Longevity, 
    dplyr::select(TPM70_ratio, TPM90_ratio, num_alleles_mean, logbreed_season, BreedingType, Abundance, life_span_years,
                 latitude, Longevity, Age_primiparity, Generation_time, logAbundance) %>% 
    gather(-TPM70_ratio, -BreedingType, key = "var", value = "value") %>%  
    ggplot(aes(x = TPM70_ratio, y = value, color = BreedingType)) + geom_point()+ geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ var, scales = "free") + theme_bw()


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



## model 1: genetic diversity vs. het excess -------------------------------------------------------

# make sure here that the df isn't a tibble
# (1) check the relationship between Het excess and the other genetic variables
# stats_mod %>% 
#     dplyr::select(c("obs_het_mean","TPM90_ratio", "TPM70_ratio", "num_alleles_mean", 
#         "mean_allele_range", "prop_low_afs_mean"
#                 )) %>% #"Generation_time", "Abundance", "logAbundance", "SSD", "BreedingType"
#     ggduo()

# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)
stats_mod_gen <- 
    stats_mod %>% 
    mutate(num_alleles_mean = ((num_alleles_mean - mean(num_alleles_mean)) / (2*sd(num_alleles_mean))), 
        obs_het_mean = (obs_het_mean - mean(obs_het_mean) / (2*sd(obs_het_mean))),
        prop_low_afs_mean = (prop_low_afs_mean - mean(prop_low_afs_mean) / (2*sd(prop_low_afs_mean)))) %>% 
    as.data.frame()


# model 1, check convergence
mod_gen <- MCMCglmm(TPM80_ratio ~ num_alleles_mean +  obs_het_mean + prop_low_afs_mean, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_gen,nitt=1100000,burnin=100000,thin=1000)
mod_gen2 <- MCMCglmm(TPM80_ratio ~ num_alleles_mean +  obs_het_mean + prop_low_afs_mean, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_gen,nitt=1100000,burnin=100000,thin=1000)
mod_gen3 <- MCMCglmm(TPM80_ratio ~ num_alleles_mean +  obs_het_mean + prop_low_afs_mean, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_gen,nitt=1100000,burnin=100000,thin=1000)

# diagnostics chain convergence
plot(mcmc.list(mod_gen$Sol, mod_gen2$Sol, mod_gen3$Sol))
# gelman rubin criterion
gelman.diag(mcmc.list(mod_gen$Sol, mod_gen2$Sol, mod_gen3$Sol))
# summary
summary(mod_gen)
# 
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
R2_gen <- mcmcR2::partR2(mod_gen, partvars = c("num_alleles_mean", "obs_het_mean", "prop_low_afs_mean"),
                 data = stats_mod_gen, inv_phylo = inv_phylo, prior = prior, 
                 nitt = 1100000, burnin = 100000, thin = 1000)

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






## Model for Hypothesis (2) - Genetic diversity and demography -------------------------------------

# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)

stats_mod_div <- 
    stats_mod %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logAbundance = ((logAbundance - mean(logAbundance)) / (2*sd(logAbundance))),
        logSSD = (logSSD - mean(logSSD) / (2*sd(logSSD))),
        logbreed_season = ((logbreed_season - mean(logbreed_season, na.rm = TRUE)) / (2*sd(logbreed_season, na.rm = TRUE))))

# plot

stats_mod_div <- stats_mod_div %>% mutate(BreedingType = relevel(BreedingType, ref = "land"))

mod1 <- MCMCglmm(num_alleles_mean ~ logAbundance + BreedingType + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_div,nitt=1100000,burnin=100000,thin=1000)

mod2 <- MCMCglmm(num_alleles_mean ~ logAbundance + BreedingType + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_div,nitt=1100000,burnin=100000,thin=1000)

mod3 <- MCMCglmm(num_alleles_mean ~ logAbundance + BreedingType + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_div,nitt=1100000,burnin=100000,thin=1000)

# diagnostics chain convergence
plot(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))
# gelman rubin criterion
gelman.diag(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))
# summary
summary(mod1)
# 
# model checks
plot(mod1$Sol)
plot(mod1$VCV)
autocorr(mod1$Sol)
autocorr(mod1$VCV)

# summary(lm(formula = num_alleles_mean~logAbundance + BreedingType + SSD, data = stats_mod_div))

# R2s
R2_div <- mcmcR2::partR2(mod1, partvars = c("logAbundance", "BreedingType", "SSD"),
    data = stats_mod_div, inv_phylo = inv_phylo, prior = prior, 
    nitt = 1100000, burnin = 100000, thin = 1000)

# out <- R2mcmc(mod1)
# out$partR2
R2_div$R2 %>% write_delim("data/processed/models/mod_div_R2.txt")
R2_div$SC %>% write_delim("data/processed/models/mod_div_SC.txt")

# save summary to file
sum_mod_gen <- mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply(mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_div_beta.txt")




## Model for Hypothesis (3) - Het-excess and demography --------------------------------------------

stats_mod_het <- 
    stats_mod %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logAbundance = ((logAbundance - mean(logAbundance)) / (2*sd(logAbundance))),
        logbreed_season = ((logbreed_season - mean(logbreed_season, na.rm = TRUE)) / (2*sd(logbreed_season, na.rm = TRUE))),
        logSSD = (logSSD - mean(logSSD) / (2*sd(logSSD))),
        logharem_size = (log(harem_size))) %>% 
    mutate(logharem_size = (logharem_size - mean(logharem_size) / (2*sd(logharem_size))))

stats_mod_het <- stats_mod_het %>% mutate(BreedingType = relevel(BreedingType, ref = "land"))

# modeling with SSD
mod1 <- MCMCglmm(TPM80_ratio ~ SSD + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_het,nitt=1100000,burnin=100000,thin=1000)
mod2 <- MCMCglmm(TPM80_ratio ~ SSD + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_het,nitt=1100000,burnin=100000,thin=1000)
mod3 <- MCMCglmm(TPM80_ratio ~ SSD + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_het,nitt=1100000,burnin=100000,thin=1000)

# diagnostics chain convergence
plot(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))
# gelman rubin criterion
gelman.diag(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))
# summary
summary(mod1)
# 
# model checks
plot(mod1$Sol)
plot(mod1$VCV)
autocorr(mod1$Sol)
autocorr(mod1$VCV)

summary(mod1)

R2_het <- mcmcR2::partR2(mod1, partvars = c("SSD", "BreedingType"),
    data = stats_mod_het, inv_phylo = inv_phylo, prior = prior, 
    nitt = 1100000, burnin = 100000, thin = 1000)

R2_het <- mcmcR2::partR2(mod1, partvars = c("mating_system", "BreedingType"),
    data = stats_mod_het, inv_phylo = inv_phylo, prior = prior, 
    nitt = 1100000, burnin = 100000, thin = 1000)

R2_het

out <- R2mcmc(mod1)
# out$partR2
R2_het$R2 %>% write_delim("data/processed/models/mod_het_R2.txt")
R2_het$SC %>% write_delim("data/processed/models/mod_het_SC.txt")

# save summary to file
sum_mod_gen <- mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply(mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_het_beta.txt")



### supplementary, modeling within land breeding seals ---------------------------------------------
mod_land_df <- stats_mod_div %>% 
    filter(BreedingType == "land")

mod1 <- MCMCglmm(TPM80_ratio ~ SSD, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=mod_land_df ,nitt=1100000,burnin=100000,thin=1000)
mod2 <- MCMCglmm(TPM80_ratio ~ SSD, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=mod_land_df ,nitt=1100000,burnin=100000,thin=1000)
mod3 <- MCMCglmm(TPM80_ratio ~ SSD, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=mod_land_df ,nitt=1100000,burnin=100000,thin=1000)

# diagnostics chain convergence
plot(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))
# gelman rubin criterion
gelman.diag(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))

plot(mod1$Sol)
plot(mod1$VCV)
autocorr(mod1$Sol)
autocorr(mod1$VCV)

summary(mod1)
out <- R2mcmc(mod1)
out$partR2

var_phy <- mod1$VCV[, "tip_label"] / (mod1$VCV[, "tip_label"] + mod1$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

R2_land_het <- mcmcR2::partR2(mod1, partvars = c("SSD"),
    data = mod_land_df, inv_phylo = inv_phylo, prior = prior, 
    nitt = 1100000, burnin = 100000, thin = 1000)

# out$partR2
R2_land_het$R2 %>% write_delim("data/processed/models/mod_het_land_R2.txt")
R2_land_het$SC %>% write_delim("data/processed/models/mod_het_land_SC.txt")

# save summary to file
sum_mod_het_land <- mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply(mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_het_land_beta.txt")



# modeling with harem_size -------------------------------------------------------------------------

stats_mod_het %>% filter(BreedingType == "land") %>% 
    ggplot(aes(logharem_size, TPM80_ratio)) + geom_point() + geom_smooth(method = "lm")

mod1 <- MCMCglmm(TPM80_ratio ~ logharem_size + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_het,nitt=1100000,burnin=100000,thin=1000)
summary(mod1)

summary(mod1)
summary(mod2)
out1 <- R2mcmc(mod1)
out2 <- R2mcmc(mod2)
out1$partR2
out2$partR2
part_mod1 <- partR2(mod1, partvars = c("logharem_size", "BreedingType"), data =stats_mod_het, inv_phylo = inv_phylo, prior = prior, 
    nitt = 110000, burnin = 10000, thin = 100)
part_mod2 <- partR2(mod2, partvars = c("SSD", "BreedingType"), data =stats_mod_het, inv_phylo = inv_phylo, prior = prior, 
    nitt = 110000, burnin = 10000, thin = 100)
part_mod1$R2
part_mod2$R2




ggplot(aes(TPM70_ratio,Generation_time), data = stats_mod) + geom_point() + geom_smooth(method = "lm")



# supplmentary: het-exc vs. bot --------------------------------------------------------------------

# plot
stats_mod <- stats_mod %>% mutate(BreedingType = relevel(BreedingType, ref = "land"))

mod1 <- MCMCglmm(TPM80_ratio ~ bot, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=1100000,burnin=100000,thin=1000)
mod2 <- MCMCglmm(TPM80_ratio ~ bot, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=1100000,burnin=100000,thin=1000)
mod3 <- MCMCglmm(TPM80_ratio ~ bot, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=1100000,burnin=100000,thin=1000)
summary(mod1)


# diagnostics chain convergence
plot(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))
# gelman rubin criterion
gelman.diag(mcmc.list(mod1$Sol, mod2$Sol, mod3$Sol))

plot(mod1$Sol)
plot(mod1$VCV)
autocorr(mod1$Sol)
autocorr(mod1$VCV)

summary(mod1)
out <- R2mcmc(mod1)
out$partR2


R2_het_abc <- mcmcR2::partR2(mod1, partvars = c("bot"),
    data = stats_mod, inv_phylo = inv_phylo, prior = prior, 
    nitt = 1100000, burnin = 100000, thin = 1000)

# out$partR2
R2_het_abc$R2 %>% write_delim("data/processed/models/mod_het_abc_R2.txt")
R2_het_abc$SC %>% write_delim("data/processed/models/mod_het_abc_SC.txt")

# save summary to file
R2_het_abc <- mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply(mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_het_abc_beta.txt")



# supplmentary: gen_div vs. IUCN -------------------------------------------------------------------

# plot
stats_mod_IUCN <- stats_mod %>% 
                    mutate(IUCN_binary = case_when(IUCN_rating == "vulnerable" ~ "concerned",
                                                   IUCN_rating == "near threatened" ~ "concerned",
                                                   IUCN_rating == "endangered" ~ "concerned",
                                                   IUCN_rating == "least concern" ~ "least concern"))

ggplot(stats_mod_IUCN, aes(IUCN_binary, num_alleles_mean)) + geom_boxplot() + geom_point(size = 3)

mod1 <- MCMCglmm(num_alleles_mean ~ IUCN_binary, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_IUCN,nitt=1100000,burnin=100000,thin=1000)

summary(mod1)

# save summary to file
mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply(mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_IUCN_vs_AR_beta.txt")


# supplmentary: Het-excs vs. IUCN ------------------------------------------------------------------

mod1 <- MCMCglmm(TPM80_ratio ~ IUCN_binary, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_IUCN,nitt=1100000,burnin=100000,thin=1000)

summary(mod1)

# save summary to file
mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply(mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_IUCN_vs_hetexc_beta.txt")



# supplmentary: bot vs. IUCN -----------------------------------------------------------------------

mod1 <- MCMCglmm(bot ~ IUCN_binary, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_IUCN,nitt=1100000,burnin=100000,thin=1000)

summary(mod1)

# save summary to file
mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply(mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim("data/processed/models/mod_IUCN_vs_ABCprob_beta.txt")


