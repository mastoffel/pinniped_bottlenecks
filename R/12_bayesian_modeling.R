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
library(mcmcR2)



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
                              nloc, nind, bot, logAbundance, logbreed_season) %>% 
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








## Model for Hypothesis (3) - Het-excess and demography ---------------

stats_mod_het <- 
    stats_mod %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logAbundance = ((logAbundance - mean(logAbundance)) / (2*sd(logAbundance))),
        logbreed_season = ((logbreed_season - mean(logbreed_season, na.rm = TRUE)) / (2*sd(logbreed_season, na.rm = TRUE))))


stats_mod_het <- stats_mod_het %>% mutate(BreedingType = relevel(BreedingType, ref = "land"))

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








ggplot(aes(TPM70_ratio,Generation_time), data = stats_mod) + geom_point() + geom_smooth(method = "lm")



# investigation of breeding season length vs allelic richness
stats_mod_div %>% 
    #filter(BreedingType == "land") %>% 
    ggplot(aes(num_alleles_mean, logbreed_season)) + 
    geom_point(aes(color = "BreedingType")) + 
    geom_smooth(method = "lm", se = FALSE, aes(color = "BreedingType")) + theme_martin() 
   # geom_label_repel(aes(label = common))



mod_land_df <- stats_mod_div %>% 
    filter(BreedingType == "land")

mod_land <- MCMCglmm(num_alleles_mean ~ logbreed_season, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=mod_land_df,nitt=1100000,burnin=100000,thin=1000)

summary(mod_land)











test <- MCMCglmm(num_alleles_mean ~ SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_div,nitt=110000,burnin=10000,thin=100)
summary(test)

test <- stats_mod_div %>% mutate(logbreed_season = ifelse(is.na(logbreed_season), mean(logbreed_season, na.rm = TRUE), logbreed_season))
stats_mod_div
mod2 <- MCMCglmm(num_alleles_mean ~ logAbundance + BreedingType + BreedingType * logbreed_season, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=test,nitt=1100000,burnin=100000,thin=1000)
summary(mod2)




mod2 <- MCMCglmm(num_alleles_mean ~ logAbundance + BreedingType + logbreed_season + logbreed_season * BreedingType, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_div,nitt=500000,burnin=1000,thin=500)

summary(mod1)
R2mod1 <- mcmcR2(mod1)
R2mod1

mod2 <- lm(num_alleles_mean ~ logAbundance + BreedingType + Generation_time + SSD,  data=stats_mod_div)
summary(mod2)

mod1 <- MCMCglmm(TPM70_ratio ~ logAbundance + BreedingType + Generation_time + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_div,nitt=10000,burnin=1000,thin=50)
summary(mod1)

out <- mcmcR2(mod1)


mod2 <- lm(TPM70_ratio ~ logAbundance + BreedingType + Generation_time + SSD, data = stats_mod_div)
summary(mod2)


mod1 <- MCMCglmm(bot ~ logAbundance + BreedingType + Generation_time + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_div,nitt=1000000,burnin=1000,thin=1000)
summary(mod1)


vmVarF<-numeric(9900)

for(i in 1:9900){
    
    Var<-var(as.vector(mod3$Sol[i,] %*% t(mod3$X)))
    
    vmVarF[i]<-Var
}


R2m<-vmVarF/(vmVarF+mod3$VCV[,1]+mod3$VCV[,2])

mean(R2m)



ggplot(all_stats, aes(x = TPM90_ratio, y = prop_low_afs_mean)) + geom_point(size = 4, alpha = 0.5) + geom_smooth(method = "lm")  #+ scale_x_continuous(trans = "log")
 
# An alternative way for getting the same result
mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))


R2m<-vmVarF/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])



summary(mod1)
summary(mod2)


plot(mod1$Sol)
plot(mod1$VCV)

posterior.mode(mod1$VCV)
posterior.mode(mod1$Sol)
HPDinterval(mod1$VCV)
HPDinterval(mod1$Sol)


# R2
mFixed <- mean(mod1$Sol[,2]) * mod1$X[, 2] + mean(mod1$Sol[, 3]) * mod1$X[, 3] 
mFixed
mVarF<- var(mFixed)

# An alternative way for getting the same result
mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))

# R2GLMM(m) - marginal R2GLMM
# Equ. 26, 29 and 30

# MCMCglmm - marginal
mVarF/(mVarF+sum(apply(mod1$VCV,2,mean)))

# alternative with crebile intervals
vmVarF<-numeric(198)

for(i in 1:198){
    Var<-var(as.vector(mod1$Sol[i,] %*% t(mod1$X)))
    vmVarF[i]<-Var
}

R2m<-vmVarF/(vmVarF+mod1$VCV[,1]+mod1$VCV[,2]) # include here all random effects plus errors
mean(R2m)

posterior.mode(R2m)
HPDinterval(R2m)

# R2GLMM(c) - conditional R2GLMM for full model
# Equ. 30
(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1])/(VarF + VarCorr(mF)$Container[1] + VarCorr(mF)$Population[1] + (attr(VarCorr(mF), "sc")^2))
# MCMCglmm - conditional
(mVarF+sum(apply(mmF$VCV,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV,2,mean)))
# alternative with crebile intervals
R2c<-(vmVarF+mmF$VCV[,1]+mmF$VCV[,2])/(vmVarF+mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3])
mean(R2c)
posterior.mode(R2c)
HPDinterval(R2c) 



autocorr.diag(model_simple$Sol)
plot(model_simple$VCV)
autocorr.diag(model_simple$VCV)
autocorr.plot(model_simple$VCV)
autocorr(model_simple$VCV)
heidel.diag(model_simple$VCV)
geweke.plot(model_simple$VCV)
geweke.plot(model_simple$Sol)

gelman.plot()

summary(model_simple)
plot(model_simple)


data("bird.families")
phylo.effect<-rbv(bird.families, 1, nodes="TIPS") 
phenotype<-phylo.effect+rnorm(dim(phylo.effect)[1], 0, 1)  
test.data<-data.frame(phenotype=phenotype, taxon=row.names(phenotype))


bird.families <- makeNodeLabel(bird.families)
some.families <- c("Certhiidae", "Paridae", "Gruidae","Struthionidae")

Nphylo <- drop.tip(bird.families, setdiff(bird.families$tip.label,some.families))
INphylo <- inverseA(Nphylo)
INphylo$pedigree

Aphylo <- vcv.phylo(Nphylo, cor = T)



## tryout phylogenetic models -----------------------------
test <- INphylo$Ainv
sum(rownames(test) %in% all_stats$tip_label)
data(shorebird)
seal <- comparative.data(tree_final, data.frame(all_stats), "phylo", vcv = TRUE)
mod <- pgls(num_alleles_mean~Abundance + SSD + Generation_time, seal)
summary(mod)
pgls.profile(mod)
mod1 <- lm(num_alleles_mean~ Abundance + SSD + Generation_time + abc_out, data = all_stats_div)
summary(mod1)
crunch_mod <- crunch(num_alleles_mean~Abundance + SSD + Generation_time, seal)
summary(crunch_mod)



ggplot(data = all_stats, aes(x = TPM90_ratio, y = bot)) + geom_point() + geom_smooth(method = "lm")

summary(lm(TPM90_ratio~bot, data = all_stats))

