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
#load modified phylogeny and all stats
tree_final <- read.tree("data/raw/phylogeny/higdon_mod2_28.tre")
all_stats <- read_csv("data/processed/all_stats_tree.csv") %>% 
                mutate(SSD = male_weight/female_weight) %>% 
                mutate(abc_out = ifelse(bot > 0.5, "bot", "neut")) %>% 
                mutate(BreedingType = factor(BreedingType, levels = c("ice", "land", "both")))

# produce short names for plotting
short <- c("W", "NFS", "SSL", "CSL", "GSL", "SASL", "AFS", "NZSL", "AntFS", "NZFS", "SAFS", "GFS", 
    "BS", "HoS", "GS", "HS", "ARS", "SRS", "BRS", "LRS", "MMS", "HMS", "NES", "SES", "CS", "RS", "LS", "WS")
all_stats$short <- short

# order factors according to tree
all_stats <- all_stats %>% 
    mutate(tip_label = fct_inorder(factor(tip_label)),
        species = fct_inorder(factor(species)),
        latin = fct_inorder(factor(latin)),
        common = fct_inorder(factor(common)),
        short = fct_inorder(factor(short)))

# WriteXLS(all_stats, "seal_stats_and_LH.xls")

# cdat <- comparative.data(tree_final, as.data.frame(all_stats), names.col = "tip_label")

# pairs plot for all important variables -----------------------------------------------------------
stats_mod <- all_stats %>% 
                dplyr::select(TPM70_ratio, num_alleles_mean, harem_size, SSD, BreedingType, 
                              male_weight, breeding_season_length, lactation_length, life_span_years, 
                              Abundance, Generation_time, Longevity, tip_label, mratio_mean,
                              obs_het_mean, mean_allele_range, IUCN_rating, prop_low_afs_mean,
                              nloc, nind) %>% 
                              mutate(logAbundance = log(Abundance),
                                     harem_size = log(harem_size),
                                     male_weight = log(male_weight)) %>% 
                              data.frame()



# phylogenetic mixed model -------------------------------
# phylogenetic mixed model
library("kinship2")
library(MCMCglmm)

# count GS and HS as land breeding for this
# both <- which(stats_mod$short == "GS" | stats_mod$short =="HS")
stats_mod[stats_mod$BreedingType == "both", "BreedingType"] <- "land" 
stats_mod <- stats_mod %>% mutate(BreedingType = as.factor(as.character(BreedingType))) %>% data.frame()

# modeling determinants of genetic diversity
# tree_final$tip.label <- as.character(all_stats$species)
plot(tree_final)

inv_phylo <- inverseA(tree_final, nodes="TIPS")$Ainv #,scale=TRUE
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

# prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=diag(2),nu=0.002))
# make sure that data isnt a tibble



# (1) check the relationship between Het excess and the other genetic variables
stats_mod %>% 
    dplyr::select(c("obs_het_mean", "TPM70_ratio", "num_alleles_mean", "mean_allele_range", "prop_low_afs_mean"
                )) %>% #"Generation_time", "Abundance", "logAbundance", "SSD", "BreedingType"
    ggduo()


mod_bays <- MCMCglmm(TPM70_ratio ~ num_alleles_mean + mean_allele_range + prop_low_afs_mean + obs_het_mean, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod ,nitt=100000,burnin=1000,thin=100)
summary(mod_bays)

autocorr(mod_bays$Sol)
autocorr(mod_bays$VCV)

out <- partR2(mod_bays, partvars = partvars <- c("num_alleles_mean", "obs_het_mean", "mean_allele_range", "prop_low_afs_mean"))

# leave that for now

# check the relationship between Het excess and the other genetic variables
mod0 <- MCMCglmm(TPM70_ratio ~ num_alleles_mean,
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)
summary(mod0)

mod1 <- MCMCglmm(TPM70_ratio ~ obs_het_mean,
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)
summary(mod1)

mod2 <- MCMCglmm(TPM70_ratio ~ mean_allele_range,
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)
summary(mod2)

mod3 <- MCMCglmm(TPM70_ratio ~ prop_low_afs_mean,
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)
summary(mod3)

# seems like just prop_low_afs_mean is correlated to TPM70, but fairly strong


mod1 <- MCMCglmm(num_alleles_mean ~ Abundance + BreedingType + Generation_time + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)

R2mod1 <- mcmcR2(mod1)

mod2 <- MCMCglmm(num_alleles_mean ~ Abundance + Generation_time, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)
R2mod2 <- mcmcR2(mod2)
summary(mod2)
R2mod2

R2mod1[1] - R2mod2[1]

mod3 <- MCMCglmm(num_alleles_mean ~ Generation_time, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)
R2mod3 <- mcmcR2(mod3)
summary(mod3)
R2mod3

mod4 <- MCMCglmm(num_alleles_mean ~  BreedingType + Generation_time + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)

R2mod4 <- mcmcR2(mod4)
R2mod4

stats_mod %>% dplyr::select(c("logAbundance", "BreedingType", "Generation_time", "SSD")) %>% 
        ggduo()


mod1 <- MCMCglmm(num_alleles_mean ~ logAbundance + BreedingType + Generation_time + SSD, # , #+ Abundance BreedingType
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)

summary(mod1)
out <- partR2(mod1, partvars = c( "BreedingType", "logAbundance", "Generation_time", "SSD"))#"Generation_time", "Abundance",
out

mod2 <- MCMCglmm(TPM70_ratio ~ BreedingType + Generation_time + Abundance + SSD,
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)

summary(mod2)

mod3 <- MCMCglmm(num_alleles_mean ~ TPM70_ratio,
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod,nitt=100000,burnin=1000,thin=10)

summary(mod3)


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

