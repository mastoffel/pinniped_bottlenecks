# Analysing relationship between AR and the two new life-history traits 

# phylogenetic comparative analysis
# packages
library(pacman)
p_load(tibble, dplyr, ggtree, ape, phytools, readxl, stringr, viridis, ggtree, ggthemr, cowplot,
    ggthemes, ggimage, RColorBrewer, scales, forcats, readr, caper, yhat, dplyr, ggrepel,
    GGally, tidyr, kinship2, MCMCglmm, purrr, readr)
source("R/martin.R")
library(mcmcR2)
## what should this script do:

# load data and prepare mixed models

# load (modified) phylogeney. 26 species from 10ktrees plus 3 subspecies of ringed seal
tree_final <- read.tree("data/raw/phylogeny/30_species_10ktrees_final.tre")

# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_30_modeling.csv"))

# phylogenetic mixed model preparation

# construct inverse phylo matrix and priors
inv_phylo <- inverseA(tree_final, nodes="TIPS",scale=FALSE)$Ainv #,scale=TRUE
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

nitt <- 110000
burnin <- 10000
thin <- 100


ggplot(all_stats, aes(num_alleles_mean, Generation_time)) +
    geom_point() +
    geom_smooth(method = "lm")

ggplot(all_stats, aes(num_alleles_mean, breeding_season_length)) +
    geom_point() +
    geom_smooth(method = "lm")

mod1 <- MCMCglmm(num_alleles_mean ~ Generation_time + breeding_season_length, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=all_stats,nitt=nitt,burnin=burnin,thin=thin)

sum_mod1 <- mod1 %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    rownames_to_column("components") %>% 
    mutate(post_median = apply( mod1$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode( mod1$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") 

R2_mod1 <- R2mcmc(mod1)$partR2
R2_mod1




# supplementary plot

p <- ggplot(all_stats, aes(BreedingType, Generation_time)) +
    geom_boxplot(alpha = 0.5, col = "darkgrey",  size = 0.5, width = 0.7,  outlier.shape = NA) + #aes(fill = BreedingType),
    geom_jitter(size = 3, alpha = 0.6, shape = 21, col = "black", aes(fill = BreedingType), width = 0.2) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
  #  scale_color_manual(values = c("#d8b365", "cornflowerblue")) +
    scale_fill_manual("Breeding habitat", values = c( "cornflowerblue", "#d8b365")) +
    xlab("Breeding Habitat") +
    ylab("Generation time") +
    guides(fill=FALSE)

ggsave(p, filename = "other_stuff/figures/figures_final/Sup_BreedHab_GenTime.jpg", width = 3, height = 2.5)    
    



