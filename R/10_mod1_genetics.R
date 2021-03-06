## This script both models and produces the publication ready plot for hypothesis 1:
# Genetic diversity and bottleneck signatures (het-excess, ABC pbot) are correlated.

# files needed
# (1) all_stats_30_modeling.csv

# packages
library(pacman)
p_load(tibble, dplyr, ggtree, ape, phytools, readxl, stringr, viridis, ggtree, ggthemr, cowplot,
       ggthemes, ggimage, RColorBrewer, scales, forcats, readr, caper, yhat, dplyr, ggrepel,
       GGally, tidyr, kinship2, MCMCglmm, purrr, readr)
source("R/martin.R")
## what should this script do:

# modeling
modeling <- TRUE

# plotting
plotting <- TRUE


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

# Models for for figure "do bottleneck signatures predict genetic diversity?"
## model 1: ABCprob(bot) ~ genetic diversity -------------------------------------------------------

### modeling variable at the beginning
if(modeling){
    
# standardize variables (by 1 sd, as no binary predictor is involved)
stats_mod_gen <- all_stats %>% 
            mutate(bot_stand = as.numeric(scale(bot))) %>% 
            mutate(TPM80_ratio_stand = as.numeric(scale(TPM80_ratio)))


## Genetic model (1): Allelic richness ~ ABC bottleneck probability
# model specification
run_mod <- function(iter){
    MCMCglmm(num_alleles_mean ~ bot_stand, #
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_gen,nitt=nitt,burnin=burnin,thin=thin)
}


# model name
mod_name <- "gen1_bot"

# check if model is saved
model_file_name <- paste0(mod_name, ".RData")

# run 3 independent chains for model checking
if (!file.exists(paste0("output/mcmcmodels/", model_file_name))){
    # run models
    set.seed(1234)
    models <- purrr::map(1:3, run_mod)
    saveRDS(models, file = paste0("output/mcmcmodels/", model_file_name))
}
# load models
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
library(mcmcR2)
set.seed(312)
for (r2type in c("marginal", "conditional")) {
#for (r2type in c("conditional")) {    
    model_file_name_R2 <- paste0(mod_name, "_R2_", r2type)
    
    if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){

        R2_gen <- mcmcR2::partR2(mod_gen, partvars = c("bot_stand"), type = r2type,
            data = stats_mod_gen, inv_phylo = inv_phylo, prior = prior, 
            nitt = nitt, burnin = burnin, thin = thin)
        saveRDS(R2_gen, file = paste0("output/mcmcmodels/", model_file_name_R2, ".RData"))
    }
    
    R2_gen$R2 %>% write_delim(paste0("output/mcmcmodels/", model_file_name_R2 ,".txt"))

}
R2_gen$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))
# summary
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
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))

}

## model 2: het-exc(TPM80_ratio) ~ genetic diversity ----------------------------------------------

### modeling variable at the beginning
if(modeling){
    
    # standardize variables (by 1 sd, as no binary predictor is involved)
    stats_mod_gen <- all_stats %>% 
        mutate(bot_stand = as.numeric(scale(bot))) %>% 
        mutate(TPM80_ratio_stand = as.numeric(scale(TPM80_ratio)))
    
    
    ## Genetic model (2): Allelic richness ~ het-exc (standardised)
    # model specification
    run_mod <- function(iter){
        MCMCglmm(num_alleles_mean ~ TPM80_ratio_stand, #
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=stats_mod_gen,nitt=nitt,burnin=burnin,thin=thin)
    }
    
    # model name
    mod_name <- "gen2_hetexc"
    
    # check if model is saved
    model_file_name <- paste0(mod_name, ".RData")
    
    # run 3 independent chains for model checking
    if (!file.exists(paste0("output/mcmcmodels/", model_file_name))){
        # run models
        set.seed(1234)
        models <- purrr::map(1:3, run_mod)
        saveRDS(models, file = paste0("output/mcmcmodels/", model_file_name))
    }
    # load models
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
    library(mcmcR2)
    set.seed(312)
    for (r2type in c("marginal", "conditional")) {
       
        model_file_name_R2 <- paste0(mod_name, "_R2_", r2type)
        
        if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
            
            R2_gen <- mcmcR2::partR2(mod_gen, c("TPM80_ratio_stand"), type = r2type,
                data = stats_mod_gen, inv_phylo = inv_phylo, prior = prior, 
                nitt = nitt, burnin = burnin, thin = thin)
            saveRDS(R2_gen, file = paste0("output/mcmcmodels/", model_file_name_R2, ".RData"))
        }
        
        R2_gen$R2 %>% write_delim(paste0("output/mcmcmodels/", model_file_name_R2 ,".txt"))
        
    }
    
    # SC
    R2_gen$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))
    # summary
    mod_sum <- summary(mod_gen)
    # beta coefficients
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
        write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))
    
}




# plots for paper ----------------------------------------------------------------------------------

# some controls
point_size <- 3.5
point_alpha <- 0.3

######## plot (1) het-exc vs genetic diversity #########

# re-run model, as we don't want standardized predictions for plotting
mod_plot_AR1 <- MCMCglmm(num_alleles_mean ~ TPM80_ratio, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_gen,nitt=110000,burnin=10000,thin=100)

pred_df_AR1 <- data.frame(num_alleles_mean =0,
    TPM80_ratio= seq(from = 0.1, to = 1, by = 0.05), 
    tip_label = all_stats$tip_label[1])

mod_preds_AR1 <- data.frame(predict(mod_plot_AR1, pred_df_AR1, interval = "confidence")) %>% 
             mutate(TPM80_ratio= seq(from = 0.1, to = 1, by = 0.05))

AR1_beta <- read_delim("output/mcmcmodels/gen2_hetexc_beta.txt", delim = " ")
AR1_R2 <- read_delim("output/mcmcmodels/gen2_hetexc_R2_marginal.txt", delim = " ")

p1 <- ggplot(aes(x = TPM80_ratio, y = num_alleles_mean), data = all_stats) +
    geom_line(data = mod_preds_AR1, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    scale_y_continuous(breaks = seq(from = 2, to = 10, by = 2), limits = c(1,10)) +
    scale_x_continuous(breaks = seq(from = 0.2, to = 1, by = 0.2)) +
    xlab("Heterozygosity-excess") + #Allelic richness
    ylab("Allelic richness") +
    annotate("text", x = 0.33, y = 2, label = "R^2 == '0.02 [0, 0.2]'", 
             parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    annotate("text", x = 0.33, y = 1.3, label = "beta == '-0.18 [-0.92, 0.58]'", 
             parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.9,0.1,0.25,0.9), "cm"),
        axis.line = element_line(colour = "#cccccc"),
        axis.ticks = element_line(colour = "#cccccc")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
        )
p1


if (plotting) {

######## plot (2) ABCprob(bot) vs genetic diversity #########

# re-run model, as we don't want standardized predictions for plotting
mod_plot_AR2 <- MCMCglmm(num_alleles_mean ~ bot, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_gen,nitt=110000,burnin=10000,thin=100)

pred_df_AR2 <- data.frame(num_alleles_mean =0,
    bot= seq(from = 0, to = 1, by = 0.05), 
    tip_label = all_stats$tip_label[1])

mod_preds_AR2 <- data.frame(predict(mod_plot_AR2, pred_df_AR2, interval = "confidence")) %>% 
    mutate(bot= seq(from = 0, to = 1, by = 0.05))

AR2_beta <- read_delim("output/mcmcmodels/gen1_bot_beta.txt", delim = " ")
AR2_R2 <- read_delim("output/mcmcmodels/gen1_bot_R2_marginal.txt", delim = " ")

p2 <- ggplot(aes(x = bot, y = num_alleles_mean), data = all_stats) +
    geom_line(data = mod_preds_AR2, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    scale_y_continuous(breaks = seq(from = 2, to = 10, by = 2), limits = c(1,10)) +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
    xlab("Bottleneck model probability (ABC)") + #Allelic richness
    ylab("") +
    annotate("text", x = 0.24, y = 2, label = "R^2 == '0.49 [0.2, 0.7]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    annotate("text", x = 0.24, y = 1.3, label = "beta == '-1.3 [-1.8, -0.78]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.9,0.5,0.25,0.1), "cm"),
        axis.line = element_line(colour = "#cccccc"),
        axis.ticks = element_line(colour = "#cccccc")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    )
p2

p_bot_vs_div <- plot_grid(p1, p2, ncol = 2, labels = c("A", "B"))
p_bot_vs_div

# save it!
#cowplot::save_plot(filename = "other_stuff/figures/figures_final/fig2_bot_vs_div.jpg", 
#                                plot = p_bot_vs_div, base_height = 3, base_width = 7)




### For Supplementary: bot~hetexc ------------------------------------------------------------------

# model --------------------------------------------------------------------------------------------
### modeling variable at the beginning
if(modeling){
    
    # standardize variables (by 1 sd, as no binary predictor is involved)
    stats_mod_gen <- all_stats %>% 
        mutate(bot_stand = as.numeric(scale(bot))) %>% 
        mutate(TPM80_ratio_stand = as.numeric(scale(TPM80_ratio)))
    
    # model specification
    run_mod <- function(iter){
        MCMCglmm(bot ~ TPM80_ratio_stand, #
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=stats_mod_gen,nitt=nitt,burnin=burnin,thin=thin)
    }
    
    # model name
    mod_name <- "gen3_bot_vs_hetexc"
    
    # check if model is saved
    model_file_name <- paste0(mod_name, ".RData")
    
    # run 3 independent chains for model checking
    if (!file.exists(paste0("output/mcmcmodels/", model_file_name))){
        # run models
        set.seed(1234)
        models <- purrr::map(1:3, run_mod)
        saveRDS(models, file = paste0("output/mcmcmodels/", model_file_name))
    }
    # load models
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
    library(mcmcR2)
    set.seed(3224)
    for (r2type in c("marginal", "conditional")) {
        
        model_file_name_R2 <- paste0(mod_name, "_R2_", r2type)
        
        if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
        
            R2_gen <- mcmcR2::partR2(mod_gen, c("TPM80_ratio_stand"), type = r2type,
                data = stats_mod_gen, inv_phylo = inv_phylo, prior = prior, 
                nitt = 110000, burnin = 10000, thin = 100)
            saveRDS(R2_gen, file = paste0("output/mcmcmodels/", model_file_name_R2, ".RData"))
        }
        
        R2_gen$R2 %>% write_delim(paste0("output/mcmcmodels/", model_file_name_R2 ,".txt"))
        
    }
    
    # SC
    R2_gen$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

    # summary
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
        write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))

    
}


# plot
# some controls
point_size <- 3.5
point_alpha <- 0.3

# re-run model, as we don't want standardized predictions for plotting
mod_plot_bot_sup <- MCMCglmm(bot ~ TPM80_ratio, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_gen,nitt=110000,burnin=10000,thin=100)

pred_df_bot_sup <- data.frame(bot=0,
    TPM80_ratio= seq(from = 0.1, to = 1, by = 0.05), 
    tip_label = all_stats$tip_label[1])

mod_preds_bot_sup <- data.frame(predict(mod_plot_bot_sup, pred_df_bot_sup, interval = "confidence")) %>% 
    mutate(TPM80_ratio= seq(from = 0.1, to = 1, by = 0.05))

bot_sup_beta <- read_delim("output/mcmcmodels/gen3_bot_vs_hetexc_beta.txt", delim = " ")
bot_sup_R2 <- read_delim("output/mcmcmodels/gen3_bot_vs_hetexc_R2_marginal.txt", delim = " ")

p3 <- ggplot(aes(x = TPM80_ratio, y = bot), data = all_stats) +
    geom_line(data = mod_preds_bot_sup, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
   # scale_y_continuous(breaks = seq(from = 2, to = 10, by = 2), limits = c(1,10)) +
    #scale_x_continuous(breaks = seq(from = 0.2, to = 1, by = 0.2)) +
    xlab(expression(Heterozygosity-excess ~ "("~prop[het-exc]~")")) + #Allelic richness
    ylab(expression(ABC~bottleneck~probability~(p[bot]))) +
    annotate("text", x = 0.25, y = 1, label = "R^2 == '0.32 [0.03, 0.59]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    annotate("text", x = 0.25, y = 0.9, label = "beta == '0.17 [0.04, 0.28]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.9,0.1,0.25,0.9), "cm"),
        axis.line = element_line(colour = "#cccccc"),
        axis.ticks = element_line(colour = "#cccccc")
        #axis.title.x = element_text(margin = margin(t = 10)),
        #axis.title.y = element_text(margin = margin(r = 10))
    )
p3

ggsave(filename = "other_stuff/figures/figures_final/sup_fig2_bot_vs_bot_30.jpg", height = 3,
       width = 4.4)








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

}