## This script produces the figure for the models of bottleneck signatures
# explained by life-history traits.

# phylogenetic comparative analysis
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
plotting <- FALSE

# load data and prepare mixed models

# load (modified) phylogeney. 26 species from 10ktrees plus 3 subspecies of ringed seal
tree_final <- read.tree("data/raw/phylogeny/30_species_10ktrees_final.tre")

# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_30_modeling.csv"))

# phylogenetic mixed model preparation

# construct inverse phylo matrix and priors
inv_phylo <- inverseA(tree_final, nodes="TIPS",scale=FALSE)$Ainv #,scale=TRUE
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

# standardise variables by 2 sd for further analyses
stats_mod_hetexc <- all_stats %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logharem_size = (log(harem_size)),
        life_span_years = (life_span_years - mean(life_span_years)) / (2*sd(life_span_years)),
        breeding_season_length = (breeding_season_length - mean(breeding_season_length)) / (2*sd(breeding_season_length))) %>% 
    mutate(logharem_size = (logharem_size - mean(logharem_size) / (2*sd(logharem_size)))) %>% 
    mutate(harem_size = (harem_size - mean(harem_size) / (2*sd(harem_size)))) 

stats_mod_hetexc <- stats_mod_hetexc %>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))

## model 1: het excess vs. LH with SSD -------------------------------------------------------------

# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)

# final
# nitt <- 1100000
# burnin <- 100000
# thin <- 1000
# test
 nitt <- 110000
 burnin <- 10000
 thin <- 100

# specify model
run_mod <- function(iter){
    MCMCglmm(TPM80_ratio ~ SSD + BreedingType + Generation_time + breeding_season_length, # , #+ Abundance BreedingType  + BreedingType + Generation_time
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_hetexc,nitt=nitt,burnin=burnin,thin=thin)
}

# model name
# plus means plus generation time and breeding season length
mod_name <- "hetexc_vs_lh_SSD_plus"

# check if model is saved
model_file_name <- paste0(mod_name, ".RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name))){
    # run models
    set.seed(1234)
    models <- purrr::map(1:3, run_mod)
    saveRDS(models, file = paste0("output/mcmcmodels/", model_file_name))
}
# read models from file is existing
models <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name))

# model checks according to holgers paper
plot(mcmc.list(models[[1]]$Sol, models[[2]]$Sol, models[[3]]$Sol))
gelman.diag(mcmc.list(models[[1]]$Sol, models[[2]]$Sol, models[[3]]$Sol))

# choose first model for R2 etc.
mod_hetexc <- models[[1]]

# summary
options(scipen=999)
summary(mod_hetexc)

# visually inspecting chain convergence
plot(mod_hetexc$Sol)
plot(mod_hetexc$VCV)
autocorr(mod_hetexc$Sol)
autocorr(mod_hetexc$VCV)

# variation explained by phylogeny
var_phy <- mod_hetexc$VCV[, "tip_label"] / (mod_hetexc$VCV[, "tip_label"] + mod_hetexc$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

# R2s
library(mcmcR2)
set.seed(312)
for (r2type in c("marginal", "conditional")) {
    
    # create file name
    model_file_name_R2 <- paste0(mod_name, "_R2_", r2type)
    
    # partition R2, calculate SC
    R2_hetexc <- mcmcR2::partR2(mod_hetexc, type = r2type,
        partvars = c("SSD", "BreedingType", "Generation_time", "breeding_season_length"),
        data = stats_mod_hetexc, inv_phylo = inv_phylo, prior = prior, 
        nitt = nitt, burnin = burnin, thin = thin)
    
    saveRDS(R2_hetexc, file = paste0("output/mcmcmodels/", model_file_name_R2, ".RData"))
    R2_hetexc$R2 %>% write_delim(paste0("output/mcmcmodels/", model_file_name_R2 ,".txt"))
}

# SC
R2_hetexc$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary of model (betas) to file
mod_hetexc %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_hetexc$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_hetexc$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))


## model 2: bot vs LH with SSD ---------------------------------------------------------------------

# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)

stats_mod_bot <- all_stats %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logharem_size = (log(harem_size)),
        life_span_years = (life_span_years - mean(life_span_years)) / (2*sd(life_span_years)),
        breeding_season_length = (breeding_season_length - mean(breeding_season_length)) / (2*sd(breeding_season_length))) %>% 
    mutate(logharem_size = (logharem_size - mean(logharem_size) / (2*sd(logharem_size)))) %>% 
    mutate(harem_size = (harem_size - mean(harem_size) / (2*sd(harem_size)))) 

stats_mod_bot <- stats_mod_hetexc %>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))

# specify model
run_mod <- function(iter){
    MCMCglmm(bot ~ SSD + BreedingType + Generation_time + breeding_season_length, # , #+ Abundance BreedingType  + BreedingType + Generation_time
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_bot,nitt=nitt,burnin=burnin,thin=thin)
}

# check if model is saved
# model name
mod_name <- "bot_vs_lh_SSD_plus"

# check if model is saved
model_file_name <- paste0(mod_name, ".RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name))){
    # run models
    set.seed(1234)
    models <- purrr::map(1:3, run_mod)
    saveRDS(models, file = paste0("output/mcmcmodels/", model_file_name))
}

# read model from file if existing
models <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name))

# model checks according to holgers paper
plot(mcmc.list(models[[1]]$Sol, models[[2]]$Sol, models[[3]]$Sol))
gelman.diag(mcmc.list(models[[1]]$Sol, models[[2]]$Sol, models[[3]]$Sol))

# one model
mod_bot <- models[[1]]

# summary
summary(mod_bot)

# visually inspecting chain convergence
plot(mod_bot$Sol)
plot(mod_bot$VCV)
autocorr(mod_bot$Sol)
autocorr(mod_bot$VCV)

# variation explained by phylogeny
var_phy <- mod_bot$VCV[, "tip_label"] / (mod_bot$VCV[, "tip_label"] + mod_bot$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

# R2s
library(mcmcR2)
set.seed(312)
for (r2type in c("marginal", "conditional")) {
    
    # file name 
    model_file_name_R2 <- paste0(mod_name, "_R2_", r2type)
    # calculate R2, SC
    R2_hetexc <- mcmcR2::partR2(mod_bot, type = r2type,
        partvars = c("SSD", "BreedingType", "Generation_time", "breeding_season_length"),
        data = stats_mod_hetexc, inv_phylo = inv_phylo, prior = prior, 
        nitt = nitt, burnin = burnin, thin = thin)
    # save to file
    saveRDS(R2_hetexc, file = paste0("output/mcmcmodels/", model_file_name_R2, ".RData"))
    R2_hetexc$R2 %>% write_delim(paste0("output/mcmcmodels/", model_file_name_R2 ,".txt"))
}

# SC
R2_hetexc$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary of the model (betas) to file
mod_bot %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_bot$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_bot$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))




# plots for SSD------------------------------------------------------------------------------------
if (plotting) {
    

point_size = 3.5
# boxplot 1 - bot -------------------
p1 <- ggplot(aes(BreedingType, bot), data = all_stats) +
    geom_boxplot(alpha = 0.5, col = "darkgrey",  size = 0.7, width = 0.7, aes(fill = BreedingType), outlier.shape = NA) + #
   # geom_point(size = point_size, alpha = point_alpha, aes(color = BreedingType)) + # abc_out
    geom_jitter(size = point_size, alpha = 0.6, shape = 21, col = "black", aes(fill = BreedingType), width = 0.2) +
    theme_martin(base_family = "Arial", highlight_family = "Arial") + #Arial
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab("Breeding Habitat") +
    ylab(expression(ABC~bottleneck~probability~"("~p[bot]~")")) +
    guides(fill=FALSE, color = FALSE) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0, 1.05)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
       # panel.grid.major.y = element_line(colour = "lightgrey", size = 0.2),
        plot.margin = unit(c(0.9,0.5,0.25,0.1), "cm"),
        axis.line.x = element_line(colour = "#cccccc"),
        axis.line.y = element_line(colour = "#cccccc"),
        #axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "#cccccc")) 
  
p1

# boxplot 2 - TPM80 -------------------
p2 <- ggplot(aes(BreedingType, TPM80_ratio), data = all_stats) +
    geom_boxplot(alpha = 0.5, col = "darkgrey",  size = 0.7, width = 0.7, aes(fill = BreedingType), outlier.shape = NA) + #
    # geom_point(size = point_size, alpha = point_alpha, aes(color = BreedingType)) + # abc_out
    geom_jitter(size = point_size, alpha = 0.6, shape = 21, col = "black", aes(fill = BreedingType), width = 0.2) +
    theme_martin(base_family = "Arial", highlight_family = "Arial") +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab("Breeding Habitat") +
    ylab(expression(Heterozygosity-excess ~ "("~prop[het-exc]~")")) +
    guides(fill=FALSE, color = FALSE) +
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.05)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # panel.grid.major.y = element_line(colour = "lightgrey", size = 0.2),
        plot.margin = unit(c(0.9,0.5,0.25,0.1), "cm"),
        axis.line.x = element_line(colour = "#cccccc"),
        axis.line.y = element_line(colour = "#cccccc"),
        #axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "#cccccc")
        )

p2

plot_grid(p1, p2)


# scatterplot -----------------------
mod_hetexc_plot <- MCMCglmm(TPM80_ratio ~ SSD, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_hetexc,nitt=110000,burnin=10000,thin=100)

pred_df_hetexc <- data.frame(TPM80_ratio = 0,
    SSD = seq(from = 0.5, to = 8, by = 0.1),
    tip_label = all_stats$tip_label[1])

mod_preds_hetexc <- data.frame(predict(mod_hetexc_plot, pred_df_hetexc, interval = "confidence")) %>% 
    mutate(SSD = seq(from = 0.5, to = 8, by = 0.1))

set.seed(57)
p3 <- ggplot(aes(x = SSD, y = TPM80_ratio), data = all_stats) +
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 0.2, alpha = 0.5) +
    geom_point(size = point_size, alpha = 0.7,  aes(col = bot)) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    theme_martin(base_family = "Arial", highlight_family = "Arial") +
    xlab("Sexual Size Dimorphism (SSD)") +
    ylab(expression(Heterozygosity-excess ~ "("~prop[het-exc]~")")) +
    #ylab("Heterozygosity-excess") +
    # scale_x_continuous(breaks = log(c(1,2,3,4,6,8,10,15,20,30,40,50)), labels = c(1,2,3,4,6,8,10,15,20,30,40,50)) + #breaks = c(1,2,3,4,5,6,7,8)
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 1.2, alpha = 0.5, color = "grey") + 
    scale_fill_manual(values = c("cornflowerblue", "goldenrod")) +
    #ylab("Heterozygosity-excess") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.9,0.2,0.25,0.1), "cm"),
        axis.line = element_line(colour = "#cccccc"),
        axis.ticks = element_line(colour = "#cccccc"),
        legend.position = c(0.77,0.3),
        legend.title=element_text(size=10)
        #legend.text.align = - 1
    ) +
    guides(color = guide_colorbar(
        title = expression(p[bot]~"("~ABC~")"), 
        #title = expression(paste("ABC bottleneck\nprobability (Pbot)")),
        #title = expression(atop("ABC bottleneck", paste("probability p"[bot]))),
        barwidth = 0.7, barheight = 3, 
        title.position = "right")) + #, label.position = "bottom"
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1)) +
    scale_color_distiller(palette = "RdBu",
        direction = -1, 
        #name = expression(atop("ABC bottleneck", paste("probability p"[bot]))),
        labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    geom_text_repel(label = all_stats$short,size = 2.5, alpha = 1, color = "grey50",#  aes(label = common) , 
        segment.alpha= 1,  box.padding = unit(0.4, "lines"), point.padding = unit(0.7, "lines"),
        segment.size = 0.1,  force = 1, min.segment.length = unit(0.01, "lines"))

p3
# expression(atop("ABC bottleneck", paste("probability p" [reported]))))
# expression(p[bot])
# expression(Channel~Density~(km/km^2))
# bquote('Assimilation ('*mu~ 'mol' ~CO[2]~ m^-2~s^-1*')')
# ggsave(filename = "other_stuff/figures/figures_final/fig3_hetexc_vs_lh.jpg", width = 9, height = 3)


# load model output
mod_beta <- read_delim("output/mcmcmodels/hetexc_vs_lh_SSD_plus_beta.txt", delim = " ")
mod_R2 <- read_delim("output/mcmcmodels/hetexc_vs_lh_SSD_plus_R2_marginal.txt", delim = " ")
mod_SC <- read_delim("output/mcmcmodels/hetexc_vs_lh_SSD_plus_SC.txt", delim = " ")

mod_beta2 <- read_delim("output/mcmcmodels/bot_vs_lh_SSD_plus_beta.txt", delim = " ")
mod_R22 <- read_delim("output/mcmcmodels/bot_vs_lh_SSD_plus_R2_marginal.txt", delim = " ")
mod_SC2 <- read_delim("output/mcmcmodels/bot_vs_lh_SSD_plus_SC.txt", delim = " ")

# beta coefficients
mod_out <- mod_beta[-1, c("components", "post_median", "lower", "upper")]
names(mod_out) <- c("comps", "pe", "cilow", "cihigh")
mod_out$comps <- fct_relevel(mod_out$comps, "Generation_time", "breeding_season_length", "SSD", "BreedingTypeice")
mod_out2 <- mod_beta2[-1, c("components", "post_median", "lower", "upper")]
names(mod_out2) <- c("comps", "pe", "cilow", "cihigh")
mod_out2$comps <- fct_relevel(mod_out2$comps, "Generation_time", "breeding_season_length", "SSD", "BreedingTypeice")

p4 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0, position=position_nudge(y = 0.1)) +
    geom_point(size = 2.5, shape = 21, col = "black", fill = "grey69", position=position_nudge(y = 0.1)) +
    geom_errorbarh(data = mod_out2, alpha=0.4, color="black",height = 0, position=position_nudge(y = -0.1)) +
    geom_point(data = mod_out2, size = 2.5, shape = 21, col = "black",fill = "white", position=position_nudge(y = -0.1)) +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin(base_family = "Arial", highlight_family = "Arial") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y =  element_blank(),
        plot.margin = unit(c(1,0,0.33,0.1), "cm"),
        axis.title.x=element_text(margin=margin(t=12))) +
    scale_x_continuous(breaks = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2)) +
    # scale_y_discrete(labels = c("Breeding\nhabitat",
    #     "SSD", "Generation Time", "Breeding Season Length")) +
    xlab(expression(paste("Standardized ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p4

# SC
mod_out_SC <- mod_SC[, c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC ) <- c("comps", "pe", "cilow", "cihigh")
mod_out_SC$comps <- fct_relevel(mod_out_SC$comps, "Generation_time", "breeding_season_length", "SSD", "BreedingType")
mod_out_SC2 <- mod_SC2[, c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC2) <- c("comps", "pe", "cilow", "cihigh")
mod_out_SC2$comps <- fct_relevel(mod_out_SC2$comps, "Generation_time", "breeding_season_length", "SSD", "BreedingType")
# structure coefficients
p5 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_SC) + 
    geom_errorbarh(alpha=0.4, color="black",height = 0, position=position_nudge(y = 0.1)) +
    geom_point(size = 2.5, shape = 21, col = "black", fill = "grey69", position=position_nudge(y = 0.1)) +
    geom_errorbarh(data = mod_out_SC2, alpha=0.4, color="black",height = 0, position=position_nudge(y = -0.1)) +
    geom_point(data = mod_out_SC2, size = 2.5, shape = 21, col = "black",fill = "white", position=position_nudge(y = -0.1)) +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin(base_family = "Arial", highlight_family = "Arial") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = c(0.5), size = 8),
        plot.margin = unit(c(1,0.2,0.49,0), "cm"),
        axis.text.x=element_text(margin=margin(t=6))
        ) +
    scale_y_discrete(labels = c("Generation\nTime", "Breeding\nSeason Length",
                                 "SSD", "Breeding\nhabitat"
        )) +
    xlab(expression(paste("Structure coefficient", " r(", hat(Y),",x)") )) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p5

# R2
col_legend <- "#969696"
mod_out_R2 <- mod_R2[, c("combinations", "medianR2", "lower", "upper")][1:5, ]
names(mod_out_R2) <- c("comps", "pe", "cilow", "cihigh")
mod_out_R2$comps <- factor(mod_out_R2$comps, levels = c("Generation_time", "breeding_season_length", "SSD", "BreedingType", "full model"))

mod_out_R22 <- mod_R22[, c("combinations", "medianR2", "lower", "upper")][1:5, ]
names(mod_out_R22) <- c("comps", "pe", "cilow", "cihigh")
mod_out_R22$comps <- factor(mod_out_R22$comps, c("Generation_time", "breeding_season_length", "SSD", "BreedingType", "full model"))
# mod_out_R2$comps <- rev(fct_inorder(factor(mod_out_R2$comps)))
# structure coefficients
p6 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_R2 ) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0, position=position_nudge(y = 0.1)) +
    geom_point(size = 2.5, shape = 21, col = "black", fill = "grey69", position=position_nudge(y = 0.1)) +
    geom_errorbarh(data = mod_out_R22, alpha=0.4, color="black",height = 0, position=position_nudge(y = -0.1)) +
    geom_point(data = mod_out_R22, size = 2.5, shape = 21, col = "black",fill = "white", position=position_nudge(y = -0.1)) +
    
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin(base_family = "Arial", highlight_family = "Arial") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.text.y = element_text(hjust = c(0.5), size = 8),
        plot.margin = unit(c(0.3, 0.2, 0.39, 0.3), "cm")) +
    scale_y_discrete(labels = c("Generation\nTime", "Breeding\nSeason Length",
                            "SSD", "Breeding\nhabitat", "Full model")) +
    xlab(expression(paste(R^{2}))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
    #annotate("segment", x = 0.9, xend = 0.9, y = 0.8, yend = 3.5, color = col_legend) +
    #annotate("text", x = 0.98, xend = 0.95, y = 2, yend = 3, color = col_legend, label = c("common"), angle = 270) +
    annotate("segment", x = 1.2, xend = 1.2, y = 0.8, yend = 4.2, color = col_legend) +
    annotate("text", x = 1.3, xend = 1.3, y = 2.6, yend = 2.2, color = col_legend, 
             size = 3, label = c("unique"), angle = 270) +
    annotate("segment", x = 1.2, xend = 1.2, y = 5.4, yend = 4.5, color = col_legend) +
    annotate("text", x = 1.3, xend = 1.3, y = 4.6, yend = 4.5, color = col_legend, 
        size = 3, label = c("marginal"), angle = 270)
p6


p_top <- plot_grid(p2, p1, p3, nrow = 1, rel_widths = c(1,1,2), labels = c("A", "B", "C"), label_x = 0.1)
p_top

p_bot <- plot_grid(p6, p4, p5, nrow = 1, rel_widths = c(1.6,1.35,1.6),
    labels = c("D","E","F"))
p_bot

p_final <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(2,1.5))
p_final

# jpg
ggsave('other_stuff/figures/figures_final/figures_final_editing/fig3_bot_vs_lh_plus.jpg',p_final,  width=9, height=5.5)

# pdf
ggsave('other_stuff/figures/figures_final/figures_final_editing/fig3_bot_vs_lh_plus.pdf',p_final,  width=9, height=5.5)
Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts("other_stuff/figures/figures_final/fig3_bot_vs_lh.pdf")

library(extrafont)
# install.packages("extrafont")

}






# # plot for marine mammal conference
# 
# p1 <- ggplot(aes(logharem_size, TPM80_ratio), data = all_stats) +
#     geom_point(size = 3.5, alpha = 0.3) + # abc_out
#     geom_point(size = 3.5, alpha = 0.8, shape = 21, col = "black") + #, aes(fill = BreedingType)
#     #scale_color_viridis(option = "magma", direction = -1,
#     #    name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
#     theme_martin() +
#     xlab("Harem Size") +
#     ylab("Heterozygosity-excess") +
#     scale_x_continuous(breaks = log(c(1,5,10,20,50)), labels = c(1,5,10,20,50)) + #breaks = c(1,2,3,4,5,6,7,8)
#     geom_line(data = mod_preds_hetexc, aes(y = fit), size = 1, alpha = 0.5) + 
#     #ylab("Heterozygosity-excess") +
#     theme(#panel.grid.minor = element_blank(),
#         #panel.grid.major = element_blank(),
#         #panel.grid.major.x = element_blank(),
#         plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm") ,
#         #axis.line.x = element_line(color="#cccccc", size=0.3),
#         #axis.ticks.x = element_line(color="#cccccc", size=0.3),
#         #axis.text.x = element_text(margin = margin(t = 5)),
#         #axis.title.y = element_text(margin = margin(r = 10))
#         legend.direction = "vertical",
#         legend.position = c(0.85,0.25),
#         legend.title=element_text(size=10),
#         axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))
#     ) +
#     scale_color_distiller(palette = "RdBu",
#         direction = -1,
#         name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
#     
#     guides(color = guide_colorbar(barwidth = 0.5, barheight = 5, 
#         title.position = "left")) + #, label.position = "bottom"
#     scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) 
#    # geom_text_repel(aes(label = short),size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
#    #       segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
#     #    segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines"))
# p1
# 
# ggsave("figures/SMM/p2_smm.jpg", width = 4, height = 3.5)
# 
# 
# 
# p2 <- ggplot(aes(logharem_size, TPM80_ratio), data = all_stats) +
#     geom_point(size = 3.5, alpha = 0.3) + # abc_out
#     geom_point(size = 3.5, alpha = 0.8, shape = 21, aes(fill = BreedingType)) + #, aes(fill = BreedingType)
#     #scale_color_viridis(option = "magma", direction = -1,
#     #    name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
#     theme_martin() +
#     xlab("Harem Size") +
#     ylab("Heterozygosity-excess") +
#     scale_x_continuous(breaks = log(c(1,5,10,20,50)), labels = c(1,5,10,20,50)) + #breaks = c(1,2,3,4,5,6,7,8)
#     geom_line(data = mod_preds_hetexc, aes(y = fit), size = 1, alpha = 0.5) + 
#     #ylab("Heterozygosity-excess") +
#     theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.major.x = element_blank(),
#         plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm") ,
#         axis.line.x = element_line(color="#cccccc", size=0.3),
#         axis.ticks.x = element_line(color="#cccccc", size=0.3),
#         axis.text.x = element_text(margin = margin(t = 5)),
#         #axis.title.y = element_text(margin = margin(r = 10))
#         legend.direction = "vertical",
#         legend.position = c(0.85,0.15),
#         legend.title=element_text(size=10),
#         axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))
#     ) +
#     scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
#     guides(color = guide_colorbar(barwidth = 0.5, barheight = 5, 
#         title.position = "left")) + #, label.position = "bottom"
#     scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) 
# # geom_text_repel(aes(label = short),size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
# #       segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
# #    segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines"))
# p2
# 
# ggsave("figures/SMM/p3_smm.jpg", width = 4, height = 3.5)






## OTHER models 1: het excess vs. LH with logharem_size ---------------------------------------------------

# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)

stats_mod_hetexc <- 
    all_stats %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logharem_size = (log(harem_size))) %>% 
    mutate(logharem_size = (logharem_size - mean(logharem_size) / (2*sd(logharem_size)))) %>% 
    mutate(harem_size = (harem_size - mean(harem_size) / (2*sd(harem_size)))) 

stats_mod_hetexc <- stats_mod_hetexc %>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))


second_var <- "logharem_size" # "logharem_size"

if (second_var == "logharem_size"){
    run_mod <- function(iter){
        MCMCglmm(TPM80_ratio ~ logharem_size + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=stats_mod_hetexc,nitt=1100000,burnin=100000,thin=1000)
    }
    
} else if (second_var == "SSD"){
    run_mod <- function(iter){
        MCMCglmm(TPM80_ratio ~ SSD + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=stats_mod_hetexc,nitt=1100000,burnin=100000,thin=1000)
    }
}


# check if model is saved
# model name
mod_name <- "hetexc_vs_lh_logharem_size"

# check if model is saved
model_file_name <- paste0(mod_name, ".RData")

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
mod_hetexc <- models[[1]]

# summary
summary(mod_hetexc)

# visually inspecting chain convergence
plot(mod_hetexc$Sol)
plot(mod_hetexc$VCV)
autocorr(mod_hetexc$Sol)
autocorr(mod_hetexc$VCV)

# variation explained by phylogeny
var_phy <- mod_hetexc$VCV[, "tip_label"] / (mod_hetexc$VCV[, "tip_label"] + mod_hetexc$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

# commonality analyses and R2
model_file_name_R2 <- paste0(mod_name, "_R2.RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_hetexc <- mcmcR2::partR2(mod_hetexc, partvars = c("SSD", "BreedingType"),
        data = stats_mod_hetexc, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS(R2_hetexc, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_hetexc <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))
R2_hetexc
# out <- mcmcR2::R2mcmc(mod_hetexc)
# out$partR2
R2_hetexc$R2 %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_R2" ,".txt"))
R2_hetexc$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary to file
mod_hetexc %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_hetexc$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_hetexc$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))



point_size <- 3.5
point_alpha <- 0.3


## OTHER models 2: bot vs LH with logharem_size -----------------

stats_mod_hetexc <- 
    all_stats %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logharem_size = (log(harem_size))) %>% 
    mutate(logharem_size = (logharem_size - mean(logharem_size) / (2*sd(logharem_size)))) %>% 
    mutate(harem_size = (harem_size - mean(harem_size) / (2*sd(harem_size)))) 

stats_mod_hetexc <- stats_mod_hetexc %>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))


second_var <- "logharem_size" # "logharem_size"

if (second_var == "logharem_size"){
    run_mod <- function(iter){
        MCMCglmm(bot ~ logharem_size + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=stats_mod_hetexc,nitt=1100000,burnin=100000,thin=1000)
    }
    
} else if (second_var == "SSD"){
    run_mod <- function(iter){
        MCMCglmm(bot ~ SSD + BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
            random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
            family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
            data=stats_mod_hetexc,nitt=1100000,burnin=100000,thin=1000)
    }
}


# check if model is saved
# model name
mod_name <- "bot_vs_lh_logharem_size"

# check if model is saved
model_file_name <- paste0(mod_name, ".RData")

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
mod_hetexc <- models[[1]]

# summary
summary(mod_hetexc)

# visually inspecting chain convergence
plot(mod_hetexc$Sol)
plot(mod_hetexc$VCV)
autocorr(mod_hetexc$Sol)
autocorr(mod_hetexc$VCV)

# variation explained by phylogeny
var_phy <- mod_hetexc$VCV[, "tip_label"] / (mod_hetexc$VCV[, "tip_label"] + mod_hetexc$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

# commonality analyses and R2
model_file_name_R2 <- paste0(mod_name, "_R2.RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_hetexc <- mcmcR2::partR2(mod_hetexc, partvars = c("SSD", "BreedingType"),
        data = stats_mod_hetexc, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS(R2_hetexc, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_hetexc <- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))
R2_hetexc
# out <- mcmcR2::R2mcmc(mod_hetexc)
# out$partR2
R2_hetexc$R2 %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_R2" ,".txt"))
R2_hetexc$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary to file
mod_hetexc %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_hetexc$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_hetexc$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))



# some controls
point_size <- 3.5
point_alpha <- 0.6

## extra model 5: SSD vs breeding habitat -------------
mod_SSD <- MCMCglmm(SSD  ~ BreedingType, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_bot,nitt=110000,burnin=10000,thin=1000)

out <- mcmcR2::partR2(mod_SSD, partvars = c("BreedingType"),
    data = stats_mod_bot, inv_phylo = inv_phylo, prior = prior, 
    nitt = 110000, burnin = 10000, thin = 100)


