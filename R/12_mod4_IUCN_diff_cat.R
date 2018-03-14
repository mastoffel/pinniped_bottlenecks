
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
library(purrr)
library(readr)
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

# data
stats_mod <- all_stats 

stats_mod <- stats_mod %>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))

# plot
stats_mod_IUCN <- stats_mod %>% 
    mutate(IUCN_binary = case_when(IUCN_rating == "vulnerable" ~ "concern",
        IUCN_rating == "near threatened" ~ "concern",
        IUCN_rating == "endangered" ~ "concern",
        IUCN_rating == "least concern" ~ "least concern"))

stats_mod_IUCN <- stats_mod %>% 
    mutate(IUCN_binary = case_when(IUCN_rating == "vulnerable" ~ "concern",
        IUCN_rating == "near threatened" ~ "least concern",
        IUCN_rating == "endangered" ~ "concern",
        IUCN_rating == "least concern" ~ "least concern"))
## modeling

# mod1 gendiv vs. IUCN -----------------------------------------------------------------------------
run_mod <- function(iter){
    MCMCglmm(num_alleles_mean ~ IUCN_binary, # , #+ Abundance BreedingType  + BreedingType + Generation_time
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_IUCN,nitt=1100000,burnin=100000,thin=1000)
}

# check if model is saved
# model name
mod_name <- "IUCN_gendiv"

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
mod_IUCN <- models[[1]]

# summary
summary(mod_IUCN)

# visually inspecting chain convergence
plot(mod_IUCN$Sol)
plot(mod_IUCN$VCV)
autocorr(mod_IUCN$Sol)
autocorr(mod_IUCN$VCV)

# commonality analyses and R2
model_file_name_R2 <- paste0(mod_name, "_R2.RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_genlh <- mcmcR2::partR2(mod_IUCN, partvars = c("IUCN_binary"),
        data = stats_mod, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS(R2_genlh, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_genlh<- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))
R2_genlh
# out <- mcmcR2::R2mcmc(mod_hetexc)
# out$partR2
R2_genlh$R2 %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_R2" ,".txt"))
R2_genlh$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary to file
mod_IUCN %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_IUCN$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_IUCN$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))








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

# mod2 hetexc vs. IUCN -----------------------------------------------------------------------------
run_mod <- function(iter){
    MCMCglmm(TPM80_ratio ~ IUCN_binary, # , #+ Abundance BreedingType  + BreedingType + Generation_time
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_IUCN,nitt=1100000,burnin=100000,thin=1000)
}

# check if model is saved
# model name
mod_name <- "IUCN_hetexc"

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
mod_IUCN <- models[[1]]

# summary
summary(mod_IUCN)

# visually inspecting chain convergence
plot(mod_IUCN$Sol)
plot(mod_IUCN$VCV)
autocorr(mod_IUCN$Sol)
autocorr(mod_IUCN$VCV)

# commonality analyses and R2
model_file_name_R2 <- paste0(mod_name, "_R2.RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_genlh <- mcmcR2::partR2(mod_IUCN, partvars = c("IUCN_binary"),
        data = stats_mod, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS(R2_genlh, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_genlh<- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))
R2_genlh
# out <- mcmcR2::R2mcmc(mod_hetexc)
# out$partR2
R2_genlh$R2 %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_R2" ,".txt"))
R2_genlh$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary to file
mod_IUCN %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_IUCN$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_IUCN$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))


# mod3 bot vs. IUCN ---------------------------------
run_mod <- function(iter){
    MCMCglmm(bot ~ IUCN_binary, # , #+ Abundance BreedingType  + BreedingType + Generation_time
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_IUCN,nitt=1100000,burnin=100000,thin=1000)
}

# check if model is saved
# model name
mod_name <- "IUCN_bot"

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
mod_IUCN <- models[[1]]

# summary
summary(mod_IUCN)

# visually inspecting chain convergence
plot(mod_IUCN$Sol)
plot(mod_IUCN$VCV)
autocorr(mod_IUCN$Sol)
autocorr(mod_IUCN$VCV)

# commonality analyses and R2
model_file_name_R2 <- paste0(mod_name, "_R2.RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_genlh <- mcmcR2::partR2(mod_IUCN, partvars = c("IUCN_binary"),
        data = stats_mod, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS(R2_genlh, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_genlh<- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))
R2_genlh
# out <- mcmcR2::R2mcmc(mod_hetexc)
# out$partR2
R2_genlh$R2 %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_R2" ,".txt"))
R2_genlh$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary to file
mod_IUCN %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_IUCN$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_IUCN$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))


# plotting
point_size <- 3.5

gendiv_beta <- read_delim("output/mcmcmodels/IUCN_gendiv_beta.txt", delim = " ")
gendiv_R2 <- read_delim("output/mcmcmodels/IUCN_gendiv_R2.txt", delim = " ")

p1 <-  ggplot(data = stats_mod_IUCN, aes(IUCN_binary, num_alleles_mean)) +
    geom_boxplot(alpha = 0.5, col = "darkgrey",  size = 0.5, width = 0.7,  outlier.shape = NA) + #aes(fill = BreedingType),
    geom_jitter(size = point_size, alpha = 0.6, shape = 21, col = "black", aes(fill = BreedingType), width = 0.2) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    scale_color_manual(values = c("#d8b365", "cornflowerblue")) +
    scale_fill_manual(values = c("#d8b365", "cornflowerblue")) +
    xlab(" ") +
    ylab(expression(paste(Allelic~richness~"("~A[r]~")"))) +
   # ylab(expression(Allelic~richness~(A[r]))) +
    scale_y_continuous(limits=c(2, 12), breaks = c(2,4,6,8,10)) +
    annotate("text", x = 1.5, y = 12, label = "R^2 == '0.19 [0, 0.33]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    annotate("text", x = 1.5, y = 11.1, label = "beta == '1.29 [-0.08, 2.62]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.2),
        plot.margin = unit(c(2,0.33,0.5,0.2), "cm"),
        axis.line.x = element_line(colour = "#cccccc"),
        #axis.line.y = element_line(colour = "#cccccc"),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "#cccccc"),
        legend.position = "none") 
p1    

hetexc_beta <- read_delim("output/mcmcmodels/IUCN_hetexc_beta.txt", delim = " ")
hetexc_R2 <- read_delim("output/mcmcmodels/IUCN_hetexc_R2.txt", delim = " ")

p2 <-  ggplot(data = stats_mod_IUCN, aes(IUCN_binary, TPM80_ratio)) +
    geom_boxplot(alpha = 0.5, col = "darkgrey",  size = 0.5, width = 0.7,  outlier.shape = NA) + #aes(fill = BreedingType),
    geom_jitter(size = point_size, alpha = 0.6, shape = 21, col = "black", aes(fill = BreedingType), width = 0.2) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    scale_color_manual(values = c("#d8b365", "cornflowerblue")) +
    scale_fill_manual(values = c("#d8b365", "cornflowerblue")) +
    xlab("IUCN status") +
    ylab(expression(atop("Heterozygosity-excess", paste("("~prop[het-exc]~")")))) +
    #ylab(expression(Heterozygosity-excess ~ "("~prop[het-exc]~")")) +
    annotate("text", x = 1.5, y = 1.2, label = "R^2 == '0.02 [0, 0.2]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    annotate("text", x = 1.5, y = 1.11, label = "beta == '0.0 [-0.16, 0.14]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1), limits = c(0.15, 1.2)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.2),
       # plot.margin = unit(c(0.9,0.5,0.25,0.1), "cm"),
        axis.line.x = element_line(colour = "#cccccc"),
        #axis.line.y = element_line(colour = "#cccccc"),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "#cccccc"),
        plot.margin = unit(c(0.5,0.3,0.3,0.2), "cm") ,
        axis.title.x= element_text(margin=margin(t=0.5, unit = "cm")),
        legend.position="top")
       # legend.position = "none") 
p2

bot_beta <- read_delim("output/mcmcmodels/IUCN_bot_beta.txt", delim = " ")
bot_R2 <- read_delim("output/mcmcmodels/IUCN_bot_R2.txt", delim = " ")

p3 <-  ggplot(data = stats_mod_IUCN, aes(IUCN_binary, bot)) +
    geom_boxplot(alpha = 0.5, col = "darkgrey",  size = 0.5, width = 0.7,  outlier.shape = NA) + #aes(fill = BreedingType),
    geom_jitter(size = point_size, alpha = 0.6, shape = 21, col = "black", aes(fill = BreedingType), width = 0.2) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    scale_color_manual(values = c("#d8b365", "cornflowerblue")) +
    scale_fill_manual(values = c("#d8b365", "cornflowerblue")) +
    xlab(" ") +
    #ylab(expression("ABC bottleneck\nprobability"~(p[bot]))) +
    ylab(expression(atop("ABC bottleneck", paste("probability (p"[bot]~")")))) +
    annotate("text", x = 1.5, y = 1.2, label = "R^2 == '0.1 [0, 0.38]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    annotate("text", x = 1.5, y = 1.11, label = "beta == '-0.18 [-0.4, 0.02]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0, 1.2)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.2),
        # plot.margin = unit(c(0.9,0.5,0.25,0.1), "cm"),
        axis.line.x = element_line(colour = "#cccccc"),
        #axis.line.y = element_line(colour = "#cccccc"),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "#cccccc"),
        plot.margin = unit(c(1.9,0.3,0.3,0.2), "cm") ,
        axis.title.x= element_text(margin=margin(t=0.5, unit = "cm")),
        legend.position="")
# legend.position = "none") 
p3

p_final <- plot_grid(p1,p2,p3, nrow = 1, labels = c("A", "B", "C"))
ggsave(p_final, filename = "other_stuff/figures/figures_final/IUCN.jpg", width = 9, height = 3.5)


ggsave(p_final, filename = "other_stuff/figures/figures_final/IUCN.pdf", width = 9, height = 3.5)

Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts("other_stuff/figures/figures_final/IUCN.pdf")





p2 <-  ggplot(data = stats_mod_IUCN, aes(IUCN_binary, TPM80_ratio)) +
    geom_boxplot(alpha = 0.5, col = "darkgrey",  size = 0.5, width = 0.7,  outlier.shape = NA) + #aes(fill = BreedingType),
    geom_jitter(size = point_size, alpha = 0.6, shape = 21, col = "black", aes(fill = BreedingType), width = 0.2) +
    theme_martin() +
    scale_color_manual(values = c("cornflowerblue", "#d8b365")) +
    scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
    xlab(" ") +
    ylab("Heterozygosity-excess") +
    annotate("text", x = 0.9, y = 0.95, label = "R^2 == '0.19 [0, 0.33]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    annotate("text", x = 0.9, y = 0.88, label = "beta == '1.29 [-0.08, 2.62]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    scale_y_continuous(breaks = c(0.2,0.4,0.6,0.8,1), limits = c(0.15, 1)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", size = 0.2),
        plot.margin = unit(c(0.9,0.5,0.25,0.1), "cm"),
        axis.line.x = element_line(colour = "#cccccc"),
        #axis.line.y = element_line(colour = "#cccccc"),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "#cccccc"),
        legend.position = "none") 
p2



