## This script produces the figure for the models of bottleneck signatures
# explained by life-history traits.

# phylogenetic comparative analysis
library(patchwork)
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


stats_mod_genlh <- all_stats %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        logAbundance = ((logAbundance - mean(logAbundance)) / (2*sd(logAbundance))),
        Generation_time = (Generation_time - mean(Generation_time)) / (2*sd(Generation_time)),
        SSD = (SSD - mean(SSD)) / (2*sd(SSD)),
        bot = (bot - mean(bot)) / (2*sd(bot)),
        TPM80_ratio = (TPM80_ratio - mean(TPM80_ratio)) / (2*sd(TPM80_ratio)),
        logharem_size = (log(harem_size))) %>% 
    mutate(logharem_size = (logharem_size - mean(logharem_size)) / (2*sd(logharem_size))) %>% 
    mutate(harem_size = (harem_size - mean(harem_size)) / (2*sd(harem_size)))

stats_mod_genlh <- stats_mod_genlh %>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))

## model 1: gen div vs. LH with SSD -------------------------------------------------------

# standardize by 2 sd to make estimates comparable with BreedingHabitat variable(
# Gelman (2008), Schielzeth (2014)

run_mod <- function(iter){
    MCMCglmm(num_alleles_mean ~ logAbundance + SSD + BreedingType + bot + TPM80_ratio, # , #+ Abundance BreedingType  + BreedingType + Generation_time
        random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
        family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
        data=stats_mod_genlh,nitt=1100000,burnin=100000,thin=1000)
}

# check if model is saved
# model name
mod_name <- "gendiv_vs_lh_plus_bot"

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
# new model

# one model
mod_genlh <- models[[1]]

# summary
summary(mod_genlh)

# visually inspecting chain convergence
plot(mod_genlh$Sol)
plot(mod_genlh$VCV)
autocorr(mod_genlh$Sol)
autocorr(mod_genlh$VCV)

# variation explained by phylogeny
var_phy <- mod_genlh$VCV[, "tip_label"] / (mod_genlh$VCV[, "tip_label"] + mod_genlh$VCV[, "units"])
posterior.mode(var_phy)
median(var_phy)
HPDinterval(var_phy)

# commonality analyses and R2
model_file_name_R2 <- paste0(mod_name, "_R2.RData")

if (!file.exists(paste0("output/mcmcmodels/", model_file_name_R2))){
    set.seed(324)
    R2_genlh <- mcmcR2::partR2(mod_genlh, partvars = c("SSD", "BreedingType", "logAbundance", "bot", "TPM80_ratio"),
        data = stats_mod_genlh, inv_phylo = inv_phylo, prior = prior, 
        nitt = 1100000, burnin = 100000, thin = 1000)
    saveRDS( R2_genlh, file = paste0("output/mcmcmodels/", model_file_name_R2))
}

R2_genlh<- readr::read_rds(paste0("output/mcmcmodels/", model_file_name_R2))
R2_genlh
# out <- mcmcR2::R2mcmc(mod_hetexc)
# out$partR2
R2_genlh$R2 %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_R2" ,".txt"))
R2_genlh$SC %>% write_delim(paste0("output/mcmcmodels/", mod_name, "_SC" ,".txt"))

# save summary to file
mod_genlh %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_genlh$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_genlh$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI") %>% 
    write_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"))

beta_genlh<- readr::read_delim(paste0("output/mcmcmodels/", mod_name, "_beta" ,".txt"), delim = " ")


# plots ------------------------------------------------------------------------------------
point_size = 3.5
#


######## plot (2) ABCprob(bot) vs genetic diversity #########
stats_mod_gen <- all_stats
# re-run model, as we don't want standardized predictions for plotting
set.seed(101)
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
AR2_R2 <- read_delim("output/mcmcmodels/gen1_bot_R2.txt", delim = " ")
point_alpha <- 0.4

p1 <- ggplot(aes(x = bot, y = num_alleles_mean), data = all_stats) +
    geom_line(data = mod_preds_AR2, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    scale_y_continuous(breaks = seq(from = 2, to = 10, by = 2), limits = c(1.5,10), position = "right") +
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
    xlab(expression(ABC~bottleneck~model~probability~(p[bot]))) + #Allelic richness
    ylab(expression(Allelic~richness~(A[r]))) +
    #annotate("text", x = 0.24, y = 2, label = "R^2 == '0.49 [0.2, 0.7]'", 
    #    parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    #annotate("text", x = 0.24, y = 1.3, label = "beta == '-1.3 [-1.8, -0.78]'", 
    #    parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.9,0.5,0.25,0.1), "cm"),
        axis.line = element_line(colour = "#cccccc"),
        axis.ticks = element_line(colour = "#cccccc"),
        axis.title.y.right = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 15))) +
    geom_text_repel(label = all_stats$short,size = 2.5, alpha = 1, color = "grey50",#  aes(label = common) , 
        segment.alpha= 1,  box.padding = unit(0.4, "lines"), point.padding = unit(0.7, "lines"),
        segment.size = 0.1,  force = 1, min.segment.length = unit(0.01, "lines"))


p1

set.seed(102)
mod_div <- MCMCglmm(num_alleles_mean ~ logAbundance +  BreedingType + SSD + bot + TPM80_ratio, #
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=all_stats,nitt=110000,burnin=10000,thin=100)
summary(mod_div)

# prediction long model
pred_df <- data.frame(num_alleles_mean = 0, 
    logAbundance = rep(seq(from = 5, to = 15.5, by = 0.5), each = 2),
    BreedingType = c("land", "ice"), 
    SSD = mean(all_stats$SSD), 
    bot = mean(all_stats$bot),
    TPM80_ratio = mean(all_stats$TPM80_ratio),
    tip_label = all_stats$tip_label[1])

mod_preds <- data.frame(predict(mod_div, pred_df, interval = "confidence"))

mod_preds$logAbundance <- rep(seq(from = 5, to = 15.5, by = 0.5), each = 2)
mod_preds$BreedingType <- c("land", "ice")

# addition for long model


# per land and ice breeding
pred_df <- data.frame(num_alleles_mean = 0, logAbundance = rep(seq(from = 5, to = 15.5, by = 0.5), each = 2), 
    BreedingType = c("land", "ice"), SSD = mean(all_stats$SSD), bot = mean(all_stats$bot), 
    TPM80_ratio = mean(all_stats$TPM80_ratio), tip_label = all_stats$tip_label[1])

mod_preds <- data.frame(predict(mod_div, pred_df, interval = "confidence"))

mod_preds$logAbundance <- rep(seq(from = 5, to = 15.5, by = 0.5), each = 2)
mod_preds$BreedingType <- c("land", "ice")


p2 <- ggplot(aes(logAbundance, num_alleles_mean), data = all_stats) +
    geom_point(size = 3.5, alpha = 0.7, aes(color = BreedingType)) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, col ="black") +
    geom_line(data = mod_preds, aes(y = fit, color = BreedingType), size = 1, alpha = 0.5) +
    #geom_line(stat = "smooth", method = "lm",  alpha = 0.6,  aes(color = BreedingType)) +
    #geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.08, 
    #aes(fill = BreedingType)) +
    scale_color_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    # scale_fill_manual(values = c("cornflowerblue", "#d8b365"), name = "Breeding Habitat") +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    scale_x_continuous(trans = "log", breaks = c(log(100), log(1000), log(10000), log(100000), log(1000000), log(10000000)), 
        labels = c(expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}),  expression(10^{7})),
        limits = c(4.5, 16.2)) + 
    scale_y_continuous(breaks = c(seq(from = 0, to = 10, by =2)), limits = c(1.5, 10)) +
    theme(legend.position=c(0.3, 0.8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.9,0.1,0.396,0.1), "cm"),
        axis.line = element_line(colour = "#cccccc"),
        axis.ticks = element_line(colour = "#cccccc"),
        legend.title=element_text(size=10),
        axis.title.y = element_blank()) +
    xlab("Global abundance") +
    geom_text_repel(label = all_stats$short,size = 2.5, alpha = 1, color = "grey50",#  aes(label = common) , 
        segment.alpha= 1,  box.padding = unit(0.4, "lines"), point.padding = unit(0.7, "lines"),
        segment.size = 0.1,  force = 1, min.segment.length = unit(0.01, "lines"))

p2



# load model output
mod_beta <- read_delim("output/mcmcmodels/gendiv_vs_lh_plus_bot_beta.txt", delim = " ")
mod_R2 <- read_delim("output/mcmcmodels/gendiv_vs_lh_plus_bot_R2.txt", delim = " ")
mod_SC <- read_delim("output/mcmcmodels/gendiv_vs_lh_plus_bot_SC.txt", delim = " ")


# beta coefficients
mod_out <- mod_beta[-1, c("components", "post_median", "lower", "upper")]
names(mod_out) <- c("comps", "pe", "cilow", "cihigh")

mod_out$comps <- factor(mod_out$comps, levels = c("bot", "TPM80_ratio", "SSD", "BreedingTypeice", "logAbundance"))

p3 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3.5, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(1,0,0.15,0.5), "cm"),
        axis.title.x=element_text(margin=margin(t=0.35, unit = "cm"))) +
    #scale_y_discrete(labels = c("Breeding\nhabitat\nice vs. land", "log(Abundance)", 
    #    "Sexual Size\nDimorphism")) +
    xlab(expression(paste("Effect size ", beta))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p3


mod_out_SC <- mod_SC[, c("pred", "medianSC", "lower", "upper")]
names(mod_out_SC) <- c("comps", "pe", "cilow", "cihigh")

mod_out_SC$comps <- factor(mod_out_SC$comps, levels = c("bot", "TPM80_ratio", "SSD", "BreedingType", "logAbundance"))
# structure coefficients
p4 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_SC) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3.5, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    scale_x_continuous(breaks = c(-1,-0.5,0,0.5,1), limits = c(-1, 1)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = c(0.5)),
        plot.margin = unit(c(1,0.64,0.2,0.1), "cm")) +
    scale_y_discrete(labels = c(expression(p[bot]), expression(prop[het-exc]), "SSD",  "Breeding habitat", "Abundance"
        )) +
    xlab(expression(paste("Structure coefficient", " r(", hat(Y),",x)") )) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1)
p4

plot_grid(p3, p4)
p3 + p4
col_legend <- "#969696"

mod_out_R2 <- mod_R2[, c("combinations", "medianR2", "lower", "upper")]
names(mod_out_R2) <- c("comps", "pe", "cilow", "cihigh")
mod_out_R2
mod_out_R2 <- mod_out_R2[1:6, ]
#mod_out_R2$comps <- rev(fct_inorder(factor(mod_out_R2$comps)))
mod_out_R2$comps <- factor(mod_out_R2$comps, levels = rev(c("full model","logAbundance","BreedingType","SSD",
                                            "TPM80_ratio","bot")))

# structure coefficients
p5 <- ggplot(aes(pe, comps, xmax = cihigh, xmin = cilow), data = mod_out_R2 ) + 
    # geom_point(size = 3, color = "grey69") + # abc_out
    geom_errorbarh(alpha=0.4, color="black",height = 0) +
    geom_point(size = 3.5, shape = 21, col = "black", fill = "grey69") +
    # geom_errorbarh(alpha=0.4, color="black",height = 0) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = '#333333'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(hjust = c(0.5)),
        plot.margin = unit(c(0.35, 0.5, 0.3, 0.2), "cm")) +
    scale_y_discrete(labels = rev(c("Full model","Abundance","Breeding Habitat",
        "SSD", expression(prop[het-exc]), expression(p[bot])))) +
    xlab(expression(paste(R^{2}))) +
    geom_vline(xintercept = 0, color = "black", alpha = 0.1) +
    #annotate("segment", x = 0.9, xend = 0.9, y = 0.8, yend = 3.5, color = col_legend) +
    #annotate("text", x = 0.98, xend = 0.95, y = 2, yend = 3, color = col_legend, label = c("common"), angle = 270) +
    annotate("segment", x = 1, xend = 1, y = 0.7, yend = 5.1, color = col_legend) +
    annotate("text", x = 1.1, xend = 1.1, y = 3.2, yend = 3.5, color = col_legend, label = c("unique"), angle = 270, size = 3) +
    annotate("segment", x = 1, xend = 1, y = 5.5, yend = 6.6, color = col_legend) +
    annotate("text", x = 1.1, xend = 1.1, y = 5.7, yend = 6.2, color = col_legend, label = c("marginal"), angle = 270, size = 3)
p5



p_top <- plot_grid(p1, p2, rel_widths = c(1.1,1), labels = c("A", "B"),label_fontfamily = "Lato",
    label_x = 0.1, label_y = 0.98)
p_top 
p_bot <- plot_grid(p5, p3, p4, labels = c("C", "D", "E"),label_fontfamily = "Lato", rel_widths = c(1.5, 1, 1.4),
        nrow = 1)
p_bot

p_final <- plot_grid(p_top, p_bot, ncol = 1, rel_heights = c(1.5,1))
p_final

# jpg
ggsave('other_stuff/figures/figures_final/gendiv_lh_long.jpg',p_final,  width=8, height=6)

# pdf
ggsave('other_stuff/figures/figures_final/gendiv_lh_long.pdf',p_final,  width=8, height=6)
Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts("other_stuff/figures/figures_final/gendiv_lh_long.pdf")





# plot for marine mammal conference

p1 <- ggplot(aes(logharem_size, TPM80_ratio), data = all_stats) +
    geom_point(size = 3.5, alpha = 0.3) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, col = "black") + #, aes(fill = BreedingType)
    #scale_color_viridis(option = "magma", direction = -1,
    #    name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    theme_martin() +
    xlab("Harem Size") +
    ylab("Heterozygosity-excess") +
    scale_x_continuous(breaks = log(c(1,5,10,20,50)), labels = c(1,5,10,20,50)) + #breaks = c(1,2,3,4,5,6,7,8)
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 1, alpha = 0.5) + 
    #ylab("Heterozygosity-excess") +
    theme(#panel.grid.minor = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.major.x = element_blank(),
        plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm") ,
        #axis.line.x = element_line(color="#cccccc", size=0.3),
        #axis.ticks.x = element_line(color="#cccccc", size=0.3),
        #axis.text.x = element_text(margin = margin(t = 5)),
        #axis.title.y = element_text(margin = margin(r = 10))
        legend.direction = "vertical",
        legend.position = c(0.85,0.25),
        legend.title=element_text(size=10),
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))
    ) +
    scale_color_distiller(palette = "RdBu",
        direction = -1,
        name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    
    guides(color = guide_colorbar(barwidth = 0.5, barheight = 5, 
        title.position = "left")) + #, label.position = "bottom"
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) 
# geom_text_repel(aes(label = short),size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
#       segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
#    segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines"))
p1

ggsave("figures/SMM/p2_smm.jpg", width = 4, height = 3.5)



p2 <- ggplot(aes(logharem_size, TPM80_ratio), data = all_stats) +
    geom_point(size = 3.5, alpha = 0.3) + # abc_out
    geom_point(size = 3.5, alpha = 0.8, shape = 21, aes(fill = BreedingType)) + #, aes(fill = BreedingType)
    #scale_color_viridis(option = "magma", direction = -1,
    #    name = "ABC bottleneck \nprobability %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    theme_martin() +
    xlab("Harem Size") +
    ylab("Heterozygosity-excess") +
    scale_x_continuous(breaks = log(c(1,5,10,20,50)), labels = c(1,5,10,20,50)) + #breaks = c(1,2,3,4,5,6,7,8)
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 1, alpha = 0.5) + 
    #ylab("Heterozygosity-excess") +
    theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = unit(c(0.4,0.4,0.4,0.4), "cm") ,
        axis.line.x = element_line(color="#cccccc", size=0.3),
        axis.ticks.x = element_line(color="#cccccc", size=0.3),
        axis.text.x = element_text(margin = margin(t = 5)),
        #axis.title.y = element_text(margin = margin(r = 10))
        legend.direction = "vertical",
        legend.position = c(0.85,0.15),
        legend.title=element_text(size=10),
        axis.title.x=element_text(margin=margin(t=0.5, unit = "cm"))
    ) +
    scale_fill_manual(values = c("cornflowerblue", "#d8b365")) +
    guides(color = guide_colorbar(barwidth = 0.5, barheight = 5, 
        title.position = "left")) + #, label.position = "bottom"
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1, 1.13)) 
# geom_text_repel(aes(label = short),size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
#       segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
#    segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines"))
p2

ggsave("figures/SMM/p3_smm.jpg", width = 4, height = 3.5)


