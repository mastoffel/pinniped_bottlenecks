# plots for talk


library(pacman)
p_load(patchwork, ggtree, ape, phytools, dplyr, readxl, stringr, viridis, ggtree, ggthemr,
    reshape2, cowplot, ggthemes, ggimage, RColorBrewer, scales, forcats, readr, caper,
    yhat, dplyr, ggrepel, GGally, MCMCglmm, purrr, readr)
# phylogenetic comparative analysis
source("R/martin.R")

# load data and prepare mixed models

# load (modified) phylogeney. 26 species from 10ktrees plus 3 subspecies of ringed seal
tree_final <- read.tree("data/raw/phylogeny/30_species_10ktrees_final.tre")

# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_30_modeling.csv"))

# how many percent higher is Ar in ice than land-breeding seals?
Ars <- all_stats %>% 
    group_by(BreedingType) %>% 
    summarize(Ar = mean(num_alleles_mean))
Ars 
# phylogenetic mixed model preparation

# construct inverse phylo matrix and priors
inv_phylo <- inverseA(tree_final, nodes="TIPS",scale=FALSE)$Ainv #,scale=TRUE
prior<-list(G=list(G1=list(V=1,nu=0.002)),R=list(V=1,nu=0.002))

# 
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


##### plots
point_size = 3.5
#

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
AR2_R2 <- read_delim("output/mcmcmodels/gen1_bot_R2_marginal.txt", delim = " ")
point_alpha <- 0.4

set.seed(141)
p1 <- ggplot(aes(x = bot, y = num_alleles_mean), data = all_stats) +
    geom_line(data = mod_preds_AR2, aes(y = fit), size = 1, alpha = 0.5) +
    geom_point(size = point_size, alpha = point_alpha) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    scale_y_continuous(breaks = seq(from = 2, to = 10, by = 2), limits = c(1.5,10)) + #
    scale_x_continuous(breaks = seq(from = 0, to = 1, by = 0.2)) +
    xlab(expression(Bottleneck~model~probability)) + #Allelic richness
    ylab("Genetic diversity\n(Alleles per 10 individuals)") +
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
        axis.title.y.right = element_text(angle = 90, margin = margin(t = 0, r = 0, b = 0, l = 15))) 
   # geom_text_repel(label = all_stats$short,size = 2, alpha = 1, color = "grey50",#  aes(label = common) , 
    #    segment.alpha= 1,  box.padding = unit(0.4, "lines"), point.padding = unit(0.6, "lines"),
     #   segment.size = 0.1,  force = 3, min.segment.length = unit(0.2, "lines"))

p1

ggsave(filename = "../presentations/wednesday_seminar/bot_vs_ar.jpg", width = 3.5, height = 3)




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

set.seed(110)

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
    geom_text_repel(label = all_stats$short,size = 2, alpha = 1, color = "grey50",#  aes(label = common) , 
        segment.alpha= 1,  box.padding = unit(0.4, "lines"), point.padding = unit(0.3, "lines"),
        segment.size = 0.1,  force = 1, min.segment.length = unit(0.01, "lines"))

p2



