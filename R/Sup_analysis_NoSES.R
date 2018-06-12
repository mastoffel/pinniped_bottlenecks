# phylogenetic comparative analysis
# packages
library(pacman)
p_load(tibble, dplyr, ggtree, ape, phytools, readxl, stringr, viridis, ggtree, ggthemr, cowplot,
    ggthemes, ggimage, RColorBrewer, scales, forcats, readr, caper, yhat, dplyr, ggrepel,
    GGally, tidyr, kinship2, MCMCglmm, purrr, readr)
source("R/martin.R")
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

stats_mod_hetexc <- all_stats %>% 
    mutate(Abundance = ((Abundance - mean(Abundance)) / (2*sd(Abundance))), 
        Generation_time = (Generation_time - mean(Generation_time) / (2*sd(Generation_time))),
        SSD = (SSD - mean(SSD) / (2*sd(SSD))),
        logharem_size = (log(harem_size)),
        life_span_years = (life_span_years - mean(life_span_years)) / (2*sd(life_span_years)),
        breeding_season_length = (breeding_season_length - mean(breeding_season_length)) / (2*sd(breeding_season_length))) %>% 
    #log_breeding_season_length = log(breeding_season_length),
    #log_breeding_season_length = (log_breeding_season_length - mean(log_breeding_season_length)) / (2*sd(log_breeding_season_length))) %>% 
    mutate(logharem_size = (logharem_size - mean(logharem_size) / (2*sd(logharem_size)))) %>% 
    mutate(harem_size = (harem_size - mean(harem_size) / (2*sd(harem_size)))) 

stats_mod_hetexc <- stats_mod_hetexc %>% mutate(BreedingType = as.factor(BreedingType)) %>% 
    mutate(BreedingType = relevel(BreedingType, ref = "land"))


mod_hetexc_plot <- MCMCglmm(TPM80_ratio ~ SSD, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_hetexc,nitt=nitt,burnin=burnin,thin=thin)

# R2 0.2966127 CI 0.0003757109 0.6370375
R2_hetexc <- mcmcR2::partR2(mod_hetexc_plot, type = "marginal",
    partvars = c("SSD"),
    data = stats_mod_hetexc, inv_phylo = inv_phylo, prior = prior, 
    nitt = nitt, burnin = burnin, thin = thin)

# slope: 0.06924453 CI[0,002, 0.13]
beta_hetexc <- mod_hetexc_plot %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_hetexc_plot$Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_hetexc_plot$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI")

pred_df_hetexc <- data.frame(TPM80_ratio = 0,
    SSD = seq(from = 0.5, to = 8, by = 0.1),
    tip_label = all_stats$tip_label[1])

mod_preds_hetexc <- data.frame(predict(mod_hetexc_plot, pred_df_hetexc, interval = "confidence")) %>% 
    mutate(SSD = seq(from = 0.5, to = 8, by = 0.1))
point_size = 3.5
set.seed(57)
p1 <- ggplot(aes(x = SSD, y = TPM80_ratio), data = all_stats) +
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 0.2, alpha = 0.5) +
    geom_point(size = point_size, alpha = 0.7,  aes(col = bot)) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    xlab("Sexual Size Dimorphism (SSD)") +
    ylab(expression(Heterozygosity-excess ~ "("~prop[het-exc]~")")) +
    #ylab("Heterozygosity-excess") +
    # scale_x_continuous(breaks = log(c(1,2,3,4,6,8,10,15,20,30,40,50)), labels = c(1,2,3,4,6,8,10,15,20,30,40,50)) + #breaks = c(1,2,3,4,5,6,7,8)
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 1.2, alpha = 0.5, color = "grey") + 
    scale_fill_manual(values = c("cornflowerblue", "goldenrod")) +
    annotate("text", x = 6, y = 0.15, label = "R^2 == '0.30 CI [0.00, 0.64]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
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
p1






# exclude SES
stats_mod_hetexc_noSES <- filter(stats_mod_hetexc, species != "ses")
all_stats_noSES <- filter(all_stats, species != "ses")

mod_hetexc_plot_noSES <- MCMCglmm(TPM80_ratio ~ SSD, # , #+ Abundance BreedingType  + BreedingType + Generation_time
    random=~tip_label, nodes = "TIPS", #   rcov =~us(trait):units
    family=c("gaussian"),ginverse=list(tip_label=inv_phylo),prior=prior,
    data=stats_mod_hetexc_noSES,nitt=110000,burnin=10000,thin=100)

# R2 SSD 0.2562654 CI 0.000951954881 0.7870298
R2_hetexc_noSES <- mcmcR2::partR2(mod_hetexc_plot_noSES, type = "marginal",
    partvars = c("SSD"),
    data = stats_mod_hetexc_noSES, inv_phylo = inv_phylo, prior = prior, 
    nitt = nitt, burnin = burnin, thin = thin)

# slope: 0.07740149  CI[-0.02, 0.17]
beta_hetexc_noSES <- mod_hetexc_plot_noSES  %>% 
    summary() %$%
    solutions %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("components") %>% 
    mutate(post_median = apply(mod_hetexc_plot_noSES $Sol, 2, median)) %>% 
    mutate(post_mode = posterior.mode(mod_hetexc_plot_noSES$Sol)) %>% 
    .[c(1,2,7,8,3:6)] %>% 
    rename(post_mean= post.mean,
        lower =  "l-95% CI",
        upper = "u-95% CI")



pred_df_hetexc_noSES <- data.frame(TPM80_ratio = 0,
    SSD = seq(from = 0.5, to = 8, by = 0.1),
    tip_label = all_stats$tip_label[1])

mod_preds_hetexc <- data.frame(predict(mod_hetexc_plot_noSES, pred_df_hetexc_noSES, interval = "confidence")) %>% 
    mutate(SSD = seq(from = 0.5, to = 8, by = 0.1))
point_size = 3.5
set.seed(57)

p2 <- ggplot(aes(x = SSD, y = TPM80_ratio), data = all_stats_noSES) +
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 0.2, alpha = 0.5) +
    geom_point(size = point_size, alpha = 0.7,  aes(col = bot)) + # abc_out
    geom_point(size = point_size, alpha = 0.8, shape = 21, col = "black") +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    xlab("Sexual Size Dimorphism (SSD)") +
    ylab(expression(Heterozygosity-excess ~ "("~prop[het-exc]~")")) +
    #ylab("Heterozygosity-excess") +
    # scale_x_continuous(breaks = log(c(1,2,3,4,6,8,10,15,20,30,40,50)), labels = c(1,2,3,4,6,8,10,15,20,30,40,50)) + #breaks = c(1,2,3,4,5,6,7,8)
    scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10)) +
    geom_line(data = mod_preds_hetexc, aes(y = fit), size = 1.2, alpha = 0.5, color = "grey") + 
    scale_fill_manual(values = c("cornflowerblue", "goldenrod")) +
    #ylab("Heterozygosity-excess") +
    annotate("text", x = 6, y = 0.15, label = "R^2 == '0.26 CI [0.00, 0.71]'", 
        parse = TRUE, family = "Lato", size = 3.1, colour = "#333333") +
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
    geom_text_repel(label = all_stats_noSES$short,size = 2.5, alpha = 1, color = "grey50",#  aes(label = common) , 
        segment.alpha= 1,  box.padding = unit(0.4, "lines"), point.padding = unit(0.7, "lines"),
        segment.size = 0.1,  force = 1, min.segment.length = unit(0.01, "lines"))
p2


library(patchwork)

p_final <- p1 + p2

ggsave('other_stuff/figures/figures_final/Sup_NoSES.jpg',p_final,  width=8.5, height=3.5)






