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
                              Abundance, Generation_time) %>% 
                              mutate(Abundance = log(Abundance),
                                     harem_size = log(harem_size),
                                     male_weight = log(male_weight)) %>% 
                              rename(AlleleNumber = num_alleles_mean, TPM70 = TPM70_ratio,
                                     BreedSeasonLength= breeding_season_length , LogHaremSize =  harem_size,
                                     LifeSpan = life_span_years, GenTime = Generation_time,
                                     LogAbundance = Abundance, LogMaleWeight = male_weight,
                                     LactLength = lactation_length )

# gg duo makes a nice pairs plot, not used here
# stats_mod %>%                   
#     ggduo(types = list(continuous = "smooth_loess", 
#     comboVertical =  'box', 
#     comboHorizontal = "facethist",
#     discrete = "ratio")) +
#     # ggscatmat(columns = c("AlleleNumber", "TPM70", "BreedSeasonLength")) +
#     theme_minimal()



## modeling and plotting 

# figure 1 - Allelic Richness vs. Global Abundance-------------------------------------------------
all_stats_div <- all_stats %>% 
                    dplyr::select(common,
                        num_alleles_mean, 
                        Abundance,
                        bot,
                        common,
                        short,
                        SSD, 
                        Generation_time,
                        BreedingType,
                        TPM70_ratio
                        #mean_allele_range, 
                        #obs_het_mean, 
                        #prop_low_afs_mean, 
                        #  num_alleles_sd,sd_allele_range,prop_low_afs_sd,obs_het_sd,
                    ) %>%
                    mutate(Abundance = log(Abundance))

# count GS and HS as land breeding for this
both <- which(all_stats_div$short == "GS" | all_stats_div$short =="HS")
all_stats_div[both, "BreedingType"] <- "land"

CI_div <- function(offs) {
    model <- lm(num_alleles_mean ~ 
            I(Abundance-offs), # + Generation_time, # abc_out, # + SSD + Generation_time, + bot
             data =all_stats_div)
    ests <- summary(model)$coefficients[1,1:2]

    return(c(offs,ests,ests[1]+c(-1,0,1)*1.96*ests[2]))
}

offs_div <- seq(5,16,0.1)
res_div <- as.data.frame(t(sapply(offs_div, CI_div)))
names(res_div) <- c("Abundance", "Coefficient", "Std. Error", "Upper", "Mean", "Lower")

all_names <- sapply(all_stats_div$common, function(x) sub(" ", "\n", x ))

plot_col <- magma(40)
plot_col <- rev(plot_col[c(rep(FALSE,2), TRUE)]) # get every 5th element

# create legend for Abbreviations
all_abs <- paste(all_stats$short, "=", all_stats$common)
all_abs2 <- paste(all_abs, collapse = "\n")


p_AR_Abund <- all_stats_div  %>% 
    ggplot(aes(x = Abundance, y = num_alleles_mean)) +
    geom_ribbon(data = res_div, aes(x = Abundance, y = Mean, ymin = Lower, ymax = Upper), 
        alpha = 0.2, fill = "grey") +
    geom_line(data = res_div, aes(x = Abundance, y = Mean), size = 1.5, col = "grey", alpha = 0.5) +
    geom_point(size = 5, alpha = 0.6, aes(col = bot)) + # abc_out
    geom_point(size = 5, alpha = 0.8, shape = 21, col = "black") +
    # geom_smooth(method = "lm") +
    geom_text_repel(label = short,size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
        segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
        segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines")) +
    #scale_color_manual("ABC model", values = c("red","black"), labels = c("bottleneck", "constant")) +
    scale_color_gradientn(colours=plot_col, 
       guide = "colourbar",  name = "ABC \np(bot). %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    guides(color = guide_colorbar(barwidth = 6, barheight = 0.5, 
        title.position = "top", label.position = "bottom", override.aes = list(alpha = 0.5))) +
    xlab("log(Global Abundance)") +
    ylab("Allelic Richness") +
    scale_x_continuous(breaks = c(log(1000), log(10000), log(100000), log(1000000)), 
        labels = c(expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}))) +
    #scale_y_continuous(breaks = c(0.1,2,4,6,8,10)) + 
    scale_y_continuous(breaks = c(2,4,6,8,10), limits = c(1,10)) +
    #scale_x_continuous(trans = "log") +
    # facet_wrap(~variable, scales = "free") +
    # scale_y_continuous(trans = "log") +
    #theme_classic()+
    theme_minimal() +
    theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 10), 
        axis.title.y = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.3, colour = "grey40"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
        legend.position = c(0.2,0.95),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5,
        plot.margin = unit(c(2,1,1,0.5), "cm")
        ) 
p_AR_Abund

# for legend, create empty ggplot
df <- data.frame()
legend_p <- ggplot(df) + geom_point() + xlim(3,7) + ylim(0,10) + theme_nothing() +
            annotate("text", x = 5, y = 5, label = all_abs2, size = 2.7)
legend_p

p_AR <- plot_grid(p_AR_Abund, legend_p, nrow = 1, rel_widths = c(4,1.5))
p_AR
# ggsave(filename = "abund_vs_div2.jpg", plot = p_final, width = 8, height = 5, dpi = 300)


# Figure 2: Allelic Richness vs. BreedingType-------------------------------------------------------

p_AR_Breed <- all_stats_div   %>% 
    ggplot(aes(x = BreedingType, y = num_alleles_mean)) +
    #geom_tufteboxplot(alpha = 0.8, col = "black", fill = "grey") +
    # for rectangular annotations
    #geom_boxplot(alpha = 0.1, col = "grey", fill = "lightgrey", width = 0.4) +
    #annotate("rect", xmin =-Inf, xmax = Inf, ymin=0.5, ymax = Inf, alpha = 0.6, fill = "#451077FF") +
    #annotate("rect", xmin =-Inf, xmax = Inf, ymin=-Inf, ymax = 0.5, alpha = 0.6, fill = "#FCFDBFFF") +
    
    geom_boxplot(alpha = 0.4, col = "darkgrey", fill = "lightgrey", width = 0.4) +
    # geom_jitter(size = 5, alpha = 0.6) + # abc_out
    geom_jitter(size = 5, alpha = 0.6, shape = 21, col = "black", fill = "cornflowerblue",
        width = 0.08) +
    geom_text_repel(label = short,size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
        segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
        segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines")) +
    #scale_color_manual("ABC model", values = c("red","black"), labels = c("bottleneck", "constant")) +
    scale_color_gradientn(colours=plot_col, 
        guide = "colourbar",  name = "ABC \np(bot). %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    guides(color = guide_colorbar(barwidth = 6, barheight = 0.5, 
        title.position = "top", label.position = "bottom", override.aes = list(alpha = 0.5))) +
    xlab("Breeding Habitat") +
    ylab("Allelic Richness") +
    #ylim(c(1,10)) +
    # facet_wrap(~variable, scales = "free") +
     scale_y_continuous(breaks = c(2,4,6,8,10), limits = c(1,10)) +
   # scale_y_continuous(position = "right") +
    #theme_classic()+
    theme_minimal() +
    theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 12), 
        # axis.title.y = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 12),
       # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.3, colour = "grey40"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
        legend.position = c(0.2,0.9),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title.align = 0.5,
        plot.margin = unit(c(2,0.2,1,0), "cm"))
p_AR_Breed 

# p_AR_final <- plot_grid(p_AR_Abund, legend_p, p_AR_Breed , nrow = 1, rel_widths = c(7,2,5))

p_AR_final <- plot_grid(p_AR_Abund, p_AR_Breed , nrow = 1, rel_widths = c(7,5), 
                        labels = c("A", "B"), vjust = 3, hjust = -3)
p_AR_final


ggsave(filename = "diversity.jpg", plot = p_AR_final, width = 10, height = 5, dpi = 300)

# Figure 3: Het excess vs SSD ----------------------------------------------------------------------


CI_div <- function(offs) {
    model <- lm(TPM70_ratio ~ 
            I(SSD-offs), # + Generation_time, # abc_out, # + SSD + Generation_time, + bot
        data =all_stats_div)
    ests <- summary(model)$coefficients[1,1:2]
    
    return(c(offs,ests,ests[1]+c(-1,0,1)*1.96*ests[2]))
}

offs_div <- seq(0.5,7.5,0.1)
res_div <- as.data.frame(t(sapply(offs_div, CI_div)))
names(res_div) <- c("Abundance", "Coefficient", "Std. Error", "Upper", "Mean", "Lower")

all_names <- sapply(all_stats_div$common, function(x) sub(" ", "\n", x ))

plot_col <- magma(40)
plot_col <- rev(plot_col[c(rep(FALSE,2), TRUE)]) # get every 5th element

# create legend for Abbreviations
all_abs <- paste(all_stats$short, "=", all_stats$common)
all_abs2 <- paste(all_abs, collapse = "\n")


p_hetexc_SSD <- all_stats_div  %>% 
    ggplot(aes(x = SSD, y = TPM70_ratio)) +
    geom_ribbon(data = res_div, aes(x = Abundance, y = Mean, ymin = Lower, ymax = Upper), 
        alpha = 0.2, fill = "grey") +
    geom_line(data = res_div, aes(x = Abundance, y = Mean), size = 1.5, col = "grey", alpha = 0.5) +
    geom_point(size = 5, alpha = 0.6, aes(col = bot)) + # abc_out
    geom_point(size = 5, alpha = 0.8, shape = 21, col = "black") +
    # geom_smooth(method = "lm") +
    geom_text_repel(label = short,size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
        segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
        segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines")) +
    #scale_color_manual("ABC model", values = c("red","black"), labels = c("bottleneck", "constant")) +
    scale_color_gradientn(colours=plot_col, 
        guide = "colourbar",  name = "ABC \np(bot). %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    guides(color = guide_colorbar(barwidth = 6, barheight = 0.5, 
        title.position = "top", label.position = "bottom", override.aes = list(alpha = 0.5))) +
    xlab("Sexual Size Dimorphism (Male weight / Female weight)") +
    ylab("% loci with heterozygosity excess") +
   # scale_x_continuous(breaks = c(log(1000), log(10000), log(100000), log(1000000)), 
   #        labels = c(expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}))) +
    #scale_x_continuous(trans = "log") +
    # facet_wrap(~variable, scales = "free") +
   # scale_x_continuous(trans = "log") +
    # theme_classic()+
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1,1.2)) +
    theme_minimal() +
    theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 10), 
        axis.title.y = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.3, colour = "grey40"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
        legend.position = c(0.2,0.95),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5,
        plot.margin = unit(c(2,1,1,0.5), "cm")
    )  
p_hetexc_SSD

# Figure 4: Het excess vs BreedingType -------------------------------------------------------------
p_hetexc_Breed <- all_stats_div   %>% 
    ggplot(aes(x = BreedingType, y = TPM70_ratio)) +
    #geom_tufteboxplot(alpha = 0.8, col = "black", fill = "grey") +
    # for rectangular annotations
    #geom_boxplot(alpha = 0.1, col = "grey", fill = "lightgrey", width = 0.4) +
    #annotate("rect", xmin =-Inf, xmax = Inf, ymin=0.5, ymax = Inf, alpha = 0.6, fill = "#451077FF") +
    #annotate("rect", xmin =-Inf, xmax = Inf, ymin=-Inf, ymax = 0.5, alpha = 0.6, fill = "#FCFDBFFF") +
    
    geom_boxplot(alpha = 0.4, col = "darkgrey", fill = "lightgrey", width = 0.4) +
   # geom_jitter(size = 5, alpha = 0.6) + # abc_out
    geom_jitter(size = 5, alpha = 0.6, shape = 21, col = "black", fill = "cornflowerblue",
                width = 0.08) +
    geom_text_repel(label = short,size = 3, alpha = 0.7, color = "black", #  aes(label = common) , 
        segment.alpha= 0.2, box.padding = unit(0.7, "lines"), point.padding = unit(0.3, "lines"),
        segment.size = 0.3,  force = 3, min.segment.length = unit(0.1, "lines")) +
    #scale_color_manual("ABC model", values = c("red","black"), labels = c("bottleneck", "constant")) +
    scale_color_gradientn(colours=plot_col, 
        guide = "colourbar",  name = "ABC \np(bot). %", labels=c("0", "50", "100"), breaks = c(0.05,0.5,1)) +
    guides(color = guide_colorbar(barwidth = 6, barheight = 0.5, 
        title.position = "top", label.position = "bottom", override.aes = list(alpha = 0.5))) +
    xlab("Breeding Habitat") +
    ylab("% loci with heterozygosity excess") +
    ylim(c(0.1,1)) +
    # facet_wrap(~variable, scales = "free") +
    # scale_y_continuous(trans = "log") +
    scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), limits = c(0.1,1.2)) +
    # theme_classic()+
    theme_minimal() +
    theme(text = element_text(family="Arial"),
        axis.title.x = element_text(size = 12),
        axis.text.x  = element_text(size = 10), 
       # axis.title.y = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 12),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.3, colour = "grey40"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks.x = element_line(size = 0.3, colour = "grey40"),
        legend.position = c(0.2,0.95),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5,
        plot.margin = unit(c(2,1,1,0.5), "cm")
    ) 
    
p_hetexc_Breed 
p_hetexc <- plot_grid(p_hetexc_SSD, p_hetexc_Breed, nrow = 1, rel_widths = c(7,5), 
    labels = c("A", "B"), vjust = 3, hjust = -3)
p_hetexc

# plot_full <- plot_grid(p_AR_Abund, legend_p, p_AR_Breed,
#                        p_hetexc_SSD, NULL, p_hetexc_Breed,
#                        ncol = 3, rel_widths = c(5,2,4,5,2,5),
#                        rel_heights = c(1,1,1,1,2,1))


ggsave(filename = "het_exc.jpg", plot = p_hetexc, width = 10, height = 5, dpi = 300)



all_stats %>% 
    dplyr::select(common,
        num_alleles_mean, 
        Abundance,
        abc_out,
        common
        #mean_allele_range, 
        #obs_het_mean, 
        #prop_low_afs_mean, 
        #  num_alleles_sd,sd_allele_range,prop_low_afs_sd,obs_het_sd,
    ) %>%
    mutate(Abundance = log(Abundance)) %>% 
    # melt(id = c("common", "num_alleles_mean", "abc_out", "common")) %>% 
    ggplot(aes(x = Abundance, y = num_alleles_mean)) +
    geom_point(size = 5, alpha = 0.5, aes(col = abc_out)) +
    geom_point(size = 5, alpha = 0.8, shape = 21, col = "black") +
    geom_smooth(method = "lm", se = FALSE) +
    geom_text_repel(aes(label = common), size = 3, alpha = 0.5, color = "black") +
    # facet_wrap(~variable, scales = "free") +
    # scale_y_continuous(trans = "log") +
    theme_classic()

