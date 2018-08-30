# Create Fig.1 from paper

# files needed:
# (1) seal_data_complete_rarefac10_30.csv (all pinniped data, including summary statistics)
# (2) sims_10000kbot500_*model_selection.txt (ABC model probabilities)
# (3) 30_species_10ktrees.tre (phylogenetic data)

# creates:
# all_stats_tree_30*.csv (all stats plus tree variable, for largest clusters (cl), or loci
# in HW (HW), or full datasets (""))

# plot phylogeny
library(ggtree)
library(ape)
library(phytools)
library(dplyr)
library(readxl)
library(stringr)
library(viridis)
library(ggtree)
library(reshape2)
library(cowplot)
library(ggthemes)
library(ggimage)
library(RColorBrewer)
library(scales)
library(forcats)
library(readr)
library(extrafont)
# ggthemr('dust')
source("R/martin.R")

calc_on_cluster <- FALSE
calc_on_HW <- FALSE
if (calc_on_cluster) {
    save_files <- "_cl"    
} else if (calc_on_HW) {
    save_files <- "_HW"
} else {
    save_files <- ""
}

# load all datasets
seals <- read_csv(paste0("data/processed/seal_data_complete_rarefac10_30", save_files, ".csv"))

# load model probablities from ABC
model_probs <- read_delim(paste0("data/processed/sims_10000kbot500", save_files, "_model_selection_30.txt"), 
    delim = " ", col_names = c("species", "bot", "neut"), skip = 1)
# modify species names in model_probs as they still contain _cl_
if (calc_on_cluster) {
    model_probs$species <- str_replace(model_probs$species, "_cl_[1-9]", "")
}
# are all names overlapping?
sum(seals$species %in% model_probs$species)

# join
seals <- left_join(seals, model_probs, by = "species")


# load higdon phylogeny
tree_final <- read.tree("data/raw/phylogeny/30_species_10ktrees_final.tre")
plot(tree_final)
tree_final$tip.label

# sort data like phylogeny -------------------------------------------------------------------------
# extract last name in latin
tree_lab <- unlist(lapply(strsplit(tree_final$tip.label, "_"), function(x) x[length(x)]))
tree_lab
stats_lab <- unlist(lapply(strsplit(seals$latin, " "), function(x) x[length(x)]))
# check which ones match # all (slight changes had to be done to the naming of bryonia and carconiphagus)
which(!(stats_lab %in% tree_lab)) # should be 0

# how to reorder the data to attach the tip.labels?
reorder_data_for_tiplab <- as.numeric(unlist(sapply(tree_lab, function(x) which(stats_lab == x))))
# bind tip labels
seals <- seals %>% 
    .[reorder_data_for_tiplab, ] %>% 
    bind_cols(tibble(tip_label = tree_final$tip.label), .) 

# create label including abbreviations
seals$common_abbr <- paste0(seals$common, " (", seals$short, ")")

# how is everything plotted
# if indexing tree_final$tip.label, this is the sequence how it is plotted
# plot_inds <- c(1:2, 10:12, 3:5, 9, 6:8, 13,14, 19,20, 18,17,15,16,27,28,25,26,21,22,23,24)
# 
# # indices for reordering stats according to the sequence in tree labels
# reorder_ind <- unlist(sapply(tree_lab, function(x) which(stats_lab == x)))
# 
# # plot inds needed within diversity plot later // this is the sequence in which the tree is plotted
# plot_inds <- c(1:2, 10:12, 3:5, 9, 6:8, 13,14, 19,20, 18,17,15,16,27,28,25,26,21,22,23,24)

ggtree(tree_final) + geom_tiplab() + xlim(0, 100)

# creat a data.frame where factor levels are in the ggtree plotting sequence
plotting_sequence <- c("weddell_seal", "leopard_seal" , "crabeater_seal","ross_seal", "ses", "nes",
    "hawaiian_monk_seal", "mediterranean_monk_seal", "lagoda_ringed_seal", "saimaa_ringed_seal",
    "baltic_ringed_seal", "arctic_ringed_seal",   "grey_seal_orkneys", "harbour_seal_waddensee",
    "hooded_seal", "bearded_seal", "galapagos_fur_seal", "south_american_fur_seal",
    "new_zealand_fur_seal","subantarctic_fur_seal", "antarctic_fur_seal", "guadalupe_fur_seal",
    "new_zealand_sea_lion", "south_american_sea_lion",
    "australian_fur_seal", "galapagos_sea_lion", "california_sea_lion",
    "stellers_sea_lion", "northern_fur_seal", "atlantic_walrus")


# fit everything to the sequence of plotting
plot_seq_df <- tibble("species" = plotting_sequence)
all_stats_tree <- left_join(plot_seq_df, seals, by = "species")
# make tip label first column
all_stats_tree <- all_stats_tree[c(2,1,3:ncol(all_stats_tree))]

plotting_sequence %in% all_stats_tree$species

# create factors and correct sequence
all_stats_tree <- all_stats_tree %>% 
    .[nrow(.):1, ] %>% 
    mutate(tip_label = fct_inorder(factor(tip_label)),
        species = fct_inorder(factor(species)),
        latin = fct_inorder(factor(latin)),
        common = fct_inorder(factor(common)),
        common_abbr = fct_inorder(factor(common_abbr)))

# write to file
if(!file.exists(paste0("data/processed/all_stats_tree_30", save_files, ".csv"))){
    write_excel_csv(all_stats_tree, paste0("data/processed/all_stats_tree_30", save_files, ".csv"))
}

all_stats_tree[all_stats_tree$BreedingType == "both", "BreedingType"] <- "land" 
all_stats_tree <- all_stats_tree %>% mutate(BreedingType = as.factor(as.character(BreedingType))) %>% data.frame()


# create data.frames for heatmaps -------------

# DIVERSITY
stand_div <- all_stats_tree %>% 
    dplyr::select(obs_het_mean, 
        num_alleles_mean) %>% 
    apply(2, scale) %>% 
    as_tibble() %>% 
    bind_cols(all_stats_tree[c("latin", "common", "species", "common_abbr")], .) 

# to check correct sequence just plot tree with node number
plotTree(tree_final)
tiplabels()
edgelabels()

ggtree(tree_final) + geom_text(aes(label=node))

# BOTTLENECK
bot_res <- all_stats_tree %>% 
    dplyr::select(latin, common, species, #IAM_ratio, 
        TPM70_ratio, TPM80_ratio, 
        TPM90_ratio, SMM_ratio) 

# ABC probs
abc_probs <- all_stats_tree %>% 
    dplyr::select(latin, common, species,
        bot, neut)

# group into families
#cls <- list(phocids=tree_final$tip.label[13:28],
#            otarids= tree_final$tip.label[2:12],
#            odobenids = tree_final$tip.label[1])
#tree_final <- groupOTU(tree_final, cls)

#tree_final <- groupClade(tree_final, node = 31)
#tree_final <- groupClade(tree_final, node = 1)
#img <- c("phylo_figures/AFS.jpg")
#names(img) <- "31"

# prepare data frame
all_stats_for_tree <- all_stats_tree
# all_stats_for_tree$IUCN_rating[is.na(all_stats_for_tree$IUCN_rating)] <- "data deficient"
all_stats_for_tree$IUCN_rating <- factor(all_stats_for_tree$IUCN_rating, levels = c("least concern", "near threatened", "vulnerable", "endangered"))
# check whther necessary
names(all_stats_for_tree)[1] <- "taxa"

library(magrittr)
# all_stats_for_tree %<>% mutate(color_breed = ifelse(BreedingType == "land", "#737373", "cornflowerblue"))
#test <- all_stats_for_tree %>% 
#    mutate(abund_classes = ifelse(Abundance > 1000000, "1000k", )

## node coloring problem
# convert the 1st column to 'node' number and stored in an additional column called 'node' will solve this issue.

# all_stats_for_tree$node <- NA
# tipnode <- seq_along(tree_final$tip.label)
# names(tipnode) <- tree_final$tip.label
# all_stats_for_tree$node <- tipnode[all_stats_for_tree$taxa] ## convert the tip label to tip node number
# i <- is.na(all_stats_for_tree$node)
# all_stats_for_tree$node[i] =all_stats_for_tree$taxa[i] ## fill in internal node number


## alternative ----------------------------------
# add data


# produce new data.frame with nrow == edge number for plotting in ggtree

stats_for_tree_resorted <- data.frame(taxa = tree_final$tip.label)
stats_df <- left_join(stats_for_tree_resorted, all_stats_for_tree, by = "taxa")

new_df <- data.frame(matrix(ncol = ncol(stats_df), nrow = 27))
names(new_df) <- names(stats_df)
tree_df <- rbind(stats_df, new_df)
tree_df <- data.frame(node = 1:57, tree_df)

test_df <- mutate(tree_df, 
    color = case_when(is.na(BreedingType) ~ "#bdbdbd",
        BreedingType == "ice"~ "cornflowerblue",
        BreedingType == "land" ~ "#d8b365"))

##### former ######
tree_final$edge
# create tree view
# p <- ggtree(tree_final, color = "#737373") #

table(all_stats_for_tree$BreedingType)

##### newer
p <- ggtree(tree_final)
p <- p %<+% test_df
# p + geom_tippoint() #aes(color=IUCN_rating)

p <- p +  #layout="circular" , "fan"open.angle=180  #,color = "#737373"
    aes(color = I(color)) +
    # scale_color_manual(values=c("black", "#9e9ac8","#6a51a3", "#bcbddc" )) + #color = "#737373" , color = "#737373"
    geom_tippoint(aes(size = Abundance, fill=IUCN_rating, color = color),stroke=0.5, shape=21) + #shape=21,stroke=0.5 #color = "#737373" color = "#737373"
    scale_color_manual(name = "Breeding habitat", values = c("#bdbdbd", "#d8b365", "cornflowerblue"),
        breaks = c("cornflowerblue","#d8b365" ), labels = c("ice", "land")) +
    # geom_point(aes(color = abc, size = Abundance), shape=21) +
    # geom_treescale() + 
    # ggtitle("Pinniped Phylogeny") +
    #geom_tiplab(size=5, color="black") + # tiplap 2 for circular
    #ggplot2::xlim(0, 15) +
    # scale_color_manual(values = c("#66c2a5", "#fdae61", "#f46d43", "#d53e4f", "#000000")) +
    # scale_color_manual(values = c("lightblue", "goldenrod", "red", "darkred", "grey")) +
    # scale_fill_manual(values = c("#fee8c8", "#fc8d59", "#d7301f", "#7f0000", "white")) +
    # scale_fill_manual(values = c("#FDE4A6FF", "#FB8861FF", "#B63679FF", "#000004FF", "white")) +
    # scale_fill_manual(values = c("#ffffd4", "#fed98e", "#fe9929", "#cc4c02", "white")) +
    #scale_fill_manual(values = c("#f7f7f7", "#cccccc", "#525252", "#000000", "white")) +
#scale_fill_manual(values = c("#e0ecf4", "#9ebcda", "#8c6bb1", "#4d004b", "white")) +
# scale_color_manual(values = c( "blue", "black")) +
scale_fill_manual(values = c("#f7f7f7", "#d9d9d9", "#969696", "#252525", "white")) +
    # for coloring branches
    # scale_color_manual(na.value = "#969696", values = c("#cb181d", "#969696")) +
    #scale_color_manual(values = c("#c7e9b4", "#41b6c4", "#225ea8", "#081d58", "lightgrey")) +
    scale_size_continuous(range = c(1.5, 6), trans = "sqrt", breaks=c(1000, 10000, 100000, 1000000),
        labels = c(expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}))
    ) + #range = c(0.1,5), , breaks=c(500, 1000, 10000, 100000, 1000000)
    guides(fill = guide_legend(title = "IUCN rating", title.position = "top", direction = "vertical", order = 2, override.aes = list(size=2.5)),
        size = guide_legend(title.position = "top", title = "Global abundance", direction = "horizontal", order = 1),
        color = guide_legend(title.position = "top", direction = "horizontal", order = 3)) + #color = guide_legend(title = "supported model by ABC", direction = "horizontal",label = c("bot", "const"))
    theme(plot.margin=unit(c(72, -5,32,10),"points"), #c(30,-100,20,0) unit(c(50,-50,20,0) #c(52, -5,10,10)
        legend.position= c(0.24,0.93), #legend.direction = "horizontal",
        legend.spacing = unit(5, "points"),
        legend.key.height=unit(1,"line"),
        legend.key.width = unit(1, "line"), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5,
        text=element_text(family='Arial')) #Arial
# geom_image(data = d_img, aes(image = image), size = 0.2)
p


# p1 <- flip(p,  42, 49)

# diversity
stand_div_lf <- stand_div %>% melt(id.vars = c("species","common","common_abbr"), 
    measure.vars = c("obs_het_mean", "num_alleles_mean"))
str(stand_div_lf$variable)
stand_div_lf$variable <- relevel(stand_div_lf$variable, ref = "num_alleles_mean")


#pal <- colorRampPalette(c("cornflowerblue", "white", "darkgrey"))
pal <- colorRampPalette(c("darkblue", "#f0f0f0", "#4d4d4d"))
pal <- colorRampPalette(c("#2166ac", "#ece7f2", "#525252")) ##d0d1e6 #"#f0f0f0"
pal <- colorRampPalette(c("#016c59", "#ece7f2", "#525252"))

p_div <- ggplot(stand_div_lf, aes(x = variable, y = common_abbr, fill = value)) + 
    geom_tile(color = "grey", size = 0.1) +

    scale_fill_distiller(name = "Standardized \ngenetic diversity", palette = "RdYlBu", direction = 1) +
    #scale_fill_distiller(pal(5)) +
    # scale_fill_gradientn(colours=plot_col_div, 
    #   name = "prop. \nhet-exc") +
    theme_tufte(base_family="Arial") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0, size = 12, margin = margin(t = 7)),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(15,10,20.7,10),"points"), #c(38,-400,7,-260) c(5,-150,7,-260) #15 10 22 10
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_text(hjust=0, size = 9, colour = "#525252"), #abc_cols[stand_div$abc][plot_inds]
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5, 
        text=element_text(family='Arial')) +
    # coord_fixed(ratio = 0.7) +
    scale_x_discrete(labels=c(expression(A[r]), expression(H[o])), # c("Ho", "Ar"),  ##"ARA",
        position = "bottom") +
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, 
        title.position = "top", label.position = "bottom")) 

# grid.arrange(p, p_div, nrow = 1)
plot_grid(p, p_div, nrow = 1)

# het-excess and mratio

#pal2 <- colorRampPalette(c("darkred", "#f0f0f0", "#737373"))
#pal2 <- colorRampPalette(c("#d8b365", "#f5f5f5", "#5ab4ac"))
#pal2 <- colorRampPalette(c("#ef8a62", "#f7f7f7", "#67a9cf"))
pal2 <- colorRampPalette(c("#2166ac", "#ece7f2", "#525252")) ##d0d1e6 #"#f0f0f0"
#pal2 <- colorRampPalette(c("#f0f0f0", "#737373"))
bot_res_lf <- bot_res %>% melt(id.vars = c("species", "common", "latin")) %>% 
                filter(variable == "TPM80_ratio")

p_bot <- ggplot(bot_res_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "grey", size = 0.1) +
    scale_fill_distiller(palette="RdBu",
        # name = "Prop. of loci with\n   heterozyosity-\n              excess", 
        labels=c(0, 0.5, 1), breaks = c(0,0.5,1.01),
        direction = -1, limits = c(0,1.01)) +
    # scale_fill_gradientn(colours=rev(pal2(5)),
    #    name = "% of loci with\nheterozyosity excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,0.95)) +
    guides(fill = guide_colorbar(barwidth = 5.3, barheight = 0.5, 
        title.position = "top", direction = "horizontal",
        title = "Prop. of loci with\nheterozyg.-excess", 
        #title = expression(atop("Prop. loci with", paste("heterozygosity-excess"))),
        size = 4, title.hjust = c(0.5)
        #title.theme = element_text(vjust = 2, angle = 0)
        )) +
        #title = "Prop.\nHet.-excess")) +
    theme_tufte(base_family="Arial") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.4, size = 12, margin = margin(t = 7)),
        # axis.text.x = element_blank(),
        legend.position=c(0.57, 1.208),
        # legend.position=c(0.57, 1.2337),
       # legend.margin=margin(0,0,0,0),
        plot.margin=unit(c(80.8,10,22.4,10),"points"), 
        #   plot.margin=unit(c(87.8,10,22.8,10),"points"),
        #    plot.margin=unit(c(77.8,10,20.1,10),"points"),
        #plot.margin=unit(c(8.3,10,19.5,10),"points"),# c(38,-770,7, -690) c(38,300,3, 0)
        #plot.margin=unit(c(20,10,23,10),"points"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
       # legend.title.align = -0.5, 
        text=element_text(family='Arial'),
        
        legend.justification = "top") +
    # scale_x_discrete(labels=c("PHE","M")) +
    scale_x_discrete(labels= expression(prop[het-exc]))  # "IAM",
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
 #, draw.llim = 0, draw.ulim=1 title.hjust = -2 title = "Prop. of loci with\nheterozyosity-\nexcess"
p_bot
plot_grid(p, p_div, p_bot, nrow = 1, rel_widths = c(1,0.7, 0.3))


#pal3 <- colorRampPalette(c("darkgreen", "#f0f0f0", "#737373"))
#pal3 <- colorRampPalette(c("#ef8a62", "#f7f7f7", "#67a9cf"))
pal3 <- colorRampPalette(c("#525252","#bdbdbd", "#f0f0f0"))
abc_probs_lf <- abc_probs %>% melt(id.vars = c("species", "common", "latin"))

p_abc <- ggplot(abc_probs_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "grey", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    #scale_fill_viridis(option = "magma", direction = -1,
    #    name = "ABC \nprob. %", labels=c(0, 50, 100), breaks = c(0,0.5,1)) +
    #scale_fill_gradientn(colours=rev(pal3(5)), 
    #    name = "ABC \nprob. %", labels=c(0, 50, 100), breaks = c(0,0.5,1)) +
    scale_fill_distiller(palette = "RdBu",
        name = "ABC model \nprobability %", labels=c(0, 50, 100), breaks = c(0,0.5,1),
        direction = -1) +
    guides(fill = guide_colorbar(barwidth = 4.6, barheight = 0.5, 
        title.position = "top")) +
    theme_tufte(base_family="Arial") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.4, size = 12, margin = margin(t = 7)),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(14.71,10,20.4,10),"points"), #c(38,-520,3, -660)c(38,-520,1, -660) +20
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5,
        text=element_text(family='Arial')) +
    scale_x_discrete(labels=c(expression(p[bot]),expression(p[non-bot]))) 
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
  

p_final <- plot_grid(p, p_div, p_bot, p_abc, nrow = 1, 
       rel_widths = c(0.3, 0.28, 0.07, 0.12))
#p_final <- p_final +
#    draw_plot_label(c("A", "B", "C"), c(0.62, 0.74, 0.87), c(0.995, 0.995, 0.995), size = 11)
p_final <- p_final +
    draw_plot_label(c("A", "B", "C"), c(0.67, 0.79, 0.91), c(0.995, 0.995, 0.995), size = 11)


ggplot2::ggsave(paste0("other_stuff/figures/figures_final/figures_final_editing/fig1_phylo_plot_color_30", save_files, ".jpg"), p_final,
    height = 6.1, width = 10, device = "jpg")


ggplot2::ggsave(paste0("other_stuff/figures/figures_final/figures_final_editing/fig1_phylo_plot_color_30", save_files, ".pdf"), p_final,
    height = 6.1, width = 10, device = "pdf")

Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts(paste0("other_stuff/figures/figures_final/figures_final_editing/fig1_phylo_plot_color_30", save_files, ".pdf"))




# favorites Nunito Average2San Hind2Guntur Hind2Guntur Light
# 












