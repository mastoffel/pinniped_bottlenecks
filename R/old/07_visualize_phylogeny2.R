# Create Fig.1 from paper

# files needed:
# (1) seal_data_complete_rarefac10_29.csv (all pinniped data)
# (2) sims_10000k_model_selection.txt (ABC model probabilities)
# (3) 29_species_10ktrees.tre (phylogenetic data)

# creates:
# all_stats_tree_29.csv (all stats plus tree variable)

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
if (calc_on_cluster) {
    save_files <- "_cl"    
} else {
    save_files <- ""
}

# load all datasets
seals <- read_csv(paste0("data/processed/seal_data_complete_rarefac10_29", save_files, ".csv"))

# load model probablities from ABC
model_probs <- read_delim(paste0("data/processed/sims_10000k", save_files, "_model_selection.txt"), 
    delim = " ", col_names = c("species", "bot", "neut"), skip = 1)
# modify species names in model_probs as they still contain _cl_
if (calc_on_cluster) {
    model_probs$species <- str_replace(model_probs$species, "_cl_[1-9]", "")
}
# are all names overlapping?
sum(seals$species %in% model_probs$species)

# join
seals <- left_join(seals, model_probs, by = "species")

#### stopped here
# load higdon phylogeny
tree_final <- read.tree("data/raw/phylogeny/29_species_10ktrees.tre")
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
    "new_zealand_fur_seal","subantarctic_fur_seal", "antarctic_fur_seal", "new_zealand_sea_lion", "south_american_sea_lion",
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
# if(!file.exists(paste0("data/processed/all_stats_tree_29", save_files, ".csv"))){
#     write_excel_csv(all_stats_tree, paste0("data/processed/all_stats_tree_29", save_files, ".csv"))
# }

all_stats_tree[all_stats_tree$BreedingType == "both", "BreedingType"] <- "land" 
all_stats_tree <- all_stats_tree %>% mutate(BreedingType = as.factor(as.character(BreedingType))) %>% data.frame()


# ggtree(tr, branch.length='none', layout='circular')

# prepare data frame
all_stats_for_tree <- all_stats_tree
# all_stats_for_tree$IUCN_rating[is.na(all_stats_for_tree$IUCN_rating)] <- "data deficient"
all_stats_for_tree$IUCN_rating <- factor(all_stats_for_tree$IUCN_rating, levels = c("least concern", "near threatened", "vulnerable", "endangered"))
# check whther necessary
names(all_stats_for_tree)[1] <- "taxa"

# produce new data.frame with nrow == edge number for plotting in ggtree

stats_for_tree_resorted <- data.frame(taxa = tree_final$tip.label)
stats_df <- left_join(stats_for_tree_resorted, all_stats_for_tree, by = "taxa")

new_df <- data.frame(matrix(ncol = ncol(stats_df), nrow = 27))
names(new_df) <- names(stats_df)
tree_df <- rbind(stats_df, new_df)
tree_df <- data.frame(node = 1:56, tree_df)

test_df <- mutate(tree_df, 
    color = case_when(is.na(BreedingType) ~ "#bdbdbd",
        BreedingType == "ice"~ "cornflowerblue",
        BreedingType == "land" ~ "#d8b365"))

##### former ######
##### newer
# p <- ggtree(tree_final, branch.length='none', layout='circular')
plot_heat <- tree_df[c(79, 81, 83)]

p <- ggtree(tree_final, layout = "fan", open.angle = 183)
p <- p %<+% test_df
# p + geom_tippoint() #aes(color=IUCN_rating)
p <- p +  #layout="circular" , "fan"open.angle=180  #,color = "#737373"
    aes(color = I(color)) +
    geom_tiplab2(aes(label = short), size = 3, offset = 10) +
    # scale_color_manual(values=c("black", "#9e9ac8","#6a51a3", "#bcbddc" )) + #color = "#737373" , color = "#737373"
    geom_tippoint(aes(size = Abundance,  color = color),stroke=0.5, shape=21) + #fill=IUCN_rating, shape=21,stroke=0.5 #color = "#737373" color = "#737373"
    scale_color_manual(name = "Breeding habitat", values = c("#bdbdbd", "#d8b365", "cornflowerblue"),
        breaks = c("cornflowerblue","#d8b365" ), labels = c("ice", "land")) +
    scale_fill_manual(values = c("#f7f7f7", "#d9d9d9", "#969696", "#252525", "white")) +
    scale_size_continuous(range = c(1.5, 6), trans = "sqrt", breaks=c(1000, 10000, 100000, 1000000),
        labels = c(expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}))
    ) + #range = c(0.1,5), , breaks=c(500, 1000, 10000, 100000, 1000000)
    guides(fill = guide_legend(title = "IUCN rating", title.position = "top", direction = "vertical", order = 2),
        size = guide_legend(title.position = "top", title = "Global abundance", direction = "horizontal", order = 1),
        color = guide_legend(title.position = "top", direction = "horizontal", order = 3)) + #color = guide_legend(title = "supported model by ABC", direction = "horizontal",label = c("bot", "const"))
    theme(plot.margin=unit(c(68, -5,30,10),"points"), #c(30,-100,20,0) unit(c(50,-50,20,0) #c(52, -5,10,10)
        legend.position= c(0.01,0.82), #legend.direction = "horizontal",
        legend.spacing = unit(5, "points"),
        legend.key.height=unit(1,"line"),
        legend.key.width = unit(1, "line"), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5,
        text=element_text(family='Hind Guntur Light')) +
    ggplot2::xlim(0, 50)
p
str(stats_df)
# important!
rownames(stats_df) <- stats_df$taxa

# prep
stand_div <- all_stats_tree %>% 
    dplyr::select(obs_het_mean, 
        num_alleles_mean) %>% 
    apply(2, scale) %>% 
    as_tibble() %>% 
   bind_cols(all_stats_tree[c("latin")], .) 
rownames(stand_div) <- all_stats_tree$tip_label

# BOTTLENECK
bot_res <- all_stats_tree %>% 
    dplyr::select(latin, common, species, #IAM_ratio, 
        TPM70_ratio, TPM80_ratio, 
        TPM90_ratio, SMM_ratio) 
rownames(bot_res) <- all_stats_tree$tip_label

# ABC probs
abc_probs <- all_stats_tree %>% 
    dplyr::select(latin, common, species,
        bot, neut)

col <- c(brewer.pal(7, 'RdYlBu'),
         brewer.pal(7, 'RdBu'))

col <- c(brewer.pal(9, "YlGnBu"))
# heatmap
p2 <- gheatmap(p, data = stand_div[, c(2, 3), drop = FALSE], width = 0.05,
    offset = 3, colnames_angle=90, font.size=1.5)  
    # scale_fill_distiller(name = "Standardized \ngenetic diversity", palette = "RdYlBu", direction = 1) 
p2

p3 <- gheatmap(p2, data = bot_res[, "TPM80_ratio", drop = FALSE], width = 0.05,
    offset = 5, colnames_angle=90, font.size=1.5) +
    scale_color_gradientn(colours = col)
    # scale_fill_distiller(name = "Standardized \ngenetic diversity", palette = "RdYlBu", direction = 1)
p3
# geom_image(data = d_img, aes(image = image), size = 0.2)
p


p <- ggtree(tree_final, layout = "fan", open.angle = 183)

p <- p %<+% stats_df
p + geom_tiplab2(aes(label = short))
p

col <- c("red", "blue")
gheatmap(p, stats_df[, "BreedingType", drop = FALSE], offset = 0.8, width=0.1) +
    scale_fill_manual(values = col) 















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
tree_df <- data.frame(node = 1:56, tree_df)

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
    guides(fill = guide_legend(title = "IUCN rating", title.position = "top", direction = "vertical", order = 2),
        size = guide_legend(title.position = "top", title = "Global abundance", direction = "horizontal", order = 1),
        color = guide_legend(title.position = "top", direction = "horizontal", order = 3)) + #color = guide_legend(title = "supported model by ABC", direction = "horizontal",label = c("bot", "const"))
    theme(plot.margin=unit(c(68, -5,30,10),"points"), #c(30,-100,20,0) unit(c(50,-50,20,0) #c(52, -5,10,10)
        legend.position= c(0.24,0.92), #legend.direction = "horizontal",
        legend.spacing = unit(5, "points"),
        legend.key.height=unit(1,"line"),
        legend.key.width = unit(1, "line"), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5,
        text=element_text(family='Hind Guntur Light')) 
# geom_image(data = d_img, aes(image = image), size = 0.2)
p


# p1 <- flip(p,  42, 49)

# diversity
stand_div_lf <- stand_div %>% melt(id.vars = c("species","common","common_abbr"), 
    measure.vars = c("obs_het_mean", "num_alleles_mean"))

#plot_col_div <- viridis(20)
#plot_col_div <- rev(plot_col_div[c(rep(FALSE,4), TRUE)]) # get every 5th element

#plot_col_div <- viridis(100)
#plot_col_div <- rev(plot_col_div[c(rep(FALSE,30), TRUE)])

#color names
#numColors <- length(levels(stand_div_lf$abc))
#getColors <- brewer_pal('qual')
#myPalette <- getColors(numColors)
# abc_cols <- c("#a50f15", "#525252")
# 
# plot_col <- magma(40)
# plot_col <- plot_col[c(rep(FALSE,2), TRUE)] # get every 5th element

#pal <- colorRampPalette(c("cornflowerblue", "white", "darkgrey"))
pal <- colorRampPalette(c("darkblue", "#f0f0f0", "#4d4d4d"))
pal <- colorRampPalette(c("#2166ac", "#ece7f2", "#525252")) ##d0d1e6 #"#f0f0f0"
pal <- colorRampPalette(c("#016c59", "#ece7f2", "#525252"))

p_div <- ggplot(stand_div_lf, aes(x = variable, y = common_abbr, fill = value)) + 
    geom_tile(color = "grey", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    #scale_fill_viridis(name = "Standardized \ndiversity", begin = 0, end = 1) +
    #scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
    #   name = "prop. \nhet-exc") +
    #scale_fill_gradientn(colours=c("#081d58", "#253494", "#1d91c0", "#7fcdbb", "#c7e9b4", "#ffffd9"), 
    #    name = "Standardized \ndiversity") +
    #scale_fill_gradientn(colours=pal(5), 
    #    name = "Standardized \ndiversity") +
    scale_fill_distiller(name = "Standardized \ngenetic diversity", palette = "RdYlBu", direction = 1) +
    #scale_fill_distiller(pal(5)) +
    # scale_fill_gradientn(colours=plot_col_div, 
    #   name = "prop. \nhet-exc") +
    theme_tufte(base_family="Hind Guntur Light") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 12),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,22.5,10),"points"), #c(38,-400,7,-260) c(5,-150,7,-260)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_text(hjust=0, size = 10, colour = "#525252"), #abc_cols[stand_div$abc][plot_inds]
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5, 
        text=element_text(family='Hind Guntur Light')) +
    # coord_fixed(ratio = 0.7) +
    scale_x_discrete(labels=c(expression(H[o]),expression(A[r])), # c("Ho", "Ar"),  ##"ARA",
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
bot_res_lf <- bot_res %>% melt(id.vars = c("species", "common", "latin"))

p_bot <- ggplot(bot_res_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "grey", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    #scale_fill_viridis(option = "magma", direction = -1,  
    #    name = "% of loci with\nheterozyosity excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,0.95)) +
    #scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
    #    name = "prop. \nhet-exc") +
    #scale_fill_gradientn(colours=plot_col, 
    #    name = "Prop. of loci \nin heterozyosity excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,1)) +
    scale_fill_distiller(palette="RdBu",
        name = "Prop. of loci with\nheterozyosity-excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,0.95),
        direction = -1) +
    # scale_fill_gradientn(colours=rev(pal2(5)),
    #    name = "% of loci with\nheterozyosity excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,0.95)) +
    theme_tufte(base_family="Hind Guntur Light") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,6,10),"points"),# c(38,-770,7, -690) c(38,300,3, 0)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(), 
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5, 
        text=element_text(family='Hind Guntur Light')) +
    # scale_x_discrete(labels=c("PHE","M")) +
    scale_x_discrete(labels=c("TPM 70", "TPM 80","TPM 90", "SMM ")) + # "IAM",
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 6.2, barheight = 0.5, 
        title.position = "top")) 

plot_grid(p, p_div, p_bot, nrow = 1)


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
    theme_tufte(base_family="Hind Guntur Light") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 12),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,14,10),"points"), #c(38,-520,3, -660)c(38,-520,1, -660)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        legend.title.align = 0.5,
        text=element_text(family='Hind Guntur Light')) +
    scale_x_discrete(labels=c(expression(P[bot]),expression(P[neut]))) +
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 3.6, barheight = 0.5, 
        title.position = "top")) 

p_final <- plot_grid(p, p_div, p_bot, p_abc, nrow = 1, rel_widths = c(0.32, 0.27, 0.15, 0.09))
p_final


ggplot2::ggsave(paste0("other_stuff/figures/figures_final/phylo_plot_color", save_files, ".pdf"), p_final,
    height = 6, width = 9.2, device = "pdf")

Sys.setenv(R_GSCMD = "/usr/local/bin/gs")
extrafont::embed_fonts(paste0("other_stuff/figures/figures_final/phylo_plot_color", save_files, ".pdf"))




# favorites Nunito Average2San Hind2Guntur Hind2Guntur Light
# 












