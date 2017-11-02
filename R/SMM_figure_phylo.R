# plot phylogeny
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
library(extrafont)
# ggthemr('dust')
source("R/martin.R")

# load all datasets
seals <- read_csv("data/processed/seal_data_complete_rarefac10_29.csv")
# load model probablities from ABC
model_probs <- read_delim("data/processed/sims_5000k_large_bot_model_selection.txt", 
    delim = " ", col_names = c("species", "bot", "neut"), skip = 1)
model_probs <- model_probs %>% rbind(., data.frame("species" = "subantarctic_fur_seal", "bot" = 0.4, "neut" = 0.6))
# are all names overlapping?
sum(seals$species %in% seals$species)

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
# check which ones match # all
which(!(stats_lab %in% tree_lab))

# how to reorder the data to attach the tip.labels?
reorder_data_for_tiplab <- as.numeric(unlist(sapply(tree_lab, function(x) which(stats_lab == x))))
# bind tip labels
seals <- seals %>% 
    .[reorder_data_for_tiplab, ] %>% 
    bind_cols(tibble(tip_label = tree_final$tip.label), .) 

ggtree(tree_final) + geom_tiplab() + xlim(0, 100)

# creat a data.frame where factor levels are in the ggtree plotting sequence
plotting_sequence <- c("weddell_seal", "leopard_seal" ,  "crabeater_seal","ross_seal", "ses", "nes",
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
        common = fct_inorder(factor(common)))

# write to file
if(!file.exists("data/processed/all_stats_tree.csv_29")){
    write_excel_csv(all_stats_tree, "data/processed/all_stats_tree_29.csv")
}
all_stats_tree <- read_csv("data/processed/all_stats_tree_29.csv")

all_stats_tree[all_stats_tree$BreedingType == "both", "BreedingType"] <- "land" 
all_stats_tree <- all_stats_tree %>% mutate(BreedingType = as.factor(as.character(BreedingType))) %>% data.frame()



# create data.frames  -------------
stats_df <- all_stats_tree %>% 
                rename(taxa = tip_label) %>% 
                dplyr::select(c("taxa", "species","common", "bot")) %>% 
                mutate(model = ifelse(bot > 0.5, "bot", "neut"))



p <- ggtree(tree_final, layout = "fan", open.angle = 110)
p <- p %<+% stats_df
p + geom_tiplab2(aes(label = common))
p + geom_tiplab2(aes(label = common), offset = 3) +
    xlim(0,70) 


# resort data frame for plotting in ggtree
df_tree <- data.frame("taxa" = tree_final$tip.label)
tree_sorted_df <- left_join(df_tree, stats_df)

# create matrix with all nodes
new_df <- data.frame(matrix(ncol = ncol(tree_sorted_df), nrow = 27))
names(new_df) <- names(tree_sorted_df)
tree_df <- rbind(tree_sorted_df, new_df)
tree_df <- data.frame(node = 1:56, tree_df)

col_df <- mutate(tree_df, 
    color = case_when(is.na(model) ~ "grey",
        model == "bot"~ "grey", #"cornflowerblue"
        model == "neut" ~ "grey")) 

p <- ggtree(tree_final, layout = "fan", open.angle = 180, aes(color = I(color)))
p <- p %<+% col_df
p + geom_tiplab2(aes(label = common, color = I(color)), offset =3, size = 3) +
    xlim(0,60)






tree_final$

names(col_df)[2] <- "tip.label"

left_join(data.frame("tip.label"= tree_final$tip.label), col_df, by = "tip.label")

##### newer
#tree_final$tip.label <- as.character(tree_final$tip.label)

p <- ggtree(tree_final, layout = "fan", open.angle = 110)
p + geom_tiplab2()


col_df
tree_final$edge

p + geom_tiplab2()
p <- p %<+% col_df
p + geom_tiplab2(aes(label = common))

p$facet #geom_point(aes(color = BreedingType)) +
p + geom_tiplab2(aes(label = common), offset = 3) +
    xlim(0,70) +
    aes(color = I(color))

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
    theme(plot.margin=unit(c(52, -5,13,10),"points"), #c(30,-100,20,0) unit(c(50,-50,20,0) #c(52, -5,10,10)
        legend.position= c(0.26,0.90), #legend.direction = "horizontal",
        legend.spacing = unit(5, "points"),
        legend.key.height=unit(1,"line"),
        legend.key.width = unit(1, "line"), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5,
        text=element_text(family='Lato')) 
# geom_image(data = d_img, aes(image = image), size = 0.2)
p


# p1 <- flip(p,  42, 49)

# diversity
stand_div_lf <- stand_div %>% melt(id.vars = c("species","common"), 
    measure.vars = c("obs_het_mean", "num_alleles_mean", #"mean_allele_range",
        "prop_low_afs_mean"))

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
#pal <- colorRampPalette(c("#d8b365", "#f5f5f5", "#5ab4ac"))
#pal <- colorRampPalette(c("#ef8a62", "#f7f7f7", "#67a9cf"))
#pal <- colorRampPalette(c("#8c510a", "#e0e0e0", "#01665e"))
#pal <- colorRampPalette(c("#d6604d", "#e0e0e0", "#878787"))
#pal <- colorRampPalette(c("#542788", "#e0e0e0", "#b35806"))
# cols <- ifelse(stand_div_lf$abc == "neut", "grey", "red")
pal <- colorRampPalette(c("#2166ac", "#ece7f2", "#525252")) ##d0d1e6 #"#f0f0f0"
pal <- colorRampPalette(c("#016c59", "#ece7f2", "#525252"))

p_div <- ggplot(stand_div_lf, aes(x = variable, y = common, fill = value)) + 
    geom_tile(color = "grey", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    #scale_fill_viridis(name = "Standardized \ndiversity", begin = 0, end = 1) +
    #scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
    #   name = "prop. \nhet-exc") +
    #scale_fill_gradientn(colours=c("#081d58", "#253494", "#1d91c0", "#7fcdbb", "#c7e9b4", "#ffffd9"), 
    #    name = "Standardized \ndiversity") +
    #scale_fill_gradientn(colours=pal(5), 
    #    name = "Standardized \ndiversity") +
    scale_fill_distiller(name = "Standardized \ndiversity", palette = "RdYlBu", direction = 1) +
    #scale_fill_distiller(pal(5)) +
    # scale_fill_gradientn(colours=plot_col_div, 
    #   name = "prop. \nhet-exc") +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,5,10),"points"), #c(38,-400,7,-260) c(5,-150,7,-260)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_text(hjust=0, size = 9, colour = "#525252"), #abc_cols[stand_div$abc][plot_inds]
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5, 
        text=element_text(family='Lato')) +
    # coord_fixed(ratio = 0.7) +
    scale_x_discrete(labels=c("HET", "AR", "LFA"),  ##"ARA",
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
        name = "% of loci with\nheterozyosity excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,0.95),
        direction = -1) +
    # scale_fill_gradientn(colours=rev(pal2(5)),
    #    name = "% of loci with\nheterozyosity excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,0.95)) +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 9),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,1,10),"points"),# c(38,-770,7, -690) c(38,300,3, 0)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5, 
        text=element_text(family='Lato')) +
    # scale_x_discrete(labels=c("PHE","M")) +
    scale_x_discrete(labels=c("70", "80","90", "100 ")) + # "IAM",
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 5.6, barheight = 0.5, 
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
        name = "ABC \nprob. %", labels=c(0, 50, 100), breaks = c(0,0.5,1),
        direction = -1) +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 9),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,2,10),"points"), #c(38,-520,3, -660)c(38,-520,1, -660)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5,
        text=element_text(family='Lato')) +
    scale_x_discrete(labels=c("Bot","Con")) +
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 2.5, barheight = 0.5, 
        title.position = "top")) 

p_final <- plot_grid(p, p_div, p_bot, p_abc, nrow = 1, rel_widths = c(0.4, 0.27, 0.17, 0.10))
p_final

#p_final <- plot_grid(p, p_div, p_bot, p_abc, nrow = 1, rel_widths = c(0.2, 0.16, 0.09, 0.05, 0.02))
#p_final

save_plot("figures/phylo_plot2.jpg", p_final,
    ncol = 2, # we're saving a grid plot of 2 columns
    nrow = 1, # and 2 rows
    # each individual subplot should have an aspect ratio of 1.3
    # base_aspect_ratio = 0.9,
    base_height = 6,
    base_width = 5
)










