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
# ggthemr('dust')
# prepare data
# load descriptive data
load("data/processed/seal_ss_rarefaction40ind.RData")
# load all datasets
seals <- read_excel("data/processed/seal_data_complete.xlsx")
names(seals)[11] <- "IUCN_rating"

# load model probablities from ABC
model_probs <- read.table("data/processed/sims_5000k_large_bot_model_selection.txt")

# are all names overlapping?
sum(rownames(all_sumstats_full) %in% seals$species)

ind_reorder <- NULL
for (i in seals$species){
    ind_reorder <- c(ind_reorder, which(rownames(all_sumstats_full) == i))
}

# reorder sumstats according to seal descriptives file
sumstats <- all_sumstats_full[ind_reorder, ]
modelprobs <- model_probs[ind_reorder, ]
# create new data.frame as a combination
all_stats <- cbind(seals[c(1:17, 36:56)], sumstats, modelprobs)

# load higdon phylogeny
tree_final <- read.tree("data/raw/phylogeny/higdon_mod_28.tre")
plot(tree_final)
tree_final$tip.label

# sort data like phylogeny -------------------------------------------------------------------------
tree_lab <- unlist(lapply(strsplit(tree_final$tip.label, "_"), function(x) x[length(x)]))
tree_lab
stats_lab <- unlist(lapply(strsplit(all_stats$latin, " "), function(x) x[length(x)]))
# check which ones match
stats_lab %in% tree_lab
# manipulate
stats_lab[7] <- "pusillus"
stats_lab[11] <- "bryonia"
stats_lab[18] <- "weddellii"
stats_lab[27] <- "botnica"
# everything matched now
stats_lab %in% tree_lab

# now reorder stats (not really useful as tree plotting is different from tree label sequence)
reorder_ind <- unlist(sapply(tree_lab, function(x) which(stats_lab == x)))
all_stats_tree <- all_stats[as.numeric(reorder_ind), ]

# row names in all_stats_tree apparently have to match the tree tip labels
# mutate has to be before naming rows
all_stats_tree <- mutate(all_stats_tree, abc = ifelse(bot > 0.5, "bot", "neut"))
rownames(all_stats_tree) <- tree_final$tip.label
id <- rownames(all_stats_tree)
# create data.frames for heatmaps -------------
# diversity
stand_div <- all_stats_tree[c("obs_het_mean", "num_alleles_mean",  
                            "mean_allele_range", "prop_low_afs_mean")] %>% 
    apply(2, scale) %>% 
    as.data.frame()
# for ggtree
# stand_div <- cbind(id, stand_div)

stand_div <- cbind(all_stats_tree$latin, all_stats_tree$seal_names_real, stand_div, all_stats_tree$abc)
names(stand_div)[1] <- "species"
names(stand_div)[2] <- "species_common"
names(stand_div)[ncol(stand_div)] <- "abc"

# reverse df for plotting
correct_levels_latin <- as.character(stand_div$species)
correct_levels_common <- as.character(stand_div$species_common)
# to check correct sequence just plot tree with node number
plotTree(tree_final)
tiplabels()
edgelabels()
# plot inds needed within diversity plot later
plot_inds <- c(1:2, 10:12, 3:5, 9, 6:8, 13,14, 19,20, 18,17,15,16,27,28,25,26,21,22,23,24)
plot_seq_latin <- correct_levels_latin[c(1:2, 10:12, 3:5, 9, 6:8, 13,14, 19,20, 18,17,15,16,27,28,25,26,21,22,23,24)]
stand_div$species <- factor(stand_div$species, levels = plot_seq_latin)

#plot_seq_common <- correct_levels_common[c(1:2, 10:12, 3:5, 9, 6:8, 13,14, 19,20, 18,17,15,16,27,28,25,26,21,22,23,24)] 
# stand_div$species_common <- factor(stand_div$species_common, levels = plot_seq_common)
# rownames(stand_div) <- tree_final$tip.label

p <- ggtree(tree_final) + geom_text(aes(label=node))
p
# bottleneck
bot_res <- all_stats_tree[c("IAM_ratio", "TPM70_ratio", "TPM90_ratio","TPM95_ratio", "SMM_ratio")]   #%>% #"mratio_mean"
            # mutate(mratio_mean= 1-mratio_mean) %>% 
            #apply(2, scale) %>% 
            #as.data.frame()
bot_res <- cbind(all_stats_tree$latin, bot_res)
names(bot_res)[1] <- "species"
# correct_levels <- as.character(bot_res$species)
# to check correct sequence just plot tree with node number
bot_res$species <- factor(bot_res$species, levels = plot_seq_latin)

# ABC probs
abc_probs <- all_stats_tree[c("bot", "neut")]
abc_probs <- cbind(all_stats_tree$latin, abc_probs)
names(abc_probs)[1] <- "species"
# correct_levels <- as.character(abc_probs$species)
# to check correct sequence just plot tree with node number
abc_probs$species <- factor(abc_probs$species, levels = plot_seq_latin)

# abc star
abc_star <- abc_probs_lf %>% filter(variable == "bot") %>% mutate(bot = ifelse(value > 0.5, "bot", "const"))
abc_star$species <- factor(abc_star$species, levels = plot_seq_latin)

# add nodes
# plotTree(tree_final,node.numbers=T)
# id <- rownames(all_stats_tree)
# all_stats_tree <- cbind(id, all_stats_tree)

# group into families
cls <- list(phocids=tree_final$tip.label[13:28],
            otarids= tree_final$tip.label[2:12],
            odobenids = tree_final$tip.label[1])
tree_final <- groupOTU(tree_final, cls)

#tree_final <- groupClade(tree_final, node = 31)
#tree_final <- groupClade(tree_final, node = 1)
#img <- c("phylo_figures/AFS.jpg")
#names(img) <- "31"

# create tree view
p <- ggtree(tree_final, color = "#969696") #aes(color = group)

all_stats_for_tree <- cbind(rownames(all_stats_tree), all_stats_tree)
all_stats_for_tree$IUCN_rating[is.na(all_stats_for_tree$IUCN_rating)] <- "data deficient"
all_stats_for_tree$IUCN_rating <- factor(all_stats_for_tree$IUCN_rating, levels = c("least concern", "vulnerable", "near threatened", "endangered",
    "data deficient"))
rownames(all_stats_for_tree) <- NULL
names(all_stats_for_tree)[1] <- "taxa"

#test <- all_stats_for_tree %>% 
#    mutate(abund_classes = ifelse(Abundance > 1000000, "1000k", )

# add data
p <- p %<+% all_stats_for_tree 

p + geom_tippoint() #aes(color=IUCN_rating)

#library(ggthemr)
#ggthemr(text_size = 3)

# pictures
img <- c("phylo_figures/NES.jpg", "phylo_figures/AFS.jpg", "phylo_figures/med_monk_seal.jpg")
d_img <- data.frame(x = c(0.5, 1, 2),
    y = c(0.5, 1, 2),
    image = img)

p <- p + #layout="circular" , "fan"open.angle=180
    # scale_color_manual(values=c("black", "#9e9ac8","#6a51a3", "#bcbddc" )) +
    geom_tippoint(aes(size = Abundance, fill=IUCN_rating),stroke=0.3, color = "#737373", shape=21) + #shape=21,stroke=0.5
    # ggtitle("Pinniped Phylogeny") +
    #geom_tiplab(size=5, color="black") + # tiplap 2 for circular
    #ggplot2::xlim(0, 15) +
    # scale_color_manual(values = c("#66c2a5", "#fdae61", "#f46d43", "#d53e4f", "#000000")) +
    # scale_color_manual(values = c("lightblue", "goldenrod", "red", "darkred", "grey")) +
    scale_fill_manual(values = c("#fee8c8", "#fc8d59", "#d7301f", "#7f0000", "lightgrey")) +
    #scale_color_manual(values = c("#c7e9b4", "#41b6c4", "#225ea8", "#081d58", "lightgrey")) +
    scale_size_continuous(range = c(1.5, 5), trans = "sqrt", breaks=c(1000, 10000, 100000, 1000000),
                      labels = c(expression(10^{3}), expression(10^{4}), expression(10^{5}), expression(10^{6}))
                       ) + #range = c(0.1,5), , breaks=c(500, 1000, 10000, 100000, 1000000)
    guides(fill = guide_legend(title = "IUCN rating", title.position = "top", direction = "vertical", order = 2),
           size = guide_legend(title.position = "top", title = "Global abundance", direction = "horizontal", order = 1)
           ) +
    theme(plot.margin=unit(c(52, -5,10,10),"points"), #c(30,-100,20,0) unit(c(50,-50,20,0)
          legend.position= c(0.26,0.95), #legend.direction = "horizontal",
          legend.spacing = unit(5, "points"),
        legend.key.height=unit(1,"line"),
        legend.key.width = unit(1, "line"), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5) 
    # geom_image(data = d_img, aes(image = image), size = 0.2)


            
   #c(30,-50,20,0)
    #geom_text(aes(label=node)) +
    # geom_cladelabel(node=41, label="Phocidae", angle = 270, align=T, hjust='center', offset = 5, 
    #     offset.text=1, color = "cornflowerblue",  barsize=0.5) +
    # geom_cladelabel(node=31, label="Otariidae", align=T,angle = 270, hjust='center', offset = 5,
    #     offset.text=1, barsize=0.5, color = "goldenrod") +
    # geom_cladelabel(node=1, label="Odobenidae", align=T,angle = 270, color = "darkgrey", hjust=0.7, 
    #     offset.text=6, barsize=0.5) +
    #scale_color_manual(values=c("black","darkgrey", "goldenrod", "cornflowerblue" )) +




# p1 <- flip(p,  42, 49)

# diversity
stand_div_lf <- stand_div %>% melt(id.vars = c("species","species_common", "abc"), 
                measure.vars = c("obs_het_mean", "num_alleles_mean", "mean_allele_range",
                                 "prop_low_afs_mean"))

#plot_col_div <- viridis(20)
#plot_col_div <- rev(plot_col_div[c(rep(FALSE,4), TRUE)]) # get every 5th element

plot_col_div <- viridis(100)
plot_col_div <- rev(plot_col_div[c(rep(FALSE,30), TRUE)])

#color names
#numColors <- length(levels(stand_div_lf$abc))
#getColors <- brewer_pal('qual')
#myPalette <- getColors(numColors)
abc_cols <- c("#a50f15", "#525252")

# cols <- ifelse(stand_div_lf$abc == "neut", "grey", "red")
p_div <- ggplot(stand_div_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "white", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    # scale_fill_viridis() +
    #scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
     #   name = "prop. \nhet-exc") +
    scale_fill_gradientn(colours=c("#081d58", "#253494", "#1d91c0", "#7fcdbb", "#c7e9b4", "#ffffd9"), 
        name = "Standardized \ndiversity") +
   # scale_fill_gradientn(colours=plot_col_div, 
    #   name = "prop. \nhet-exc") +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,0,10),"points"), #c(38,-400,7,-260) c(5,-150,7,-260)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_text(hjust=0, size = 9, colour = "#525252"), #abc_cols[stand_div$abc][plot_inds]
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5) +
    # coord_fixed(ratio = 0.7) +
    scale_x_discrete(labels=c("HET", "AR","ARA","LFA"), 
        position = "bottom") +
    guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5, 
            title.position = "top", label.position = "bottom")) 

# grid.arrange(p, p_div, nrow = 1)
plot_grid(p, p_div, nrow = 1)

# het-excess and mratio
plot_col <- magma(40)
plot_col <- rev(plot_col[c(rep(FALSE,2), TRUE)]) # get every 5th element

bot_res_lf <- bot_res %>% melt(id.vars = "species")
p_bot <- ggplot(bot_res_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "white", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    #scale_fill_viridis(option = "viridis", direction = -1) +
    #scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
    #    name = "prop. \nhet-exc") +
    scale_fill_gradientn(colours=plot_col, 
        name = "Prop. of loci \nin heterozyosity excess", labels=c(0, 0.5, 1.0), breaks = c(0,0.5,1)) +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 9),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,0,10),"points"),# c(38,-770,7, -690) c(38,300,3, 0)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(), 
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5) +
    # scale_x_discrete(labels=c("PHE","M")) +
    scale_x_discrete(labels=c("IAM", "70", "90","95", "SMM")) +
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 0.5, 
        title.position = "top")) 

plot_grid(p, p_div, p_bot, nrow = 1)

# abs results
plot_col_abc <- magma(40)
plot_col_abc <- rev(plot_col_abc[c(rep(FALSE,2), TRUE)])

abc_probs_lf <- abc_probs %>% melt(id.vars = "species")
p_abc <- ggplot(abc_probs_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "white", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    # scale_fill_viridis(option = "viridis", direction = -1) +
    scale_fill_gradientn(colours=plot_col_abc, 
        name = "ABC \nprob. %", labels=c(0, 50, 100), breaks = c(0,0.5,1)) +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1, size = 9),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,0,10),"points"), #c(38,-520,3, -660)c(38,-520,1, -660)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        legend.title.align = 0.5) +
    scale_x_discrete(labels=c("Bott","Const")) +
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 2.5, barheight = 0.5, 
        title.position = "top")) 

p_final <- plot_grid(p, p_div, p_bot, p_abc, nrow = 1, rel_widths = c(0.2, 0.18, 0.09, 0.05))
p_final

# some stars

# create a column of stars
p_abc_star <- ggplot(abc_star, aes(x = variable, y = species, color = bot)) + 
    geom_point(size = 2.5, alpha = 1) +
    scale_color_manual(values = c("cornflowerblue", "white")) +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_blank(),
        # axis.text.x = element_blank(),
        legend.position="none",
        plot.margin=unit(c(59,10,23,-5),"points"), #c(38,-520,3, -660)c(38,-520,1, -660)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_blank()
        )

p_final <- plot_grid(p, p_div, p_bot, p_abc, p_abc_star, nrow = 1, rel_widths = c(0.2, 0.18, 0.09, 0.05, 0.02))
p_final

save_plot("phylo_plot.pdf", p_final,
    ncol = 2, # we're saving a grid plot of 2 columns
    nrow = 1, # and 2 rows
    # each individual subplot should have an aspect ratio of 1.3
   # base_aspect_ratio = 0.9,
    base_height = 6,
    base_width = 5
)



