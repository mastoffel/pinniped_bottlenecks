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
tree <- read.tree("data/raw/phylogeny/higdon.tre")

# species not in the study
to_drop <- c("P._caspica", "P._sibirica", "P._largha", "P._groenlandica", "H._fasciata", 
             "M._tropicalis", "A._townsendi", "A._phillippii", "N._cinerea",
             "A._tropicalis")
tree2 <- drop.tip(tree, to_drop)

# add nodes
plotTree(tree2,node.numbers=T)
# some command
# plot(tree)
# nodelabels()
# tiplabels()

# add galapagos sea lion and ringed seals
tree_final <- tree2 %>% 
    bind.tip("Z._wollebaeki", edge.length = 0, where=which(tree2$tip.label=="Z._californianus")) %>%
    bind.tip("P._hispida_hispida", where=which(.$tip.label=="P._hispida"), position=1) %>%
    bind.tip("P._hispida_saimensis" ,where=which(.$tip.label=="P._hispida"), position=1) %>%
    bind.tip("P._hispida_ladogensis", where=which(.$tip.label=="P._hispida"), position=1) %>%
    bind.tip("P._hispida_botnica", where=which(.$tip.label=="P._hispida"), position=1) %>%
    drop.tip("P._hispida")
plotTree(tree_final)


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

stand_div <- cbind(all_stats_tree$latin, stand_div)
names(stand_div)[1] <- "species"
# reverse df for plotting
correct_levels <- as.character(stand_div$species)
# to check correct sequence just plot tree with node number
plot_seq <- correct_levels[c(1:2, 10:12, 3:5, 9, 6:8, 13,14, 19,20, 18,17,15,16,27,28,25,26,21,22,23,24)]
stand_div$species <- factor(stand_div$species, levels = plot_seq)
# rownames(stand_div) <- tree_final$tip.label

# bottleneck
bot_res <- all_stats_tree[c("IAM_ratio", "TPM70_ratio", "TPM90_ratio","TPM95_ratio", "SMM_ratio")] #%>% #"mratio_mean"
            # mutate(mratio_mean= 1-mratio_mean) %>% 
            apply(2, scale) %>% 
            as.data.frame()
bot_res <- cbind(all_stats_tree$latin, bot_res)
names(bot_res)[1] <- "species"
correct_levels <- as.character(bot_res$species)
# to check correct sequence just plot tree with node number
bot_res$species <- factor(bot_res$species, levels = plot_seq)

# ABC probs
abc_probs <- all_stats_tree[c("bot", "neut")]
abc_probs <- cbind(all_stats_tree$latin, abc_probs)
names(abc_probs)[1] <- "species"
correct_levels <- as.character(abc_probs$species)
# to check correct sequence just plot tree with node number
abc_probs$species <- factor(abc_probs$species, levels = plot_seq)


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
p <- ggtree(tree_final)

all_stats_for_tree <- cbind(rownames(all_stats_tree), all_stats_tree)
all_stats_for_tree$IUCN_rating[is.na(all_stats_for_tree$IUCN_rating)] <- "data deficient"
all_stats_for_tree$IUCN_rating <- factor(all_stats_for_tree$IUCN_rating, levels = c("least concern", "vulnerable", "near threatened", "endangered",
    "data deficient"))
rownames(all_stats_for_tree) <- NULL
names(all_stats_for_tree)[1] <- "taxa"

p <- p %<+% all_stats_for_tree 

p + geom_tippoint(aes(color=IUCN_rating))

library(ggthemr)
ggthemr(text_size = 3)


p <- p + #layout="circular" , "fan"open.angle=180
    geom_tippoint(aes(size = Abundance, color=IUCN_rating)) +
    # ggtitle("Pinniped Phylogeny") +
    #geom_tiplab(size=5, color="black") + # tiplap 2 for circular
    #ggplot2::xlim(0, 15) +
    scale_color_manual(values = c("#66c2a5", "#fdae61", "#f46d43", "#d53e4f", "#000000")) +
   # scale_color_manual(values = c("#440154FF", "#3E4A89FF", "#25848EFF", "#38B977FF", "#000000")) +
    scale_size_continuous(range = c(0.1,3), trans = "log", breaks=c(1000, 10000, 100000, 1000000),
                      labels = c("1k", "10k", "100k", "1000k")
                       ) + #range = c(0.1,5), , breaks=c(500, 1000, 10000, 100000, 1000000)
    guides(colour = guide_legend(title = "IUCN Rating", title.position = "top", direction = "vertical", order = 2),
           size = guide_legend(title.position = "top", title = "Global Abundance", direction = "horizontal", order = 1)
           ) +
    theme(plot.margin=unit(c(50, 10,10,10),"points"), #c(30,-100,20,0) unit(c(50,-50,20,0)
          legend.position= c(0.26,0.95), #legend.direction = "horizontal",
          legend.spacing = unit(0, "points"),
        legend.key.height=unit(0.7,"line"),
        legend.key.width = unit(0.7, "line"), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title.align = 0.5) 


            
   #c(30,-50,20,0)
    #geom_text(aes(label=node)) +
    # geom_cladelabel(node=41, label="Phocidae", angle = 270, align=T, hjust='center', offset = 5, 
    #     offset.text=1, color = "cornflowerblue",  barsize=0.5) +
    # geom_cladelabel(node=31, label="Otariidae", align=T,angle = 270, hjust='center', offset = 5,
    #     offset.text=1, barsize=0.5, color = "goldenrod") +
    # geom_cladelabel(node=1, label="Odobenidae", align=T,angle = 270, color = "darkgrey", hjust=0.7, 
    #     offset.text=6, barsize=0.5) +
    #scale_color_manual(values=c("black","darkgrey", "goldenrod", "cornflowerblue" )) +
# p <- inset(p, img)
p


# p1 <- flip(p,  42, 49)

# diversity
stand_div_lf <- stand_div %>% melt(id.vars = "species")

#plot_col_div <- viridis(20)
#plot_col_div <- rev(plot_col_div[c(rep(FALSE,4), TRUE)]) # get every 5th element

plot_col_div <- magma(100)
plot_col_div <- rev(plot_col_div[c(rep(FALSE,2), TRUE)])

p_div <- ggplot(stand_div_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "white", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    # scale_fill_viridis() +
    #scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
     #   name = "prop. \nhet-exc") +
    scale_fill_gradientn(colours=c("#081d58", "#253494", "#1d91c0", "#7fcdbb", "#c7e9b4", "#ffffd9"), 
        name = "standardized \ndiversity") +
    #scale_fill_gradientn(colours=plot_col_div, 
    #   name = "prop. \nhet-exc") +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,0,10),"points"), #c(38,-400,7,-260) c(5,-150,7,-260)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_text(hjust=0),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title.align = 0.5) +
    # coord_fixed(ratio = 0.7) +
    scale_x_discrete(labels=c("AR","HET","ARA","LFA"), 
        position = "bottom") +
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, 
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
        name = "prop. of loci \nin het.excess") +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,0,10),"points"),# c(38,-770,7, -690) c(38,300,3, 0)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(), 
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title.align = 0.5) +
    # scale_x_discrete(labels=c("PHE","M")) +
    scale_x_discrete(labels=c("IAM", "70", "90","95", "SMM")) +
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 5, barheight = 0.5, 
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
        name = "abc \nprob.") +
    theme_tufte(base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text.x = element_text(angle = 70, hjust = 1),
        # axis.text.x = element_blank(),
        legend.position="top",
        plot.margin=unit(c(0,10,0,10),"points"), #c(38,-520,3, -660)c(38,-520,1, -660)
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        axis.text.y = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title.align = 0.5) +
    scale_x_discrete(labels=c("Bott","Const")) +
    #labs(x=NULL, y=NULL) +
    # coord_fixed(ratio = 0.7) +
    guides(fill = guide_colorbar(barwidth = 2.5, barheight = 0.5, 
        title.position = "top")) 

p_final <- plot_grid(p, p_div, p_bot, p_abc, nrow = 1, rel_widths = c(0.2, 0.2, 0.1, 0.05))
p_final


save_plot("phylo_plot.pdf", p_final,
    ncol = 2, # we're saving a grid plot of 2 columns
    nrow = 1, # and 2 rows
    # each individual subplot should have an aspect ratio of 1.3
    base_aspect_ratio = 1
)



