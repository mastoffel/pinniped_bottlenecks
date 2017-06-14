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
# reorder_ind <- unlist(sapply(tree_lab, function(x) which(stats_lab == x)))
# all_stats_tree <- all_stats[as.numeric(reorder_ind), ]

# row names in all_stats_tree apparently have to match the tree tip labels
#rownames(all_stats_tree) <- tree_final$tip.label
#id <- rownames(all_stats_tree)

# create data.frames for heatmaps -------------
# diversity
stand_div <- all_stats_tree[c("num_alleles_mean", "obs_het_mean" , 
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
bot_res <- all_stats_tree[c("TPM70_ratio", "mratio_mean")] %>% 
            mutate(TPM70_ratio = 1-TPM70_ratio) %>% 
            apply(2, scale) %>% 
            as.data.frame()
bot_res <- cbind(id, bot_res)

# ABC probs
abc_probs <- all_stats_tree[c("bot", "neut")]
rownames(abc_probs) <- tree_final$tip.label


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

p <- ggtree(tree_final, aes(color=group)) + #layout="circular" , "fan"open.angle=180
    ggtitle("Pinniped Phylogeny") +
    geom_tiplab(size=2, color="black") + # tiplap 2 for circular
    ggplot2::xlim(0, 15) +
    geom_text(aes(label=node)) +
    # geom_cladelabel(node=41, label="Phocidae", angle = 270, align=T, hjust='center', offset = 5, 
    #     offset.text=1, color = "cornflowerblue",  barsize=0.5) +
    # geom_cladelabel(node=31, label="Otariidae", align=T,angle = 270, hjust='center', offset = 5,
    #     offset.text=1, barsize=0.5, color = "goldenrod") +
    # geom_cladelabel(node=1, label="Odobenidae", align=T,angle = 270, color = "darkgrey", hjust=0.7, 
    #     offset.text=6, barsize=0.5) +
    scale_color_manual(values=c("black","darkgrey", "goldenrod", "cornflowerblue" ))
p
p1 <- flip(p, 42, 49)

# diversity
stand_div_lf <- stand_div %>% melt(id.vars = "species")

p_div <- ggplot(stand_div_lf, aes(x = variable, y = species, fill = value)) + 
    geom_tile(color = "white", size = 0.1) +
    # labs(x = "microsat mutation model", y = "") +
    scale_fill_viridis(option="magma") +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        # axis.text.x = element_text(angle = 50, hjust = 1),
        axis.text.x = element_blank(),
        legend.position="non",
        plot.margin=unit(c(30,0,20,0),"points"),
        axis.title.x=element_blank()) +
    coord_fixed(ratio = 0.7) 
grid.arrange(p, p_div, nrow = 1)


bot_res_lf <- bot_res %>% melt(id.vars = "id") %>% mutate(.panel='bottleneck')
p3 <- facet_plot(p2 + xlim_tree(20), panel = "bottleneck", data = bot_res_lf, geom_tile, 
    aes(x= as.numeric(as.factor(variable)), fill = value)) +  
    theme_tree2() 
p3





plot_grid(p, p_bot_ratio, nrow = 1, rel_widths = c(0.3, 0.7))


gheatmap(p, stand_div, offset = 15, width=0.5, font.size=3, colnames_angle=-45, hjust=0) +
    scale_fill_viridis(option = "magma") 



 facet_plot(p, 'heatmap', data = test_stats, geom_tile, aes(x=BreedingType))


facet_plot(p, 'heatmap', data = all_stats_tree, geom_tile, aes(x=BreedingType))


    

p2 <- facet_plot(p, panel="dot", data=test, geom=geom_point, mapping = aes(x=mratio_mean))


gheatmap(p,  all_stats_tree[23:27], offset = 5, width=0.5, font.size=3) +
        + scale_fill_manual(values=c("steelblue", "firebrick", "darkgreen"))

df <- as.data.frame(seals[c("allelic_richness", "mean_het")])
pp <- p %>%
    gheatmap(df, offset=8, width=0.6, colnames=FALSE) 
pp


# p <- ggtree(tree_final, layout="circular" ) + #layout="circular" 
#     ggtitle("Pinniped Phylogeny") +
#     geom_tiplab2(size=2, color="black") +
#     ggplot2::xlim(0, 30) +
#     geom_hilight(node=41,  alpha=.6, fill = "cornflowerblue") +
#     geom_hilight(node=31,  alpha=.6, fill = "goldenrod") +
#     geom_hilight(node=1,  alpha=.6, fill = "darkgreen")


