## modify the higdon phylogeny

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

# some command
plotTree(tree)
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

write.tree(tree_final, file = "data/raw/phylogeny/higdon_mod_28.tre", append = FALSE,
    digits = 10, tree.names = FALSE)