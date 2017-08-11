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


# load all datasets
seals <- read_csv("data/processed/seal_data_complete_rarefac10.csv")

# load model probablities from ABC
model_probs <- read_delim("data/processed/sims_5000k_large_bot_model_selection.txt", 
    delim = " ", col_names = c("species", "bot", "neut"), skip = 1)

# are all names overlapping?
sum(seals$species %in% seals$species)

# join
seals <- left_join(seals, model_probs, by = "species")


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