# Downloaded phylogeny from 10k trees project http://10ktrees.nunn-lab.org/Carnivora/dataset.html
# added ringed seal subspecies manually

library(ape)
library(ggtree)
tree <- read.nexus(file = "data/raw/phylogeny/Pinnipeds25_tree10k.nex")
plot(tree)

# create ringed seal phylogeny
ringed_seal_tree <- read.tree(text = "(((Phoca_hispida_saimensis:0.5,Phoca_hispida_ladogensis:0.5):0.5,Phoca_hispida_botnica:1):0.5, Phoca_hispida_hispida:1.5);")
plot(ringed_seal_tree)
tree_full <- bind.tree(tree, ringed_seal_tree, where = 24, position = 0)
plot(tree_full)
nodelabels() # shorten Node 44
tiplabels()
# remove some part of the edge to align species along the tree
tree_full$edge.length[30] <- tree_full$edge.length[30] - 1.5
plot(tree_full)
nodelabels() # shorten Node 44
tiplabels()
length(tree$tip.label)
?drop.tip

write.tree(tree_full, file = "data/raw/phylogeny/28_species_10ktrees.tre")

ggtree(tree, tiplabels())
ggtree(tree_full) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

tree_test <- add.species.to.genus(tree, "Pusa hispida saimanensis", where = "Pusa hispida")


