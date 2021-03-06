# plots for the talk

library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(stringr)

# load phyologenetic sequence

# load all datasets
phylo <- read_excel("data/raw/phylogeny_overview.xlsx", sheet = 2, col_names = F)
seals <- read_excel("data/processed/all_data_seals.xlsx")
seals


# get phylo order
reorder <- unlist(lapply(phylo$X3, function(x) which(str_detect(seals$species, x))))
seals <- seals[reorder, ]
seals[, 2:6] <- lapply(seals[, 2:6], as.numeric)

# get krueger data and order
krueger_data <- read_excel("data/processed/seal_data_krueger.xlsx", sheet = 1, col_names = T)

# decide at some point whether to include ringed seal subspecies
krueger_data <- krueger_data[1:38, ]
krueger_data <- krueger_data[!is.na(krueger_data$dataset_name), ]

# prepare krueger data for matching
seals$species %in% krueger_data$dataset_name 
# reorder krueger data
reorder <- unlist(lapply(seals$species, function(x) which(str_detect(krueger_data$dataset_name, x))))
krueger_data <- krueger_data[reorder, ]
# bind seals and krueger_data
seals <- cbind(seals, krueger_data)

# rename the hookers sea lion
seal_names_real <- phylo$X2[1:28]
seal_names_real[10] <- "New Zealand Sea Lion"

# delete doubled species column
seals <- seals[-which(duplicated(names(seals)))]
# delete rownames
rownames(seals) <- 1:nrow(seals)


## heatmap for ratio
# plot p_vals
bot_ratio <- melt(seals[, c(1, which(str_detect(names(seals), "_ratio")))], id.vars = "species")
bot_ratio$species <- factor(bot_ratio$species, levels = bot_ratio$species[28:1]) # 


p <- ggplot(bot_ratio, aes(x= variable, y = species, fill = value)) + 
    #facet_grid(.~dataset) + 
    geom_tile(color = "white", size = 0.1) +
    labs(x = "microsat mutation model", y = "") +
    # scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#1d91c0", "#225ea8", "#0c2c84", "#081d58")) +
    #scale_fill_gradientn(colours=c("#fff7fb", "#d0d1e6", "#67a9cf", "#02818a", "#014636")) +
    #scale_fill_gradientn(colours=c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#cb181d", "#a50f15", "#67000d")) +
    # scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#41b6c4", "#225ea8", "#253494", "#081d58")) +
    # for p_val
    # scale_fill_gradientn(colours=c("#081d58",  "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1", "#ffffd9"), 
    #     name = "p-val") +
    
    scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
        name = "prop. \nhet-exc") +
    
    scale_x_discrete(labels=c("IAM", "TPM70", "TPM90", "TPM95", "SMM")) +
    scale_y_discrete(labels = rev(seal_names_real)) +
    # scale_fill_gradient(name = "p-value", label = comma,  breaks = c(0.05, 0.5)) +
    theme_tufte(base_family="Helvetica") +
    coord_equal() +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9), 
        axis.text.x = element_text(angle = 50, hjust = 1))

ggplot2::ggsave(p, 
    filename = "heatmap_ratio.tiff", 
    #  path = "D:/recent R",
    width = 9,
    scale = 0.6,
    height = 11, units = "in",
    dpi = 300)


# plot p_vals
bot_p_vals <- melt(seals[, c(1, which(str_detect(names(seals), "Exc")))], id.vars = "species")
bot_p_vals$species <- factor(bot_p_vals$species, levels = bot_p_vals$species[28:1])

p <- ggplot(bot_p_vals, aes(x= variable, y = species, fill = value)) + 
    #facet_grid(.~dataset) + 
    geom_tile(color = "white", size = 0.1) +
    labs(x = "microsat mutation model", y = "") +
    # scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#1d91c0", "#225ea8", "#0c2c84", "#081d58")) +
    #scale_fill_gradientn(colours=c("#fff7fb", "#d0d1e6", "#67a9cf", "#02818a", "#014636")) +
    #scale_fill_gradientn(colours=c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#cb181d", "#a50f15", "#67000d")) +
    # scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#41b6c4", "#225ea8", "#253494", "#081d58")) +
    # for p_val
    scale_fill_gradientn(colours=c("#081d58",  "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1", "#ffffd9"), 
        name = "p-val") +
    scale_x_discrete(labels=c("IAM", "TPM70", "TPM90", "TPM95", "SMM")) +
    # scale_y_discrete(labels = rev(seal_names_real)) +
    scale_y_discrete(labels = element_blank()) +
    # scale_fill_gradient(name = "p-value", label = comma,  breaks = c(0.05, 0.5)) +
    theme_tufte(base_family="Helvetica") +
    coord_equal() +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9), 
        axis.text.x = element_text(angle = 50, hjust = 1))

ggplot2::ggsave(p, 
    filename = "heatmap_pval.tiff", 
   #  path = "D:/recent R",
    width = 9,
    scale = 0.6,
    height = 11, units = "in",
    dpi = 300)



# harem_size / sexual size dimorphism vs bottleneck
library(ggrepel)
names(seals)
seals$weight_ratio <- seals$male_weight / seals$female_weight

p1 <- ggplot(seals, aes(x= weight_ratio, y = TPM70_ratio, label = seal_names_real)) + 
    #facet_grid(.~dataset) + 
    geom_smooth(method = "lm", size = 0.5, color = "grey", fill = "lightblue") +
    geom_point(size = 2, alpha = 0.5) +
    labs(x = "male / female weight ratio", y = "bottleneck footprint \n (heterozygosity excess ratio)") +
    geom_text_repel(size = 3) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0)))

ggplot2::ggsave(p1, 
    filename = "bottleneck_vs_sizeratio.tiff", 
    #  path = "D:/recent R",
    width = 14,
    scale = 0.6,
    height = 9, units = "in",
    dpi = 300)


p2 <- ggplot(seals, aes(x= weight_ratio, y = TPM70_ratio, label = seal_names_real)) + 
    #facet_grid(.~dataset) + 
    geom_smooth(method = "lm", size = 0.5, color = "grey", fill = "lightblue") +
    geom_point(aes(size = harem_size), alpha = 0.5) +
    labs(x = "male / female size ratio", y = "bottleneck footprint \n (heterozygosity excess ratio)") +
    geom_text_repel(size = 3) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid'),
          axis.title.y=element_text(margin=margin(0,20,0,0)),
          axis.title.x=element_text(margin=margin(20,0,0,0)))
          
ggplot2::ggsave(p2, 
    filename = "bottleneck_vs_sizeratio2.tiff", 
    #  path = "D:/recent R",
    width = 14,
    scale = 0.6,
    height = 9, units = "in",
    dpi = 300)


# some tryouts with ollies data

ggplot(seals, aes(x= mean_SSD, y = TPM70_ratio, label = seal_names_real)) + 
    #facet_grid(.~dataset) + 
    geom_smooth(method = "lm", size = 0.5, color = "grey", fill = "lightblue") +
    geom_point(aes(size = harem_size), alpha = 0.5) +
    labs(x = "male / female size ratio", y = "bottleneck footprint \n (heterozygosity excess ratio)") +
    geom_text_repel(size = 3) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0)))


        
# inbreeding

ggplot(seals, aes(x= Var.112, y = TPM70_ratio, label = seal_names_real)) + 
    #facet_grid(.~dataset) + 
    geom_smooth(method = "lm", size = 0.5, color = "grey", fill = "lightblue") +
    geom_point(aes(size = harem_size), alpha = 0.5) +
    labs(x = "male / female size ratio", y = "bottleneck footprint \n (heterozygosity excess ratio)") +
    geom_text_repel(size = 3) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0)))


library(dplyr)
seals_g2 <- seals
seals_g2$real_names_species <- phylo$X2[1:28]
seals_g2$real_names_species[10] <- "New Zealand Sea Lion"
seals_g2$species <- factor(seals_g2$species, levels = seals_g2$species[order(seals_g2$g2)])
seals_g2$real_names_species <- factor(seals_g2$real_names_species, levels = seals_g2$real_names_species[order(seals_g2$g2)])
names(seals_g2)[7] <- "IUCN_rating"
seals_g2$IUCN_rating <- factor(seals_g2$IUCN_rating, levels = c("endangered", "vulnerable", "near threatened", "least concern", "data deficient"))


seals_g2[is.na(seals_g2$IUCN_rating), "IUCN_rating"] <- "data deficient"
library(ggplot2)

p3 <- ggplot(seals_g2, aes(x=real_names_species, y=g2,  label = sample_size, color = IUCN_rating)) + 
    geom_errorbar(aes(ymin=CIlow, ymax=CIup), width=.01) +
    geom_point(aes(size = Abundance), alpha = 0.5) +
    scale_size_continuous(breaks=c(1000, 10000, 100000, 1000000, 5000000),
        guide = guide_legend(title = "worldwide abundance")) +
    scale_color_manual(values = c("#d53e4f", "#f46d43", "#fdae61", "#66c2a5", "#000000"),
        guide = guide_legend(title = "IUCN red list rating")) +
    geom_text(angle = 0, vjust = 1.4, hjust = -0.5, size = 3) +
    #scale_x_discrete(name = aes(real_names_species)) +
    theme_classic() +
    theme(axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid'),
        axis.title.y=element_text(margin=margin(0,20,0,0)),
        axis.title.x=element_text(margin=margin(20,0,0,0))) +
    #scale_x_discrete(labels = real_names_species) +
    geom_hline(yintercept = 0) +
    labs(x = "species", y = "variance in inbreeding (identity disequilibrium)") +
    theme(legend.position = c(0.7, 0.5)) +
    coord_flip()  
  
ggplot2::ggsave(p3, 
    filename = "inbreeding.tiff", 
    #  path = "D:/recent R",
    width = 12,
    scale = 0.6,
    height = 12, units = "in",
    dpi = 300)



## diversity
names(seals)
genetic_diversity <- seals %>% select_("species", "mean_het", "Hexp", "Hobs", "TPM70_ratio", "TPM90_ratio") 
gen_div <- melt(genetic_diversity, id.vars = "species")

p <- ggplot(gen_div, aes(x= variable, y = species, fill = value)) + 
    #facet_grid(.~dataset) + 
    geom_tile(color = "white", size = 0.1) +
    labs(x = "microsat mutation model", y = "") +
    # scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#1d91c0", "#225ea8", "#0c2c84", "#081d58")) +
    #scale_fill_gradientn(colours=c("#fff7fb", "#d0d1e6", "#67a9cf", "#02818a", "#014636")) +
    #scale_fill_gradientn(colours=c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#cb181d", "#a50f15", "#67000d")) +
    # scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#41b6c4", "#225ea8", "#253494", "#081d58")) +
    # for p_val
    # scale_fill_gradientn(colours=c("#081d58",  "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "#edf8b1", "#ffffd9"), 
    #     name = "p-val") +
    
    scale_fill_gradientn(colours=c( "#ffffd9","#c7e9b4", "#7fcdbb", "#1d91c0", "#253494", "#081d58"), 
        name = "prop. \nhet-exc") +
    
    scale_x_discrete(labels=c("IAM", "TPM70", "TPM90", "TPM95", "SMM")) +
    scale_y_discrete(labels = rev(seal_names_real)) +
    # scale_fill_gradient(name = "p-value", label = comma,  breaks = c(0.05, 0.5)) +
    theme_tufte(base_family="Helvetica") +
    coord_equal() +
    theme(plot.title=element_text(hjust=0),
        axis.ticks=element_blank(),
        axis.text=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=9), 
        axis.text.x = element_text(angle = 50, hjust = 1))