# scatterplot for presentation
library(reshape2)
head(all_stats)
library(ggrepel)
?melt

stats_melted <- all_stats %>% melt(id.vars = c("seal_names_real", "family", "BreedingType", "abc", 
                                              "mean_SSD", "TPM70_ratio", "num_alleles_mean"),
                                   measure.vars = c("harem_size", 
                                            "mratio_mean"))

library(ggthemr)
ggthemr(palette = "fresh", spacing = 3, text_size = 10)

p <- ggplot(data = all_stats, aes(x = mean_SSD, y = TPM70_ratio, fill = BreedingType, size = harem_size)) +
    geom_smooth(method = "lm", size = 0.5, color = "grey", fill = "lightblue", se = FALSE) +
    geom_point(alpha = 0.7, shape = 21)  +
    # geom_point(shape = 21, col = "black") +
    scale_size_continuous(range = c(1.5, 6), breaks = c(1, 5, 10, 30)) +
    scale_fill_manual(values= c("#FDE725FF", "#22A884FF", "#440154FF")) +
    xlab("Sexual Size Dimorphism") +
    ylab("Prop. loci in \nHeterozygosity Excess") +
    guides(fill = guide_legend(title = "Breeding Location"),
        size = guide_legend(title = "Harem Size"))

p <- ggplot(data = all_stats, aes(x = mean_SSD, y = TPM70_ratio, fill = BreedingType, size = harem_size, label = seal_names_real)) +
    # geom_smooth(method = "lm", size = 0.5, color = "grey", fill = "lightblue", se = FALSE) +
    geom_point(alpha = 0.7, shape = 21)  +
    ggrepel::geom_text_repel(size = 2, alpha = 0.3, box.padding = unit(0.5, "lines"), segment.alpha = 0.2,
        segment.size = 0.2) +
   # geom_point(shape = 21, col = "black") +
    scale_size_continuous(range = c(1.5, 6), breaks = c(1, 5, 10, 30)) +
    scale_fill_manual(values= c("#FDE725FF", "#22A884FF", "#440154FF")) +
    xlab("Sexual Size Dimorphism") +
    ylab("Bottleneck strength\n(Prop. loci with \nheterozygosity excess)") +
    guides(fill = guide_legend(title = "Breeding Location"),
           size = guide_legend(title = "Harem Size"))

p
ggsave("scatterplot2.jpg", plot = p, width = 7, height = 5)

ggplot(data = all_stats, aes(y = TPM70_ratio, x = BreedingType)) +
    geom_boxplot() 


mod <- lm(data = all_stats, formula = TPM70_ratio~mean_SSD + BreedingType + harem_size)
mod <- lm(data = all_stats, formula = TPM70_ratio~mean_SSD)
summary(mod)
