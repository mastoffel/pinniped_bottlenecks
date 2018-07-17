# evaluate post.predictive checks

library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
# simulated based on posteriors
all_checks <- read_delim("output/model_evaluation/check5_postpred/sims_10000kbot500_post_pred_checks2.txt", delim = " ")
head(all_checks)

# compare these summary statistics
sumstats <- c("num_alleles_mean", 
    "exp_het_mean", 
    "mratio_mean", 
    "prop_low_afs_mean",
    "mean_allele_range")

# empirical summary stats
all_stats <- read_csv("data/processed/all_stats_30_modeling.csv") %>% 
    dplyr::select(species, common, sumstats, bot)

# long format for plotting
all_checks_long <- all_checks %>% 
    left_join(all_stats[c("species", "common")], by = "species") %>% 
    dplyr::select(common, sumstats) %>% 
    gather(sumstat, value, -common) 

# observed sumstats long format
all_sumstats_full_long <- all_stats %>% 
    gather(sumstat, value, -species, -common, -bot)

# lookup_table <- paste(all_stats$species, " = ", all_stats$common)

sumstat_names <- c(
    exp_het_mean = "Expected\nhetetozygosity",
    mean_allele_range = "Allelic range",
    mratio_mean = "M-ratio",
    num_alleles_mean = "Allelic richness",
    prop_low_afs_mean  = "Prop. of low\nfrequency alleles"
)

# bottlenecked or not bottlenecked
all_stats$mod <- ifelse(all_stats$bot > 0.5, "Bottleneck", "Non-bottleneck")

all_data <- all_checks_long %>% 
    left_join(all_stats[c("common", "mod")])

p <- ggplot(all_data, aes(value)) +
    geom_histogram(aes(fill = mod)) +
    geom_vline(aes(xintercept = value), all_sumstats_full_long) +
    facet_grid(common ~ sumstat, scales = "free", labeller = labeller(
        sumstat = sumstat_names
    )) +
    theme_martin() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0),
        axis.text= element_text(size = 7),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.title=element_blank(),
        axis.line.x = element_line(color="black", size = 0.5)) 
    


ggsave(filename = "post_pred_checks.jpg", plot = p, width = 10, height = 10)

ggplot(all_checks, aes(value)) +
    geom_histogram() +
    facet_wrap(species ~ sumstat)