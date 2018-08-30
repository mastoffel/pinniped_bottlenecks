library(abc)
library(readr)
library(dplyr)
library(purrr)
load("output/model_evaluation/check3_modeval/sims_gof_bot.RData")
load("output/model_evaluation/check3_modeval/sims_gof_neut.RData")

# load stats to know sequence of species
all_stats <- read_delim("data/processed/all_sumstats_40ind_30_server_sequence.txt", delim = " ") 

all_stats_full <- read_csv("data/processed/all_stats_30_modeling.csv") %>% 
    dplyr::select(species, common, bot)

# load goodness-of-fits
all_gfit_bot <- unlist(map(sims_10000kbot500good_of_fit_bot, function(x) x[[1]]))
all_gfit_neut <- unlist(map(sims_10000kbot500good_of_fit_neut, function(x) x[[1]]))
all_gfit <- c(all_gfit_bot, all_gfit_neut)

all_gfit_df <- data.frame(species = rep(rep(all_stats$species, each = length(all_gfit)/2 / 30), 2), gfits = all_gfit) %>% 
                    mutate(mod = rep(c("bot", "neut"), each = length(all_gfit_bot)))

# load observed value
obs_gfit_bot <- unlist(map(sims_10000kbot500good_of_fit_bot, function(x) x[[2]]))
obs_gfit_neut <- unlist(map(sims_10000kbot500good_of_fit_neut, function(x) x[[2]]))
obs_gfit <- c(obs_gfit_bot, obs_gfit_neut)

obs_gfit_df <- data.frame(species = rep(rep(all_stats$species, each = length(obs_gfit)/2 / 30), 2), gfits = obs_gfit) %>% 
    mutate(mod = rep(c("bot", "neut"), each = length(obs_gfit_bot)))

obs_gfit_plot <- obs_gfit_df %>% 
    left_join(all_stats_full, by = "species") %>% 
    mutate(bot = ifelse(bot > 0.5, "bot", "neut")) %>%
    filter(mod == bot)


# get bottleneck model goodness of fit for bottlenecked species and neutral model gof for neutral species
gfit_plot <- all_gfit_df %>% 
                left_join(all_stats_full, by = "species") %>% 
                mutate(bot = ifelse(bot > 0.5, "bot", "neut")) %>%
                filter(mod == bot)
source("R/martin.R")
p <- ggplot(gfit_plot, aes(gfits)) +
    geom_histogram(aes(fill = mod), bins = 15) +
    geom_vline(aes(xintercept = gfits ), obs_gfit_plot) +
    facet_wrap(~ common, scales = "free") +
    theme_martin() +
    scale_fill_manual(values = c("#fc8d62","#8da0cb")) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle = 0),
        axis.text= element_text(size = 9.5),
        axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.title=element_blank(),
        axis.line.x = element_line(color="black", size = 0.5)) 
p
