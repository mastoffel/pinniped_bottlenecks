# visualise output from 4 model analysis

library(readr)
library(tidyr)
library(abc)
library(readxl)
library(dplyr)
library(ggplot2)
library(viridis)
source("R/martin.R")
mod4_ms <- read_delim(file = "data/processed/sims_1000k_4mods_expand_model_selection_30.txt", delim = " ",
           col_names = c("species", "bot_iceage", "bot", "neut_iceage", "neut"), skip = 1)

mod4_ms_lf <- gather(mod4_ms, model, probability, bot_iceage:neut)

?melt
ggplot(mod4_ms_lf, aes(x = model, y = species, fill = probability)) +
    geom_tile(color = "grey", size = 0.1) +
        scale_fill_viridis(name = "Model probabilities", direction = -1) +
    theme_martin()

# table
library(knitr)
library(kableExtra)
all_stats_origin <- read_xlsx("data/raw/table_data.xlsx")
all_stats_table <- all_stats_origin %>% 
    dplyr::select("species", "common", "latin") %>% 
    left_join(mod4_ms, by = "species") %>% 
    dplyr::select(-species) %>% 
    dplyr::rename(`Common name` = common,
        `Scientific name` = latin,
        `LGM + Bottleneck` =  bot_iceage,
        `Bottleneck` = bot,
        `LGM + Non-bottleneck` = neut_iceage,
        `Non-bottleneck` = neut) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(`Scientific name` = cell_spec(`Scientific name`, format = "latex", italic = TRUE))

#options(knitr.table.format = "latex")

align_tab1 <- c("l", "l", rep(x = "c", ncol(all_stats_table) - 2))

kable(all_stats_table , format = "latex",  escape = F, 
    booktabs = TRUE, align = align_tab1, digits = 3, linesep = "") %>% 
    kable_styling(latex_options = c( "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
   #add_header_above(c(" " = 2, "Conservation, demography, ecology and life-history data" = 6,
    #    "Genetic data" = 5), italic = TRUE) %>% 
    kable_as_image("other_stuff/figures/figures_final/new_figures_revision_2/supplementary_figures_4mod/model_probs_expand.png", keep_pdf = TRUE)


# new plot
# mod4_ms <- read_delim(file = "data/processed/sims_1000k_4mods_model_selection_30.txt", delim = " ",
#     col_names = c("species", "bot_iceage", "bot", "neut_iceage", "neut"), skip = 1)

# transform into logical
mod4_ms_logi <- mod4_ms
mod4_ms_logi[2:5] <- data.frame(t(apply(mod4_ms[2:5], 1, function(x) max(x) == x)))

mod2_ms <- read_delim(file = "data/processed/sims_10000kbot500_model_selection_30.txt", delim = " ",
    col_names = c("species", "bot2", "neut2"), skip = 1)
mod2_ms_logi <- mod2_ms
mod2_ms_logi[2:3] <- data.frame(t(apply(mod2_ms[2:3], 1, function(x) max(x) == x)))

mod_ms_all <- mod4_ms_logi %>% 
    left_join(mod2_ms_logi, by = "species")

n_bot <- sum(mod_ms_all$bot2)
n_neut <- sum(mod_ms_all$neut2)

mod_ms_all[2:ncol(mod_ms_all)] <- lapply(mod_ms_all[2:ncol(mod_ms_all)], as.logical)

# which were bottlenecked in both scenarios
n_bot_bot <- apply(mod_ms_all[2:ncol(mod_ms_all)], 1, function(x) (x["bot_iceage"] | x["bot"]) & x["bot2"])
sum(n_bot_bot)

# which favored the recent bottleneck scenario in both models ## 8
#n_bot_bot <- apply(mod_ms_all[2:ncol(mod_ms_all)], 1, function(x) (x["bot"]) & x["bot2"])
#sum(n_bot_bot)

# which were neutral in both  scenarios
n_bot_bot <- apply(mod_ms_all[2:ncol(mod_ms_all)], 1, function(x) (x["neut_iceage"] | x["neut"]) & x["neut2"])
sum(n_bot_bot)




# how is the new classification
n_bot_bot <- apply(mod_ms_all[2:ncol(mod_ms_all)], 1, function(x) (x["bot"]) & x["bot2"])
sum(n_bot_bot)
n_iceagebot_bot <- apply(mod_ms_all[2:ncol(mod_ms_all)], 1, function(x) (x["bot_iceage"]) & x["bot2"])
sum(n_iceagebot_bot)
n_neut_neut <- apply(mod_ms_all[2:ncol(mod_ms_all)], 1, function(x) (x["neut"]) & x["neut2"])
sum(n_neut_neut)

# classification
bot_original <- c("neutral" = 0, "iceage_neutral" = 0, "iceage_bottleneck" = 3, "bottleneck" = 8)
neut_original <- c("neutral" = 4, "iceage_neutral" = 9, "iceage_bottleneck" = 6, "bottleneck" = 0)

library(tibble)
class_df <- data.frame(bot_original, neut_original) %>% 
                rownames_to_column("new_model") %>% 
                gather(key = model, value = freq, -new_model)

library(ggplot2)
ggplot(class_df, aes(model, freq, color = new_model)) + 
    geom_point()

ggplot(class_df, aes(x=model, y=freq, fill = new_model)) + 
    geom_bar(stat="identity", position = "dodge", width = 1.3) +
    theme_martin() +
    scale_fill_manual(name = "New model classification", 
        values= c("#a6611a", "#dfc27d", "#80cdc1", "#018571"),
        labels = c("Bottleneck", "Ice age + Bottleneck", "Ice age + Neutral", "Neutral")) +
    xlab("Original model classification") +
    ylab("Frequency") +
    scale_x_discrete(labels = c("Bottleneck", "Neutral")) +
    scale_y_continuous(breaks = c(seq(from = 0, to = 9, by = 2))) 
    

ggsave(filename = "other_stuff/figures/figures_final/figures_4mod/classification_plot.jpg",
    width = 5, height = 2.6)


