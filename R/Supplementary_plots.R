# Making additional plots for the supplementary


# Confusion Matrix
# Supplementary figure 1
library(forcats)
conf_mat <- data.frame("Bottleneck" = c(15, 85), "Neutral" = c(89, 11)) %>% 
    gather(key = `Classified into`, value = Frequency) %>% 
    mutate(`Simulated under` = c("Neutral", "Bottleneck", "Neutral", "Bottleneck")) %>% 
    mutate(`Simulated under` = fct_relevel(`Simulated under`, c("Neutral", "Bottleneck")))

p <- ggplot(conf_mat, aes(x=`Classified into`, y=Frequency, fill = `Simulated under`)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = seq(from = 0, to = 100, by = 20)) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    scale_fill_manual(values = c("grey", "gray27"))

ggsave(filename = "other_stuff/figures/figures_final/figures_revised_manuscript/Sup_fig_conf_mat.jpg",
    height = 3, width = 4.5)



# Ne bot cross-validation

library(ggplot2)
source("R/martin.R")
load("output/model_evaluation/check4_params/cv_param_it100_parallel_rej_sims_10000kbot500_bot_30.RData")

all_cv_df <- data.frame("estimate" = all_cv$estim[[1]]$nbot, "true" = all_cv$true$nbot)

p <- ggplot(all_cv_df, aes(true, estimate)) +
    geom_point(size = 2.3, alpha = 0.3) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    geom_line(data = data.frame(x = 1:500, y = 1:500), mapping = aes(x, y)) +
    scale_y_continuous(limits = c(0, 400)) +
    xlab("True value") +
    ylab("Parameter estimate")

ggsave(p, file = "other_stuff/figures/figures_final/figures_revised_manuscript/Sup_CV_plot.jpg",
       width = 4.5, height = 4)


# Confusion Matrix for 4 models
# Supplementary figure 1
library(forcats)
library(readr)
library(dplyr)
library(tidyr)
conf_mat <- read_delim("output/model_evaluation/check2_models/sims_1000k_4mods_expand_model_missclass_freq.txt",
    delim = " ")
names(conf_mat) <- c("var1", "var2", "freq")

library(ggplot2)
source("R/martin.R")
p <- ggplot(conf_mat, aes(x=var1, y= freq, fill = var2)) +
    geom_bar(stat = "identity") +
    scale_y_continuous(breaks = seq(from = 0, to = 100, by = 20)) +
    scale_x_discrete(labels = c("LGM +\nBottleneck ", "Bottleneck", "LGM +\nNeutral", "Neutral")) +
    theme_martin(base_family = "Hind Guntur Light", highlight_family = "Hind Guntur Light") +
    xlab("Simulated under") +
    ylab("Frequency") +
    #scale_fill_viridis_d(name = "Classified into", option = "viridis",
    #    labels = c("Bottleneck +\nIce Age ", "Bottleneck", "Neutral+\nIce Age", "Neutral"))
    scale_fill_manual(name = "Classified into", 
        values = rev(c("#d9d9d9", "#969696", "#525252", "#252525")),
        labels = c("LGM +\nBottleneck ", "Bottleneck", "LGM +\nNeutral", "Neutral"))
p
ggsave(filename = "other_stuff/figures/figures_final/figures_revised_manuscript/Sup_fig_conf_mat_4mod_expand.jpg",
    height = 3, width = 5.5)
