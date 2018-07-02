# make a nice plot if nbot cross-validation

library(ggplot2)
source("R/martin.R")
load("output/model_evaluation/check4_params/cv_param_it100_parallel_rej_sims_10000kbot500_bot_30.RData")

all_cv_df <- data.frame("estimate" = all_cv$estim[[1]]$nbot, "true" = all_cv$true$nbot)

p <- ggplot(all_cv_df, aes(true, estimate)) +
    geom_point(size = 2.3, alpha = 0.3) +
    theme_martin() +
    geom_line(data = data.frame(x = 1:500, y = 1:500), mapping = aes(x, y)) +
    scale_y_continuous(limits = c(0, 400))

ggsave(p, file = "other_stuff/figures/figures_final/figures_revised_manuscript/Sup_CV_plot.jpg",
       width = 4.5, height = 4)
