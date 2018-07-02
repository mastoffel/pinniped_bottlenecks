# additional visualisations

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
