# create prior distribution plots for schematic

library(ggplot2)
source("R/martin.R")
ne <- data.frame(vals = rlnorm(100000, 10.5, 1))

ggplot(ne, aes(vals)) +
    # stat_density(geom="line", adjust = 3, fill = "orange", alpha = 0.7) +
    geom_density(adjust = 3, fill = "#f68e25", alpha = 1, size = 0.1) +
    theme_martin() +
    theme(axis.line = element_line(color="black", size = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank()) +
    expand_limits(x = 0, y = 0) +
  #  xlab(expression(N[e])) +
    xlab("") +
    ylab("") +
    #scale_x_continuous(expand = c(0, 0)) +
    xlim(c(0, 500000)) 

ggsave(filename = "../bottleneck/other_stuff/figures/figures_final/new_figures_revision2/Ne.jpg",
       width = 4, height = 3)


uniform_ne <- data.frame(vals = rep(1:5, each = 100000), x = 1:500000)

ne <- data.frame(vals = rlnorm(100, 10.5, 1))

ggplot(ne, aes(vals)) +
    # stat_density(geom="line", adjust = 3, fill = "orange", alpha = 0.7) +
   # geom_abline(intercept = 0, slope = 0, color = "orange", alpha = 0.7) +
    theme_martin() +
    ylim(c(-1, 1)) +
   # geom_ribbon(aes(ymin=-1, ymax=0), fill = "#f68e25", alpha = 0.2) +
    geom_ribbon(aes(ymin=-0.01, ymax=0.15), fill = "#f68e25") +
    theme(axis.line = element_line(color="black", size = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank()) +
    expand_limits(x = 0, y = 0) +
    #  xlab(expression(N[e])) +
    xlab("") +
    ylab("") +
    #scale_x_continuous(expand = c(0, 0)) +
    xlim(c(0, 500000)) 

ggsave(filename = "../bottleneck/other_stuff/figures/figures_final/new_figures_revision2/Ne_unif.jpg",
    width = 4, height = 3)





ggplot(ne, aes(vals)) +
    # stat_density(geom="line", adjust = 3, fill = "orange", alpha = 0.7) +
    geom_density(adjust = 3, fill = "#f68e25", alpha = 1, size = 0.1) +
    theme_martin() +
    theme(axis.line = element_line(color="black", size = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        rect = element_rect(fill = "transparent")) +
    expand_limits(x = 0, y = 0) +
    #  xlab(expression(N[e])) +
    xlab("") +
    ylab("") +
    #scale_x_continuous(expand = c(0, 0)) +
    xlim(c(0, 500000)) 
ggsave(filename = "../bottleneck/other_stuff/figures/figures_final/new_figures_revision2/Ne_nobg.png",
    width = 4, height = 3, bg = "transparent")



ggplot(ne, aes(vals)) +
    # stat_density(geom="line", adjust = 3, fill = "orange", alpha = 0.7) +
    # geom_abline(intercept = 0, slope = 0, color = "orange", alpha = 0.7) +
    theme_martin() +
    ylim(c(-1, 1)) +
    # geom_ribbon(aes(ymin=-1, ymax=0), fill = "#f68e25", alpha = 0.2) +
    geom_ribbon(aes(ymin=-0.01, ymax=0.15), fill = "#f68e25") +
    theme(axis.line = element_line(color="black", size = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        rect = element_rect(fill = "transparent")) +
    expand_limits(x = 0, y = 0) +
    #  xlab(expression(N[e])) +
    xlab("") +
    ylab("") +
    #scale_x_continuous(expand = c(0, 0)) +
    xlim(c(0, 500000)) 

ggsave(filename = "../bottleneck/other_stuff/figures/figures_final/new_figures_revision2/Ne_unif.png",
    width = 4, height = 3, bg = "transparent")




# Ne LGM

ne <- data.frame(vals = rlnorm(100000, 9, 1))

ggplot(ne, aes(vals)) +
    # stat_density(geom="line", adjust = 3, fill = "orange", alpha = 0.7) +
    geom_density(adjust = 3, fill = "#f68e25", alpha = 1, size = 0.1) +
    theme_martin() +
    theme(axis.line = element_line(color="black", size = 2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        rect = element_rect(fill = "transparent")) +
    expand_limits(x = 0, y = 0) +
    #  xlab(expression(N[e])) +
    xlab("") +
    ylab("") +
    #scale_x_continuous(expand = c(0, 0)) +
    xlim(c(0, 500000)) 

ggsave(filename = "../bottleneck/other_stuff/figures/figures_final/new_figures_revision2/Ne_LGM.png",
    width = 4, height = 3, bg = "transparent")















df <- data.frame()
ggplot(df) +
    # stat_density(geom="line", adjust = 3, fill = "orange", alpha = 0.7) +
   # geom_density(fill = "#f68e25") +
    theme_martin() +
    theme(axis.line = element_line(color="black", size = 0.6),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text = element_blank()) +
    geom_abline(intercept = 0, slope = 0) +
    expand_limits(x = 0, y = 0) +
    #  xlab(expression(N[e])) +
    xlab("") +
    ylab("") +
    #scale_x_continuous(expand = c(0, 0)) +
    xlim(c(0, 500000)) 

