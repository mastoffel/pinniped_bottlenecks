# plot for visualizing abc

library(ggplot2)
library(ggthemr)
library(ggExtra)
library(cowplot)
ggthemr(palette = "fresh", layout = "clean", spacing = 3)

set.seed(1)
# Generate the data 
n = 500   	  # Sample size
alpha = 3   # Simulated (i.e. "true") causal effect of predictor on response
beta1 = 0.7   # Simulated (i.e. "true") causal effect of predictor on response
pred1 = rnorm(n, mean=10, sd=1)
epsilon = rnorm(n, mean=0, sd=0.6)
resp = alpha + ((pred1-10) * beta1) + epsilon
md = data.frame(resp, pred1)
plot(md$resp, md$pred1)

md$col <- ifelse(md$pred1 > 10.5 & md$pred1 < 12.5, "tol", "nottol")

# ggplot

ggplot(md, aes(x = pred1, y = resp)) +
    geom_point(size = 3, alpha = 0.5) +
    ylim(0, 7) +
    xlab("summary statistic \n(allele number, M-ratio)") +
    ylab("model parameter \n(Ne, mu)") +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = "abc1.jpg", width = 6, height = 4)

ggplot(md, aes(x = pred1, y = resp)) +
    geom_point(size = 3, alpha = 0.5) +
    xlab("summary statistic \n(allele number, M-ratio)") +
    ylab("model parameter \n(Ne, mu)") +
    ylim(0, 7) +
    #geom_vline(xintercept = 11.5,  linetype = "longdash") +
    geom_segment(aes(x = 11.5, y = 0, xend = 11.5, yend = 6), linetype = "longdash", size = 0.3,
        colour = "goldenrod") +
    annotate("text", x = 11.5, y = 6.2, label = "atop(bold(S(emp)))", parse  =TRUE, 
        colour = "goldenrod")  +
    theme(plot.margin=unit(c(1,1,1,1),"cm"))
ggsave(filename = "abc2.jpg", width = 6, height = 4)

p1 <- ggplot(md, aes(x = pred1, y = resp, colour = factor(col))) +
    geom_point(size = 3, alpha = 0.5) +
    ylim(0, 7) +
    xlab("summary statistic \n(allele number, M-ratio)") +
    ylab("model parameter \n(Ne, mu)") + 
    #geom_vline(xintercept = 11.5,  linetype = "longdash") +
    geom_segment(aes(x = 11.5, y = 0, xend = 11.5, yend = 6), linetype = "longdash", size = 0.3,
        colour = "black") +
    annotate("text", x = 11.5, y = 6.2, label = "atop(bold(S(emp)))", parse  =TRUE, 
        colour = "black") +
    geom_segment(aes(x = 11.5, y = 0, xend = 12.5, yend = 0), size = 0.3,
        colour = "goldenrod", arrow = arrow(length = unit(0.3, "cm"))) +
    geom_segment(aes(x = 11.5, y = 0, xend = 10.5, yend = 0), size = 0.3,
        colour = "goldenrod", arrow = arrow(length = unit(0.3, "cm"))) +
    annotate("text", x = 11, y = 0, label = "atop(bold(-tol))", parse  =TRUE, 
        colour = "goldenrod")  +
    annotate("text", x = 12, y = 0, label = "atop(bold(+tol))", parse  =TRUE, 
        colour = "goldenrod")  +
    geom_segment(aes(x = 12.5, y = 0, xend = 12.5, yend = 7), size = 0.3,
        colour = "goldenrod", alpha = 0.5) +
    geom_segment(aes(x = 10.5, y = 0, xend = 10.5, yend = 7), size = 0.3,
        colour = "goldenrod", alpha = 0.5) +
    scale_colour_manual(values = c("grey", "goldenrod")) +
    theme(legend.position="none") +
    theme(plot.margin=unit(c(1,0,1,1),"cm"))
ggsave(filename = "abc3.jpg", width = 6, height = 4)

md_hist <- md[md$col == "tol", ]
p2 <- ggplot(md_hist, aes(x = resp)) +
    geom_histogram(bins = 30, fill = "goldenrod", colour = "black", alpha = 0.8,
        size = 0.2) +
    coord_flip() +
    ylab("summary statistic \n(allele number, M-ratio)") +
    no_axes() +
    xlim(0, 7) +
    theme(plot.margin=unit(c(0.6,1,2.3,0),"cm"))

plot_grid(p1, p2, rel_widths = c(2,1))
ggsave(filename = "abc4.jpg", width = 7.3, height = 4)




## Last plot on black background

p1 <- ggplot(md, aes(x = pred1, y = resp, colour = factor(col))) +
    geom_point(size = 3, alpha = 0.5) +
    ylim(0, 7) +
    xlab("summary statistic \n(allele number, M-ratio)") +
    ylab("model parameter \n(Ne, mu)") + 
    #geom_vline(xintercept = 11.5,  linetype = "longdash") +
    geom_segment(aes(x = 11.5, y = 0, xend = 11.5, yend = 6), linetype = "longdash", size = 0.3,
        colour = "cornflowerblue") +
    annotate("text", x = 11.5, y = 6.2, label = "atop(bold(S(emp)))", parse  =TRUE, 
        colour = "cornflowerblue") +
    geom_segment(aes(x = 11.5, y = 0, xend = 12.5, yend = 0), size = 0.3,
        colour = "goldenrod", arrow = arrow(length = unit(0.3, "cm"))) +
    geom_segment(aes(x = 11.5, y = 0, xend = 10.5, yend = 0), size = 0.3,
        colour = "goldenrod", arrow = arrow(length = unit(0.3, "cm"))) +
    annotate("text", x = 11, y = 0, label = "atop(bold(-tol))", parse  =TRUE, 
        colour = "goldenrod")  +
    annotate("text", x = 12, y = 0, label = "atop(bold(+tol))", parse  =TRUE, 
        colour = "goldenrod")  +
    geom_segment(aes(x = 12.5, y = 0, xend = 12.5, yend = 7), size = 0.3,
        colour = "goldenrod", alpha = 0.5) +
    geom_segment(aes(x = 10.5, y = 0, xend = 10.5, yend = 7), size = 0.3,
        colour = "goldenrod", alpha = 0.5) +
    scale_colour_manual(values = c("grey", "goldenrod")) +
    theme(plot.margin=unit(c(1,0,1,1),"cm"),
        legend.position="none",
          plot.background = element_rect(fill = 'black', colour = 'black'),
        panel.background = element_rect(fill = 'black', colour = 'black'),
          axis.text = element_text(colour = "white"),
          axis.line = element_line(colour = "white"),
          axis.ticks = element_line(colour = "white"),
          axis.title =  element_text(colour = "white"))
p1
#ggsave(filename = "abc3.jpg", width = 6, height = 4)

md_hist <- md[md$col == "tol", ]
p2 <- ggplot(md_hist, aes(x = resp)) +
    geom_histogram(bins = 30, fill = "goldenrod", colour = "black", alpha = 0.8,
        size = 0.2) +
    coord_flip() +
    ylab("summary statistic \n(allele number, M-ratio)") +
    no_axes() +
    xlim(0, 7) +
    theme(plot.margin=unit(c(0.6,1,2.3,0),"cm"),
        plot.background = element_rect(fill = 'black', colour = 'black'),
        panel.background = element_rect(fill = 'black', colour = 'black'))

plot_grid(p1, p2, rel_widths = c(2,1))
#ggsave(filename = "abc4.jpg", width = 7.3, height = 4)
