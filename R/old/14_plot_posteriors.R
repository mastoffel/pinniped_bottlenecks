library(ggplot2)
library(ggthemr)
library(cowplot)
library(dplyr)


ggthemr(palette = "pale", layout = "clean", spacing = 1, line_weight = 1, text_size = 14)

# load all posteriors
abc_full <- read.table("abc_full_3k.txt", header = TRUE)
# subset important species
abc_sub <- abc_full %>% 
                filter(species %in% c("galapagos_fur_seal" , "south_american_fur_seal",
                                    "antarctic_fur_seal" , "new_zealand_fur_seal",
                                     "australian_fur_seal", "hawaiian_monk_seal",
                                     "mediterranean_monk_seal" , "nes" , "saimaa_ringed_seal",
                                     "lagoda_ringed_seal", "baltic_ringed_seal")) 

 abc_sub <- abc_full %>% 
     filter(species %in% c("antarctic_fur_seal",
         "mediterranean_monk_seal" , "nes",  "hawaiian_monk_seal",
         "lagoda_ringed_seal")) 

seal_names <- c(
    `antarctic_fur_seal` = "Antarctic \nFur Seal",
    `hawaiian_monk_seal` = "Hawaiian \nMonk Seal",
    `mediterranean_monk_seal` = "Mediterranean \nMonk Seal",
    `nes` = "Northern \nElephant Seal",
    `lagoda_ringed_seal` = "Ladoga \nRinged Seal"
)

# subset parameter
abc_sub_Nbot <- abc_sub %>% filter(pars %in% c("N_bot"))

p1 <- ggplot(abc_sub_Nbot, aes(posterior)) +
    geom_histogram(aes(y=..density..), position="identity", 
        colour = "cornflowerblue",fill="white", alpha=0.3, binwidth = 3) +
    geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
    facet_wrap(~species, ncol = 1, scales = "free", labeller = as_labeller(seal_names),
        strip.position = c("left")) +
    xlim(c(0,300)) +
    xlab("N (bottleneck)") +
    ylab("") +
    theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text.y = element_text( angle = 180))


# subset parameter
abc_sub_start_bot <- abc_sub %>% filter(pars %in% c("start_bot"))

p2 <- ggplot(abc_sub_start_bot, aes(posterior)) +
    geom_histogram(aes(y=..density..), position="identity", 
        colour = "cornflowerblue",fill="white", alpha=0.3, binwidth = 0.3) +
    geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
    facet_wrap(~species, ncol = 1, scales = "free",
        strip.position = c("left")
    ) +
    #xlim(c(0,300)) +
    xlab("start bottleneck (gen. ago)") +
    ylab("") +
    theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text.y = element_blank())

# subset parameter
abc_sub_end_bot <- abc_sub %>% filter(pars %in% c("end_bot"))

p3 <- ggplot(abc_sub_end_bot, aes(posterior)) +
    geom_histogram(aes(y=..density..), position="identity", 
        colour = "cornflowerblue",fill="white", alpha=0.3, binwidth = 0.3) +
    geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
    facet_wrap(~species, ncol = 1, scales = "free",
        strip.position = c("left")
    ) +
    #xlim(c(0,300)) +
    xlab("end bottleneck") +
    ylab("") +
    theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text.y = element_blank())

plot_grid(p1,p2,p3, ncol = 3, rel_widths = c(1.7, 1, 1))

ggsave(
    filename = "N_bot_all.jpg",
    width = 10.2,
    height = 11, units = "in",
    dpi = 300)





## plotting mutation posteriors

# subset parameter
abc_sub_mu <- abc_sub %>% filter(pars %in% c("mu"))

p1 <- ggplot(abc_sub_mu, aes(posterior)) +
    geom_histogram(aes(y=..density..), position="identity", 
        colour = "cornflowerblue",fill="white", alpha=0.3, binwidth = 0.0001) +
    geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
    facet_wrap(~species, ncol = 1, scales = "free", labeller = as_labeller(seal_names),
        strip.position = c("left")) +
    xlim(c(0,0.002)) +
    xlab("mutation rate") +
    ylab("") +
    theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text.y = element_text( angle = 180))

# subset parameter
abc_sub_p_single <- abc_sub %>% filter(pars %in% c("p_single"))

p2 <- ggplot(abc_sub_p_single, aes(posterior)) +
    geom_histogram(aes(y=..density..), position="identity", 
        colour = "cornflowerblue",fill="white", alpha=0.3, binwidth = 0.01) +
    geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
    facet_wrap(~species, ncol = 1, scales = "free",
        strip.position = c("left")
    ) +
    #xlim(c(0,300)) +
    xlab("prop. multistep mutations") +
    ylab("") +
    theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        strip.text.y = element_blank())














# subset parameter
abc_sub_Nbot <- abc_sub %>% filter(pars %in% c("N_bot"))

 p1 <- ggplot(abc_sub_Nbot, aes(posterior)) +
    geom_histogram(aes(y=..density..), position="identity", 
        colour = "cornflowerblue",fill="white", alpha=0.3) +
    geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
    xlim(c(0,100)) +
    facet_grid(.~species) +
    theme(plot.margin=unit(c(20,10,10,10),"points"),
        strip.text.x = element_blank()) +
    xlab("N (bottleneck)") +
    ylab("")

 ggsave(p1,
     filename = "N_bot.jpg",
     width = 10,
     height = 4, units = "in",
     dpi = 300)
    

 # subset parameter
 abc_sub_startbot <- abc_sub %>% filter(pars %in% c("start_bot"))
 
 p2 <- ggplot(abc_sub_startbot, aes(posterior)) +
     geom_histogram(aes(y=..density..), position="identity", 
         colour = "cornflowerblue",fill="white", alpha=0.3) +
     geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
     #xlim(c(0,800)) +
     facet_grid(.~species) +
     theme(plot.margin=unit(c(20,10,10,10),"points"),
         strip.text.x = element_blank()) +
     xlab("start bottleneck (in generations)") +
     ylab("")
 
 ggsave(p2,
     filename = "start_bot.jpg",
     width = 10,
     height = 4, units = "in",
     dpi = 300)
 
 # subset parameter
 abc_sub_endbot <- abc_sub %>% filter(pars %in% c("end_bot"))
 
 p3 <- ggplot(abc_sub_endbot , aes(posterior)) +
     geom_histogram(aes(y=..density..), position="identity", 
         colour = "cornflowerblue",fill="white", alpha=0.3) +
     geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
     #xlim(c(0,800)) +
     facet_grid(.~species) +
     theme(plot.margin=unit(c(20,10,10,10),"points"),
         strip.text.x = element_blank()) +
     xlab("end bottleneck (in generations)") +
     ylab("")
 
 ggsave(p3,
     filename = "end_bot.jpg",
     width = 10,
     height = 4, units = "in",
     dpi = 300)
 
 
 # subset parameter
 abc_sub_p_single <- abc_sub %>% filter(pars %in% c("p_single"))
 
 p4 <- ggplot(abc_sub_p_single , aes(posterior)) +
     geom_histogram(aes(y=..density..), position="identity", 
         colour = "cornflowerblue",fill="white", alpha=0.3) +
     geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
     #xlim(c(0,800)) +
     facet_grid(.~species, scales = "free") +
     theme(plot.margin=unit(c(20,10,10,10),"points"),
         strip.text.x = element_blank()) +
     xlab("proportion single step mutations") +
     ylab("")
 
 ggsave(p4,
     filename = "single_p.jpg",
     width = 10,
     height = 4, units = "in",
     dpi = 300)
 
 
 # subset parameter
 abc_sub_mu <- abc_sub %>% filter(pars %in% c("mu"))
 
 p5 <- ggplot(abc_sub_mu , aes(posterior)) +
     geom_histogram(aes(y=..density..), position="identity", 
         colour = "cornflowerblue",fill="white", alpha=0.3) +
     geom_density(fill = "cornflowerblue", color="black", alpha = 0.1) +
     #xlim(c(0,800)) +
     facet_grid(.~species, scales = "free") +
     theme(plot.margin=unit(c(20,10,10,10),"points"),
         strip.text.x = element_blank()) +
     xlab("microsatellite mutation rate") +
     ylab("")
 
 ggsave(p5,
     filename = "mutationrate.jpg",
     width = 10,
     height = 4, units = "in",
     dpi = 300)
 
 
p6 <-  plot_grid(p1, p2, p3, p4, p5, ncol = 1)
ggsave(p6,
    filename = "full_abc.jpg",
    width = 12,
    height = 30, units = "in",
    dpi = 300)
 