p <- ggtree(tree_final, aes(color=group)) + #layout="circular" , "fan"open.angle=180
    ggtitle("Pinniped Phylogeny") +
    geom_tiplab(size=2, color="black") + # tiplap 2 for circular
    # ggplot2::xlim(0, 40) +
    # geom_cladelabel(node=41, label="Phocidae", angle = 270, align=T, hjust='center', offset = 5, 
    #     offset.text=1, color = "cornflowerblue",  barsize=0.5) +
    # geom_cladelabel(node=31, label="Otariidae", align=T,angle = 270, hjust='center', offset = 5,
    #     offset.text=1, barsize=0.5, color = "goldenrod") +
    # geom_cladelabel(node=1, label="Odobenidae", align=T,angle = 270, color = "darkgrey", hjust=0.7, 
    #     offset.text=6, barsize=0.5) +
    scale_color_manual(values=c("black","darkgrey", "goldenrod", "cornflowerblue" ))

stand_div_lf <- stand_div %>% melt(id.vars = "id")
p2 <- facet_plot(p + xlim_tree(20), panel = "diversity", data = stand_div_lf, geom_tile, 
    aes(x= as.numeric(as.factor(variable)), fill = value)) +  
    scale_fill_viridis(option = "inferno") +
    # scale_fill_distiller(palette = "BuPu") + 
    theme_tree2()
#xlim_expand(c(0,1), panel = "heatmap")

bot_res_lf <- bot_res %>% melt(id.vars = "id") %>% mutate(.panel='bottleneck')
p3 <- facet_plot(p2 + xlim_tree(20), panel = "bottleneck", data = bot_res_lf, geom_tile, 
    aes(x= as.numeric(as.factor(variable)), fill = value)) +  
    theme_tree2() 
p3