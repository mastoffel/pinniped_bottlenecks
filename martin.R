theme_martin <- function() {
    
    theme_minimal(base_size = 12,  base_family='Lato Light') +        #%+replace% 
        theme(
            panel.background  = element_blank(),
            text = element_text(color='#333333') # color='#333333' 
        )
    
}
