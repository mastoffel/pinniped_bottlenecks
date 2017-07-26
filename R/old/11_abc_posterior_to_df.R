# plotting posteriors from abc
library(abc)

## prepare the posteriors for plotting


# put all adjusted values for all parameters and species into one data.frame
posterior_to_df <- function(pop_size){
    # load data
    to_load <- paste0("data/processed/abc_estimates/0_001/abc_", pop_size, ".RData")
    abcdata_name <- load(to_load )
    abcdata <- get(abcdata_name)
    rm(list=setdiff(ls(), c("abcdata", "pop_size")))
    
    # extract all posteriors into data.frame
    parameter_df <- abcdata[[1]]
    abc_res <- abcdata[[2]]
    
    # just take neural net results
    neuralnet <- which(parameter_df$methods == "neuralnet")
    
    extract_abc <- function(ind){
        df <- data.frame(parameter_df[ind, ], abc_res[[ind]]$adj.values)
        names(df) <- c("methods", "species", "pars", "posterior")
        df
    }
    # extract to data.frame
    all_abc <- do.call(rbind, lapply(neuralnet, extract_abc))
}

# all abc to df

abc_full <- do.call(rbind, lapply(c("5k", "50k", "500k"), posterior_to_df))
write.table(abc_full, file = "abc_full_3k.txt", row.names = FALSE)

