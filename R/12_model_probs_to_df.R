# plotting posteriors from abc

# first plot heatmap of model probabilities
pop_sizes <- c("5k", "50k", "500k")
# load all posteriors
abc_full <- read.table("abc_full_3k.txt", header = TRUE)

# load 
abc_model_select <- function(pop_size){
    # load data
    to_load <- paste0("data/processed/abc_estimates/", pop_size, "_model_selection.txt")
    model_select <- read.table(to_load, header = TRUE)
}
model_probs <- do.call(rbind, lapply(pop_sizes, abc_model_select))

library(xlsx)
write.xlsx(model_probs, "data/processed/model_probs_300k.xlsx", row.names = TRUE)

