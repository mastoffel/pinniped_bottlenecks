# script to check all microsatellite datasets

library(readxl)
# sheet numbers to load
dataset_names <- excel_sheets("data/processed/seal_data_largest_clust_and_pop.xlsx")

load_dataset <- function(dataset_names) {
    read_excel("data/processed/seal_data_largest_clust_and_pop.xlsx", sheet = dataset_names)
}
# load all datasets
all_seals <- lapply(dataset_names, load_dataset)
names(all_seals) <- dataset_names


all_seals[[2]]
# check coding of msats
lapply(all_seals, function(x) apply(x[, 4:ncol(x)], 2, range, na.rm = T))


    
    
    
    

lapply(all_seals)