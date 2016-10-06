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

# weird: 
# bearded_seal$L_LS8_b, bearded_seal$L_LS20_b, 
# arctic_ringed_seal$L_LS8_a arctic_ringed_seal$L_LS8_b
# ross_seal$L_LS8_b and a, ross_seal$L_LS20_a and b, 
# weddell_seal$L_LS8_a and b, L_LS20_a and b


#
# genotypes <- all_seals[[1]][4:ncol(genotypes)]
# check_repeat_size <- function(genotypes){
#         lapply(seq(from = 1, to = ncol(genotypes), by = 2), function(x) table(genotypes[c(x:(x+1))]))    
#     
#     
#     
    
    
    
    

lapply(all_seals)