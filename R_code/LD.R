# linkage disequilibrium
library(readxl)

# sheet numbers to load
dataset_names <- excel_sheets("data/seal_data_largest_clust_and_pop_all_hw.xlsx")

load_dataset <- function(dataset_names) {
    read_excel("data/seal_data_largest_clust_and_pop_all_hw.xlsx", sheet = dataset_names)
}
# load all datasets
all_seals <- lapply(dataset_names, load_dataset)
names(all_seals) <- dataset_names

source("LDgenepop_mod.R")
# try with strataG
library(strataG)

test_LD <- function(genotypes) {
    msats <- new("gtypes", gen.data = genotypes[, 4:ncol(genotypes)], ploidy = 2)
    # using modified version of LDgenepop
    ld <- LDgenepop_mod(msats, exec = "/Users/martin/programs/Genepop")
    out <- (p.adjust(ld$p.value, method = "bonferroni") < 0.01)
}

all_LD <- lapply(all_seals, test_LD)
loci_in_ld_per_ind <- lapply(all_LD, function(x) sum(x) / length(x))
max(unlist(loci_in_ld_per_ind))
