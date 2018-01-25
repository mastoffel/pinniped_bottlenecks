library(readxl)
library(devtools)
# install_github("hadley/readxl")
library(readxl)
# sheet numbers to load
dataset_names <- excel_sheets("data/all_seals_full_plus_pop_cluster.xls")

load_dataset <- function(dataset_names) {
    read_excel("data/all_seals_full_plus_pop_cluster.xls", sheet = dataset_names)
}

# load all datasets
all_seals <- lapply(dataset_names, load_dataset)
names(all_seals) <- dataset_names

# sample_size, number_of_loci and total_genotypes
data_info <- function(x){
    sample_size <- nrow(x)
    number_of_loci <- (ncol(x)-2)/2
    total_genotypes <- sample_size * number_of_loci
    populations <- length(names(table(x[2])))
    out <- data.frame(sample_size, number_of_loci, total_genotypes, populations)
}

all_info <- do.call(rbind, lapply(all_seals, data_info))
write.table(all_info, "info.txt", row.names = FALSE, sep = "\t")

test <- all_seals[[4]]
test$pop
table(test$pop)
lapply(all_seals, function(x) table(x$pop))
lapply(all_seals, function(x) names(x)[1])
