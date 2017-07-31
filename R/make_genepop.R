# reads in an excel file with multiple sheets of microsatellite data and makes genepop files for each

library(devtools)
# devtools::install_github("hadley/readxl")
library(readxl)

# sheet numbers to load
dataset_names <- excel_sheets("../data/processed/seal_data_largest_clust_and_pop.xlsx")

dataset_names_hw <- excel_sheets("../data/processed/seal_data_largest_clust_and_pop_all_hw.xlsx")

load_dataset <- function(dataset_names) {
        read_excel("../data/processed/seal_data_largest_clust_and_pop.xlsx", sheet = dataset_names)
}

load_dataset_hw <- function(dataset_names) {
        read_excel("../data/processed/seal_data_largest_clust_and_pop_all_hw.xlsx", sheet = dataset_names)
}

# load all datasets
all_seals <- lapply(dataset_names, load_dataset)
names(all_seals) <- dataset_names

all_seals_hw <- lapply(dataset_names_hw, load_dataset_hw)
names(all_seals_hw) <- dataset_names_hw

# subset for biggest clusters

all_seals <- all_seals[1:39]
dataset_names <- dataset_names[1:39]

all_seals_hw <- all_seals_hw[1:39]
dataset_names_hw <- dataset_names_hw[1:39]

# MAKE GENEPOP FILE

make_genepop <- function(x) {
x <- x[-c(1,3)]
colnames(x)[1] <- "POP"
x[1] <- "1,"

# get locus names
loci <- matrix(c(colnames(x)[-1],colnames(x[1])))
loci <- unique(gsub("_a|_b", "", loci))

# convery genotypes into 3 digits
library(stringr)

threedigits <- function(x) {
        x[is.na(x)] <- "000"
        xchar <- as.character(x)
        str_length(xchar)
        new_x <- str_pad(xchar, 3, pad = "0")
}

genotypes <- apply(x[-1], 2, threedigits)

# one column per locus
short_geno <- matrix(nrow = nrow(genotypes), ncol = ncol(genotypes) / 2)
short_geno <- data.frame(short_geno, stringsAsFactors = FALSE)
length_data <- ncol(genotypes)
genotypes <- data.frame(genotypes, stringsAsFactors = FALSE)
col_num <- 1
for (i in seq(from = 1, to = length_data, by = 2)) {
        short_geno[col_num] <- paste0(genotypes[[i]], genotypes[[i+1]])
        col_num <- col_num + 1
}

# join everything together
genotypes <- cbind(x[1], short_geno)
pop <- matrix(ncol = length(genotypes),
                  nrow = nrow(loci) + 1)

pop[1,1] <- "Title Line: genotypes"
pop[2:nrow(pop),1] <- loci
pop <- as.data.frame(pop)
colnames(genotypes) <- colnames(pop)
file <- rbind(pop, genotypes)

}

out <- lapply(all_seals, make_genepop)
names(out) <- dataset_names

out_hw <- lapply(all_seals_hw, make_genepop)
names(out_hw) <- dataset_names_hw

# write files

lapply(1:length(out), function(i) write.table(out[[i]], 
                                               file = paste("../data/processed/genepop_files/", names(out)[i], "_genepop.txt", sep = ""),
                                               na = "", row.names = F, col.names = F, quote = F))

lapply(1:length(out_hw), function(i) write.table(out_hw[[i]], 
                                             file = paste("../data/processed/genepop_files/HW/", names(out_hw)[i], "_HW_genepop.txt", sep = ""),
                                             na = "", row.names = F, col.names = F, quote = F))
