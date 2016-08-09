# create combined dataset with Ollies data in the sequence of Phylogeny

# plots for the talk
# install_github("mastoffel/sealABC", dependecies = TRUE)
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(stringr)

library(hierfstat)
?hierfstat
# load phyologenetic sequence

# load all datasets
phylo <- read_excel("data/raw/phylogeny_overview.xlsx", sheet = 2, col_names = F)
seals <- read_excel("data/processed/all_data_seals.xlsx")

# get phylo order
reorder <- unlist(lapply(phylo$X3, function(x) which(str_detect(seals$species, x))))
seals <- seals[reorder, ]
seals[, 2:6] <- lapply(seals[, 2:6], as.numeric)

# get krueger data and order
krueger_data <- read_excel("data/processed/seal_data_krueger.xlsx", sheet = 1, col_names = T)

# decide at some point whether to include ringed seal subspecies
krueger_data <- krueger_data[1:38, ]
krueger_data <- krueger_data[!is.na(krueger_data$dataset_name), ]

# prepare krueger data for matching
seals$species %in% krueger_data$dataset_name 
# reorder krueger data
reorder <- unlist(lapply(seals$species, function(x) which(str_detect(krueger_data$dataset_name, x))))
krueger_data <- krueger_data[reorder, ]
# bind seals and krueger_data
seals <- cbind(seals, krueger_data)
names(seals)
seals <- seals[-c(2:4)] # duplicated data, rather take from krueger data frame

#seals <- seals[-which(duplicated(names(seals)))]

# rename the hookers sea lion
seal_names_real <- phylo$X2[1:28]
seal_names_real[10] <- "New Zealand Sea Lion"

names(seals)
seals$seal_names_real <- seal_names_real

names(seals)

seals_final <- seals[c("species" ,"seal_names_real", "latin" ,"family", "harem_size", "mean_SSD", "male_weight", "female_weight",
    "male_length", "female_length", "IUCN rating", "Abundance",  "g2", "CIlow",  "CIup", "p_val", "allelic_richness",
    "mean_het", "sd_het", "num_alleles_mean",  "num_alleles_sd",  "mean_allele_size_sd", "sd_allele_size_sd",
    "mean_allele_range", "sd_allele_range",  "exp_het_mean" , "exp_het_sd", "obs_het_mean","obs_het_sd",
    "mratio_mean", "mratio_sd" , "IAM_Wilc_Exc",   "TPM70_Wilc_Exc" , "TPM90_Wilc_Exc" ,   "TPM95_Wilc_Exc",
    "SMM_Wilc_Exc" , "IAM_ratio",  "TPM70_ratio", "TPM90_ratio",   "TPM95_ratio", "SMM_ratio",  "birth_mass",
    "breeding_season_length" , "lactation_length", "life_span_years" ,    "latitude" ,    "rel_birth_weight",
    "Longevity", "Generation_time", "Age_primiparity", "nloc", "nind")]

library(xlsx)
write.xlsx(seals_final, "data/processed/seal_data_complete.xlsx", row.names = FALSE)
# 