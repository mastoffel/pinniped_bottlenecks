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
phylo <- read_excel("data/raw/phylogeny_overview.xlsx", sheet = 2, col_names = F)[1:28, ]
names(phylo) <- c("latin", "common", "species")
# seals <- read_csv("data/processed/data_noLH_rarefac10.csv")
seals <- read_csv("data/processed/all_data_seals_rarefac10.csv")

# join datasets
seals2 <- left_join(phylo, seals, by = "species")

# get krueger data and order
krueger_data <- data.frame(read_excel("data/processed/seal_data_krueger.xlsx", sheet = 1, col_names = T))%>% 
                 dplyr::select(dataset_name, birth_mass:Age_primiparity) %>% 
                rename(species = dataset_name)

# join data to full dataset
seals3 <- dplyr::left_join(seals2, krueger_data, by = "species")
# rearrange to have life-history data first
seals_rearranged <- seals3[c(1:11, 76:86, 12:75)]
seals_rearranged$common[seals_rearranged$species == "new_zealand_sea_lion"] <- "New Zealand Sea Lion"
names(seals_rearranged)[9] <- "IUCN_rating"

write_excel_csv(seals_rearranged, "data/processed/seal_data_complete_rarefac10.csv")
# 





# seals_final <- seals[c("species" ,"seal_names_real", "latin" ,"family", "harem_size", "mean_SSD", "male_weight", "female_weight",
# "male_length", "female_length", "IUCN rating", "Abundance", "BreedingType", "g2", "CIlow",  "CIup", "p_val", "allelic_richness",
# "mean_het", "sd_het", "num_alleles_mean",  "num_alleles_sd",  "mean_allele_size_var", "sd_allele_size_var",
# "mean_allele_range", "sd_allele_range",  "exp_het_mean" , "exp_het_sd", "obs_het_mean","obs_het_sd",
# "mratio_mean", "mratio_sd" , "prop_low_afs_mean", "prop_low_afs_sd", "het_excess", "IAM_Wilc_Exc",  "TPM70_Wilc_Exc" , "TPM90_Wilc_Exc" ,   "TPM95_Wilc_Exc",
# "SMM_Wilc_Exc" , "IAM_ratio",  "TPM70_ratio", "TPM90_ratio",   "TPM95_ratio", "SMM_ratio",  "birth_mass",
# "breeding_season_length" , "lactation_length", "life_span_years" ,  "latitude" ,    "rel_birth_weight",
# "Longevity", "Generation_time", "Age_primiparity", "nloc", "nind")]
