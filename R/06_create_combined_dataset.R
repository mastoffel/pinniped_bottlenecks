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
phylo <- read_excel("data/raw/phylogeny_overview_29.xlsx", sheet = 2, col_names = F)[1:29, ]
names(phylo) <- c("latin", "common", "species")

# All genetic and some life history variables
seals <- read_csv("data/processed/all_data_seals_rarefac10_29.csv")

# join datasets
seals2 <- left_join(phylo, seals, by = "species")

# rest of the life history data
krueger_data <- data.frame(read_excel("data/processed/seal_data_krueger.xlsx", sheet = 1, col_names = T))%>% 
                dplyr::select(dataset_name, birth_mass:Age_primiparity) %>% 
                rename(species = dataset_name)

# join data to full dataset by dataset names
seals3 <- dplyr::left_join(seals2, krueger_data, by = "species")

# rearrange to have life-history data first
seals_rearranged <- seals3[c(1:11, 76:86, 12:75)]

seals_rearranged$common

# produce short names for plotting
short <- c("W", "NFS", "GFS", "SAFS", "AntFS", "SAntFS", "NZFS", "AFS", "GSL", "CSL", "NZSL", "SASL", "SSL", "HMS", "MMS", "NES", "SES",
            "CS", "WS", "LS", "RoS", "BS", "HoodS", "HS", "RS", "SRS", "LRS", "BRS", "GS")

# add abreviations
seals_rearranged$short <- short

# write to excel
write_excel_csv(seals_rearranged, "data/processed/seal_data_complete_rarefac10_29.csv")
# 



