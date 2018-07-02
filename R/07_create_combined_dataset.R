# create combined dataset with ecology/life-history/demography data in the sequence of Phylogeny

# files needed:
# (1) phylogeny_overview_30.xlsx        # phylogeny
# (2) all_data_seals_rarefac10_30.csv   # summary statistics
# (3) seal_data_krueger.xlsx            # life history

library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ggthemes)
library(stringr)
library(hierfstat)
?hierfstat

# on cluster or full data?
calc_on_cluster <- FALSE
calc_on_hw <- TRUE

if (calc_on_cluster) {
    save_files <- "_cl"    
} else if (calc_on_hw) {
    save_files <- "_HW"   
} else {
    save_files <- ""
}

# load all datasets
phylo <- read_excel("data/raw/phylogeny_overview_30.xlsx", sheet = 2, col_names = F)[1:30, ]
names(phylo) <- c("latin", "common", "species")

# All genetic and some life history variables
seals <- read_csv(paste0("data/processed/all_data_seals_rarefac10_30", save_files, ".csv"))

# join datasets
seals2 <- left_join(phylo, seals, by = "species")

# rest of the life history data from Krueger et al. 2014
krueger_data <- data.frame(read_excel("data/processed/seal_data_krueger.xlsx", sheet = 1, col_names = T))%>% 
                dplyr::select(dataset_name, birth_mass:Age_primiparity) %>% 
                dplyr::rename(species = dataset_name)

# join data to full dataset by dataset names
seals3 <- dplyr::left_join(seals2, krueger_data, by = "species")

# rearrange to have life-history data first
seals_rearranged <- seals3[c(1:11, 76:86, 12:75)]

seals_rearranged$common

# produce short names for plotting
short <- c("W", "NFS", "GFS", "SAFS", "AntFS", "SAntFS", "GuaFS", "NZFS", "AFS", "GSL", "CSL", "NZSL", "SASL", "SSL", "HMS", "MMS", "NES", "SES",
            "CS", "WS", "LS", "RoS", "BS", "HoodS", "HS", "RS", "SRS", "LRS", "BRS", "GS")

# add abreviations
seals_rearranged$short <- short

# write to excel
write_excel_csv(seals_rearranged, paste0("data/processed/seal_data_complete_rarefac10_30", save_files, ".csv"))
# 



