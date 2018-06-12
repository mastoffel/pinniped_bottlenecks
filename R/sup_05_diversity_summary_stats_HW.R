# Calculates rarified summary statistics

# files needed:
# (1) seal_data_largest_clust_and_pop_29.xlsx
# (2) bottleneck_results_29.txt or bottleneck_results_29_cl.txt
# (3) overview_29.xlsx

# plots
library(stringr)
library(sealABC)
library(readxl)
library(inbreedR)
library(ggplot2)
library(dplyr)
library(adegenet)
library(strataG)
# library(hierfstat)
library(strataG)
library(readr)
library(stringr)


shortcut_save <- "HW"

all_seals <- read_excel_sheets("data/processed/seal_data_largest_clust_and_pop_all_hw_30.xlsx")

#### load bottleneck data
bottleneck <- read_delim("data/processed/bottleneck_results_HW_30.txt", col_names = TRUE, delim = " ")

#### load additional data
# sheet numbers to load // just load second sheet
harem_data <- read_excel("data/processed/overview_30.xlsx", sheet = 2)
harem_data <- harem_data[!(is.na(harem_data$species)), ]

### diversity
# g2 
# check if calculated
g2_file <- paste0("data/processed/g2_summary_", shortcut_save, "_data.txt")
if(!file.exists(g2_file)){
    library(inbreedR)
    calc_g2s <- function(genotypes){
        g2_microsats(convert_raw(genotypes[, 4:ncol(genotypes)]), nboot = 10, nperm = 10) # increase at some point
    }
    g2s <- lapply(all_seals, calc_g2s)
    # save(g2s, file = "data/processed/full_data_all_g2s_29.RData")
    # load("data/processed/full_data_all_g2s_29.RData")
    g2s_full <- g2s
    
    # put it all in one data frame
    options(scipen = 999)
    g2_summary <- as.data.frame(do.call(rbind, lapply(g2s_full, function(x) c(x$g2, x$CI_boot, x$p_val))))
    names(g2_summary ) <- c("g2", "CIlow", "CIup", "p_val")
    g2_summary$species <- row.names(g2_summary)
    # reorder columns
    g2_summary <- g2_summary[c("species", "g2", "CIlow", "CIup", "p_val")]
    # write to txt
    readr::write_delim(g2_summary, path = g2_file, col_names = TRUE)
}
# load g2s
g2s <- read_delim(g2_file, delim = " ")

bottleneck$id

# calculate other summary statistics based on resampling 10 individuals
?mssumstats
sumstats_file <- paste0("data/processed/sumstats_30_", shortcut_save, ".txt")
if(!file.exists(sumstats_file)){
    set.seed(1122)
    sumstats <- do.call(rbind, lapply(all_seals, function(x) mssumstats(x, start_geno = 4, mratio = "loose", 
        rarefaction = TRUE, nresamp = 1000, nind = 10)))
    
    readr::write_delim(sumstats, path = sumstats_file, col_names = TRUE)
}

sumstats <- read_delim(sumstats_file, delim = " ")


# bind
diversity_stats <- cbind(g2s, sumstats)
# save(g2_summary, file = "g2_summary")

# plot diversity
ggplot(diversity_stats, aes(x = mratio_mean, y = species)) +
    geom_point() +
    geom_errorbarh(aes(xmax = mratio_mean.CIhigh, xmin = mratio_mean.CIlow))

# rename id variable in bottleneck table
names(bottleneck)[1] <- "species"
# put together bottleneck results
# check that all names are equal
sum(bottleneck$species %in% diversity_stats$species)
# match names

# join them together
###### double check here when running, subsetting for non cl, non pop species###
all_seal_data <- left_join(diversity_stats[1:30, ], bottleneck, by = "species")
names(all_seal_data)

# remove cl names from all_seal_data to join with harem data
all_seal_data$species <- str_replace(all_seal_data$species, "_cl_[1-9]", "")

# # load all datasets
# harem_data <- read_excel("data/processed/overview.xlsx", sheet = 2)
seals <- left_join(harem_data, all_seal_data, by = "species")
# recode IUCN
names(seals)[7] <- "IUCN_rating"
# check how many NAs
lapply(all_seals, function(x) rowSums(is.na(x)))

# get loci and individuals
desc_seals <- do.call(rbind, lapply(all_seals, function(x) out <- data.frame(nloc = (ncol(x)-3)/2, nind = nrow(x))))
desc_seals$species <- rownames(desc_seals)
seals <- left_join(seals, desc_seals, by = "species")

# write all stats to file
sumstats_data <- paste0("data/processed/all_data_seals_rarefac10_30_", shortcut_save, ".csv")
write_excel_csv(seals, sumstats_data)


