# plots
library(stringr)
library(sealABC)
library(readxl)
library(inbreedR)
library(ggplot2)
library(dplyr)
library(xlsx)
library(adegenet)
library(PopGenReport)
library(strataG)
library(hierfstat)
library(strataG)
library(readr)
all_seals <- read_excel_sheets("data/processed/seal_data_largest_clust_and_pop_29.xlsx")

# set to TRUE is largest clusters instead of full datasets
calc_on_cluster <- FALSE

if (calc_on_cluster){
    
    # get 28 datasets, biggest clusters
    ids <- c("antarctic_fur_seal", "galapagos_fur_seal", "stellers_sea_lion_cl_1",
        "grey_seal_orkneys", "harbour_seal_waddensee_cl_1", "galapagos_sea_lion",
        "south_american_fur_seal_cl_1", "hooded_seal", "mediterranean_monk_seal",
        "hawaiian_monk_seal", "bearded_seal_cl_1", "crabeater_seal",
        "leopard_seal", "arctic_ringed_seal", "ross_seal",
        "weddell_seal_cl_2", "northern_fur_seal_cl_1", "atlantic_walrus_cl_2",
        "nes_cl_2", "ses_cl_4", "california_sea_lion", "south_american_sea_lion",
        "new_zealand_sea_lion", "saimaa_ringed_seal_cl_1", "lagoda_ringed_seal",
        "baltic_ringed_seal", "new_zealand_fur_seal", "australian_fur_seal")
    
    # filter for 28 species datasets
    all_seals <- all_seals[ids]
    # rename all seals
    names(all_seals) <- str_replace(names(all_seals), "_cl_[1-9]", "")
    
} else {
    # take all full datasets
    ids <- names(all_seals)[1:29] ### number of species now 29
    ids %in% names(all_seals)
    # filter for 28 species datasets
    all_seals <- all_seals[ids]
}


#### load bottleneck data
bottleneck <- read_delim("data/processed/bottleneck_results_29.txt", col_names = TRUE, delim = " ")

# rename all seals
# rename_species <- function(bottle_tests){
#     bottle_tests$id <- str_replace(bottle_tests$id, "_cl_[1-9]", "")
#     bottle_tests$id <- str_replace(bottle_tests$id, "_HW", "")
#     bottle_tests
# }
# bottleneck <- rename_species(bottleneck)


#### load additional data
# sheet numbers to load // just load second sheet
harem_data <- read_excel("data/processed/overview_29.xlsx", sheet = 2)
harem_data <- harem_data[!(is.na(harem_data$species)), ]


### diversity
# g2 
# check if calculated
if(!file.exists("data/processed/g2_summary_full_data.txt")){
    library(inbreedR)
    calc_g2s <- function(genotypes){
        g2_microsats(convert_raw(genotypes[, 4:ncol(genotypes)]), nboot = 1000, nperm = 1000)
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
    readr::write_delim(g2_summary, path = "data/processed/g2_summary_full_data_29.txt", col_names = TRUE)
}
# load g2s
g2s <- read_delim("data/processed/g2_summary_full_data_29.txt", delim = " ")

bottleneck$id

# calculate other summary statistics based on resampling 10 individuals
?mssumstats

if(!file.exists("data/processed/sumstats_29.txt")){
    
sumstats <- do.call(rbind, lapply(all_seals, function(x) mssumstats(x, start_geno = 4, mratio = "loose", 
            rarefaction = TRUE, nresamp = 1000, nind = 10)))

readr::write_delim(sumstats, path = "data/processed/sumstats_29.txt", col_names = TRUE)
}
sumstats <- read_delim("data/processed/sumstats_29.txt", delim = " ")


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
all_seal_data <- left_join(diversity_stats, bottleneck, by = "species")
names(all_seal_data)

# load harem data
# sheet numbers to load
# dataset_names <- excel_sheets("data/processed/overview.xlsx")
# # load all datasets
# harem_data <- read_excel("data/processed/overview.xlsx", sheet = 2)
seals <- left_join(harem_data, all_seal_data, by = "species")

# check how many NAs
lapply(all_seals, function(x) rowSums(is.na(x)))

# get loci and individuals
desc_seals <- do.call(rbind, lapply(all_seals, function(x) out <- data.frame(nloc = (ncol(x)-3)/2, nind = nrow(x))))
seals <- cbind(seals, desc_seals)

if(!file.exists("data/processed/all_data_seals_rarefac10_29.csv")){
write_excel_csv(seals, "data/processed/all_data_seals_rarefac10_29.csv")
}





