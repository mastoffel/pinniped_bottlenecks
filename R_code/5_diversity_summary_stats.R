# plots

#### load datasets
library(readxl)
# sheet numbers to load
dataset_names <- excel_sheets("data/processed/seal_data_largest_clust_and_pop_all_hw.xlsx")
load_dataset <- function(dataset_names) {
    read_excel("data/processed/seal_data_largest_clust_and_pop_all_hw.xlsx", sheet = dataset_names)
}
# load all datasets
all_seals <- lapply(dataset_names, load_dataset)
names(all_seals) <- dataset_names

# get 28 datasets, biggest clusters
ids <- c("antarctic_fur_seal",
    "galapagos_fur_seal",
    "stellers_sea_lion_cl_2",
    "grey_seal_orkneys",
    "harbour_seal_waddensee_cl_1",
    "galapagos_sea_lion",
    "south_american_fur_seal_cl_2",
    "hooded_seal" ,
    "mediterranean_monk_seal",
    "hawaiian_monk_seal",
    "bearded_seal_cl_1",
    "crabeater_seal",
    "leopard_seal" ,
    "arctic_ringed_seal",
    "ross_seal",
    "weddell_seal_cl_2",
    "northern_fur_seal_cl_1",
    "atlantic_walrus_cl_1",
    "nes_cl_4",
    "ses_cl_1",
    "california_sea_lion",
    "south_american_sea_lion",
    "new_zealand_sea_lion_cl_1",
    "saimaa_ringed_seal_cl_2",
    "lagoda_ringed_seal",
    "baltic_ringed_seal",
    "new_zealand_fur_seal" ,
    "australian_fur_seal")

# filter for 28 species datasets
all_seals <- all_seals[ids]
# rename all seals
names(all_seals) <- str_replace(names(all_seals), "_cl_[1-9]", "")


#### load bottleneck data
load("R/bottleneck_results.RData")
# rename all seals
rename_species <- function(bottle_tests){
    bottle_tests$id <- str_replace(bottle_tests$id, "_cl_[1-9]", "")
    bottle_tests$id <- str_replace(bottle_tests$id, "_HW", "")
    bottle_tests
}
bottleneck <- lapply(list(bottle_tests, bottle_tests_ratio), rename_species)
names(bottleneck) <- c("p_val", "ratio")

bottleneck_p_val <- bottleneck$p_val
bottleneck_ratio <- bottleneck$ratio


#### load additional data
# sheet numbers to load
harem_data <- read_excel("data/processed/overview.xlsx", sheet = 2)
harem_data <- harem_data[!(is.na(harem_data$species)), ]


# g2 and heterozygosity
library(inbreedR)

calc_g2s <- function(genotypes){
    g2_microsats(convert_raw(genotypes[, 4:ncol(genotypes)]), nboot = 1000, nperm = 1000)
}

g2s <- lapply(all_seals, calc_g2s)
# save(g2s, file = "all_g2s.RData")
load("R/all_g2s.RData")
options(scipen = 999)
g2_summary <- as.data.frame(do.call(rbind, lapply(g2s, function(x) c(x$g2, x$CI_boot, x$p_val))))
names(g2_summary ) <- c("g2", "CIlow", "CIup", "p_val")
g2_summary$species <- row.names(g2_summary)

library(dplyr)
g2_summary$species <- factor(g2_summary$species, levels = g2_summary$species[order(g2_summary$g2, 
    decreasing = TRUE)])
g2_summary <- g2_summary %>% mutate(p_val_ht = as.factor(as.numeric(p_val < 0.05)))

library(ggplot2)

ggplot(g2_summary, aes(x=species, y=g2, label = species, color = p_val_ht)) + 
    geom_errorbar(aes(ymin=CIlow, ymax=CIup), width=.1) +
    geom_line() +
    geom_point() +
    geom_text(angle = 90, vjust = 1.4, hjust = -0.5, size = 3) +
    theme(axis.text.x = element_blank())

library(inbreedR)
calc_hets <- function(genotypes){
    MLH(convert_raw(genotypes[, 4:ncol(genotypes)]))
}

all_hets <- lapply(all_seals, calc_hets)
all_hets_plus_name <- lapply(c(1:length(all_hets)), function(x) data.frame(het = all_hets[[x]], species = names(all_hets)[[x]]))
all_hets_df <- do.call(rbind, all_hets_plus_name)


all_hets_df$species <- factor(all_hets_df$species, levels = all_hets_df$species[order(all_hets_df$g2, 
    decreasing = TRUE)])

ggplot(all_hets_df, aes(x=reorder(as.factor(species), -het, FUN = median), y=het, label = as.factor(species))) + 
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


# make data.frame with all data
g2_summary
# save(g2_summary, file = "g2_summary")

all_seal_data <- g2_summary
all_seal_data <- all_seal_data[, c(5,1,2,3,4,6)]

# empirical heterozygosity
mean_mlh <- unlist(lapply(all_hets, mean, na.rm = TRUE))
var_mlh <- unlist(lapply(all_hets, var, na.rm = TRUE))

# expected heterozygosity
exp_het <- function(genotypes) {
    msats <- new("gtypes", gen.data = genotypes[, 4:ncol(genotypes)], ploidy = 2)
    # using modified version of LDgenepop
    He <- exptdHet(msats)
    out <- mean(He, na.rm = TRUE)
}

all_exp_het <- unlist(lapply(all_seals, exp_het))
all_exp_het - het

# add to data
all_seal_data$mean_mlh <- mean_mlh
all_seal_data$var_mlh <- var_mlh
all_seal_data$exp_mlh <- all_exp_het
all_seal_data$species <- as.character(all_seal_data$species)

# put together bottleneck results
bottleneck <- cbind(bottleneck_p_val, bottleneck_ratio)
bottleneck[[6]] <- NULL
bottleneck$id

# match names
bottleneck$id[[10]] <- "arctic_ringed_seal"
reorder <- unlist(lapply(all_seal_data$species, function(x) which(str_detect(bottleneck$id, x))))
bottleneck <- bottleneck[reorder, ]

# put everything together
all_seal_info <- cbind(all_seal_data, bottleneck)
all_seal_info$id <- NULL

# load harem data
library(readxl)
# sheet numbers to load
dataset_names <- excel_sheets("data/overview.xlsx")
# load all datasets
harem_data <- read_excel("data/overview.xlsx", sheet = 2)


# put everything together
seals <- cbind(harem_data[1:28, ], all_seal_info)
seals[[9]] <- NULL

# get sample size and locus number
names(all_seals)
seals$sample_size <- unlist(lapply(all_seals, nrow))
seals$loci_number <- unlist(lapply(all_seals, function(x) (ncol(x) - 3) / 2 ))

save(seals, file = "seal_summary_data.RData")










# allelic richness
library(strataG)

allelic_richness <- function(genotypes) {
    msats <- new("gtypes", gen.data = genotypes[, 4:ncol(genotypes)], ploidy = 2)
    # using modified version of LDgenepop
    ar <- allelicRichness(msats)
    out <- mean(ar, na.rm = TRUE)
}

AR <- unlist(lapply(all_seals, allelic_richness))
names(AR) <- names(all_seals)
AR
lagoda <- all_seals$lagoda_ringed_seal
msats_lagoda <- new("gtypes", gen.data = lagoda[, 4:ncol(lagoda)], ploidy = 2)
allelicRichness(msats_lagoda)



