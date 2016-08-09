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
all_seals <- read_excel_sheets("data/processed/seal_data_largest_clust_and_pop.xlsx")

# get 28 datasets, biggest clusters
ids <- c("antarctic_fur_seal", "galapagos_fur_seal", "stellers_sea_lion_cl_2",
    "grey_seal_orkneys", "harbour_seal_waddensee_cl_2", "galapagos_sea_lion",
    "south_american_fur_seal_cl_2", "hooded_seal", "mediterranean_monk_seal",
    "hawaiian_monk_seal", "bearded_seal_cl_2", "crabeater_seal",
    "leopard_seal", "arctic_ringed_seal", "ross_seal",
    "weddell_seal_cl_1", "northern_fur_seal_cl_1", "atlantic_walrus_cl_1",
    "nes_cl_1", "ses_cl_1", "california_sea_lion", "south_american_sea_lion",
    "new_zealand_sea_lion", "saimaa_ringed_seal_cl_2", "lagoda_ringed_seal",
    "baltic_ringed_seal", "new_zealand_fur_seal", "australian_fur_seal")

# filter for 28 species datasets
all_seals <- all_seals[ids]
# rename all seals
names(all_seals) <- str_replace(names(all_seals), "_cl_[1-9]", "")


#### load bottleneck data
bottleneck <- read_excel("data/processed/bottleneck_results.xlsx")

# rename all seals
rename_species <- function(bottle_tests){
    bottle_tests$id <- str_replace(bottle_tests$id, "_cl_[1-9]", "")
    bottle_tests$id <- str_replace(bottle_tests$id, "_HW", "")
    bottle_tests
}
bottleneck <- rename_species(bottleneck)


#### load additional data
# sheet numbers to load
harem_data <- read_excel("data/processed/overview.xlsx", sheet = 2)
harem_data <- harem_data[!(is.na(harem_data$species)), ]


### diversity
# g2 
library(inbreedR)
calc_g2s <- function(genotypes){
    g2_microsats(convert_raw(genotypes[, 4:ncol(genotypes)]), nboot = 1000, nperm = 1000)
}
# g2s <- lapply(all_seals, calc_g2s)
# save(g2s, file = "all_g2s.RData")
load("all_g2s.RData")

# put it all in one data frame
options(scipen = 999)
g2_summary <- as.data.frame(do.call(rbind, lapply(g2s, function(x) c(x$g2, x$CI_boot, x$p_val))))
names(g2_summary ) <- c("g2", "CIlow", "CIup", "p_val")
g2_summary$species <- row.names(g2_summary)


# heterozygosity
calc_hets <- function(genotypes){
    MLH(convert_raw(genotypes[, 4:ncol(genotypes)]))
}

all_hets <- lapply(all_seals, calc_hets)
all_hets_plus_name <- lapply(c(1:length(all_hets)), function(x) data.frame(het = all_hets[[x]], species = names(all_hets)[[x]]))
all_hets_df <- do.call(rbind, all_hets_plus_name)


### relationship between heteorzygosity and number of loci or number of individuals?
sum_hets <- all_hets_df %>% 
                group_by(species) %>%
                summarise(mean_het = mean(het, na.rm = TRUE),
                          sd_het = sd(het, na.rm = TRUE)) 
# get loci and individuals
desc_seals <- do.call(rbind, lapply(all_seals, function(x) out <- data.frame(nloc = (ncol(x)-3)/2, nind = nrow(x))))

# put het, nloc and nind in one data.frame
sumstats_diversity <- cbind(sum_hets , desc_seals)

# plotting
# ggplot(data = sumstats_diversity, aes(x = mean_het, y = nloc)) + geom_point()
# # model including both nloc and nind
# mod <- glm(data = sumstats_diversity, mean_het ~ nloc + nind)
# summary(mod)
# plot(mod)


diversity_stats <- cbind(g2_summary, sumstats_diversity)
# save(g2_summary, file = "g2_summary")

# add mssumstats to data
seal_stats <- do.call(rbind, lapply(all_seals, function(x) mssumstats(x[4:ncol(x)], type = "microsats")))
diversity_stats <- cbind(diversity_stats, seal_stats)

# put together bottleneck results
bottleneck$id

# match names
reorder <- unlist(lapply(as.character(diversity_stats$species), function(x) which(str_detect(bottleneck$id, x))))
bottleneck <- bottleneck[reorder, ]

# put everything together
all_seal_data <- cbind(diversity_stats, bottleneck)
all_seal_data$id <- NULL
names(all_seal_data)

# load harem data

# sheet numbers to load
dataset_names <- excel_sheets("data/processed/overview.xlsx")
# load all datasets
harem_data <- read_excel("data/processed/overview.xlsx", sheet = 2)

# put everything together
seals <- cbind(harem_data[1:28, ], all_seal_data)
# delete doubled species column
seals <- seals[-which(duplicated(names(seals)))]
# delete rownames
rownames(seals) <- 1:nrow(seals)

# check how many NAs
lapply(all_seals, function(x) rowSums(is.na(x)))

# calculate allelic richness
?allele2locus

seals_new <- lapply(all_seals, function(x) {
    names(x) <- str_replace_all(names(x), ".", "_")
    x
})

calc_allelic_richness <- function(genotypes, min.alleles = 80) {
   # join_alleles <- sealABC::allele2locus(genotypes[4:ncol(genotypes)], "/")
    g_types_geno <- strataG::df2gtypes(genotypes, ploidy = 2, id.col = NULL, strata.col = NULL,
        loc.col = 4)
    genind_geno <- gtypes2genind(g_types_geno)
    all_rich <- allel.rich(genind_geno, min.alleles)
    out <- as.numeric(all_rich$mean.richness)
}
# idea: maybe always subsample lowest number of loci (5?) , calculate allelic richness and then take average
all_richness <- unlist(lapply(seals_new, calc_allelic_richness, 60))

seals$allelic_richness <- all_richness

write.xlsx(seals, file = "data/processed/all_data_seals.xlsx", row.names = FALSE)

# save(seals, file = "seal_summary_data.RData")










