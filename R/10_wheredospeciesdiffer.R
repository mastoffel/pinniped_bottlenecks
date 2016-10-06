library(sealABC)
library(reshape2)
library(ggrepel)
library(ggplot2)
library(readxl)
library(ggthemr)
library(sealABC)
library(dplyr)

ggthemr('fresh', spacing = 1, text_size = 15)

# ggthemr_reset()
all_seals <- read_excel_sheets("data/processed/seal_data_largest_clust_and_pop.xlsx")
seals <- read_excel("data/processed/seal_data_complete.xlsx")
names(seals)[11] <- "IUCN_rating"
seal_names <- seals$seal_names_real




all_sumstats <- do.call(rbind, lapply(all_seals, function(x) mssumstats(x[4:ncol(x)], type = "microsats", data_type = "empirical")))
names_all_sumstats <- row.names(all_sumstats)

all_sumstats %<>% 
    slice(1:28) %>%
    mutate(species = names_all_sumstats[1:28])

ggplot(all_sumstats, aes(num_alleles_mean, het_excess)) +
    geom_point(size = 3, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))



arrow <- arrow(length = unit(0.4, "cm"), type = "closed")

ggplot(seals, aes(num_alleles_mean, mratio_mean_new)) +
    geom_point(size = 3, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))
    geom_label_repel(aes(label = seal_names), size = 2) 
   



# try calculating proportion of low alleles with frequencies under 0.05
calc_allelefreq <- function(geno){
    genotypes <- geno[4:ncol(geno)]
    g_types_geno <- strataG::df2gtypes(genotypes, ploidy = 2, id.col = NULL, strata.col = NULL,
        loc.col = 1)
    out <- strataG::alleleFreqs(g_types_geno)
}

prop_low_af <- function(afs){
    low_afs <- (afs[, "freq"] / sum(afs[, "freq"])) < 0.05
    prop_low <- sum(low_afs) / length(low_afs)
}

mean_low_af <- function(all_afs){
    out <- mean(unlist(lapply(all_afs, prop_low_af)))
}


all_mean_afs <- lapply(all_af, mean_low_af)

all_af <- lapply(all_seals, calc_allelefreq)

seal_afs <- data.frame(species = names(all_mean_afs), afs = as.numeric(all_mean_afs))

ggplot(seal_afs[1:28, ], aes(species, afs)) +
    geom_point(size = 3, alpha = 0.5) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))




# calculate mRatio with strataG

calc_mratio <- function(geno){
    genotypes <- geno[4:ncol(geno)]
        g_types_geno <- strataG::df2gtypes(genotypes, ploidy = 2, id.col = NULL, strata.col = NULL,
            loc.col = 1)
    out <- strataG::mRatio(g_types_geno)
}

all_mratios <- lapply(all_seals, calc_mratio)
all_mratios <- lapply(all_mratios, mean, na.rm = TRUE)

seals$mratio_mean_new <- all_mratios[1:28]


nes <- all_seals$nes
nes_pop <- all_seals$nes_pop
nes_clust <- all_seals$nes_cl_2

plot(names(all_mratios )[1:28], as.numeric(all_mratios ))

nes_ss <- mssumstats(nes[4:ncol(nes)], type = "microsats")
nespopss <- mssumstats(nes_pop[4:ncol(nes_pop)], type = "microsats")
nesclust_ss <- mssumstats(nes_clust [4:ncol(nes_clust )], type = "microsats")

nes_all <- rbind(nes_ss, nespopss, nesclust_ss)

lapply(nes[4:ncol(nes)], range, na.rm = TRUE)



# some checks
all_seals <- all_seals[1:28]

seals_gtypes <- lapply(all_seals, function(x) strataG::df2gtypes(x[4:ncol(x)], ploidy = 2, id.col = NULL, strata.col = NULL, loc.col = 1))

# calculate proportion fo low frequency alleles for one locus
prop_low_af <- function(g_types_geno){
    low_afs <- (afs[, "freq"] / sum(afs[, "freq"])) < 0.05
    prop_low <- sum(low_afs) / length(low_afs)
}

all_afs <- lapply(seals_gtypes, function(x) {
    afs <- strataG::alleleFreqs(x)
    prop_low_af <- function(g_types_geno){
        low_afs <- (afs[, "freq"] / sum(afs[, "freq"])) < 0.05
        prop_low <- sum(low_afs) / length(low_afs)
    }
    
})
