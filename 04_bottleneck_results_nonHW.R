# process emilys bottleneck output files 
# full data (incl non-HW loci)

library(data.table)  # faster fread() and better weekdays()
library(dplyr)       # consistent data.frame operations
library(purrr)       # consistent & safe list/vector munging
library(tidyr)       # consistent data.frame cleaning
library(lubridate)   # date manipulation
library(ggplot2)     # base plots are for Coursera professors
library(scales)      # pairs nicely with ggplot2 for plot label formatting
library(gridExtra)   # a helper for arranging individual ggplot objects
library(ggthemes)    # has a clean theme for ggplot2
library(viridis)     # best. color. palette. evar.
library(knitr)       # kable : prettier data.frame output
library(stringr)
library(reshape2)
library(xlsx)
library(ggthemes)

# load data with original population names ---------------------------------------------------------
library(readxl)

bottleneck_out <- read_excel("data/processed/out_bottleneck_stats.xls")
names(bottleneck_out)[1] <- "id"

# extract pure names
bottleneck_out$id <- sapply(strsplit(bottleneck_out$id, "_genepop"), `[[`, 1)

# check which columns should be numeric
charcols <- str_detect(names(bottleneck_out), "Def.Exc") | str_detect(names(bottleneck_out), "id") | str_detect(names(bottleneck_out), "Mode_Shift")
bottleneck_out[!charcols] <- lapply(bottleneck_out[!charcols], as.numeric)

# split up column with number of loci in het excess
str(bottleneck_out)
exc_cols <- str_detect(names(bottleneck_out), "Def.Exc")

split_up <- function(x){
    df <- apply(data.frame(str_split_fixed(x, "vs", 2)), 2, as.numeric)
    out <- as.data.frame(df)
}

sep_cols <- do.call(cbind, lapply(bottleneck_out[exc_cols], split_up))

names(sep_cols) <- str_replace(names(sep_cols), "X1", "het_def")
names(sep_cols) <- str_replace(names(sep_cols), "X2", "het_exc")
names(sep_cols) <- str_replace(names(sep_cols), "Def.Exc.", "")

# add to original data.frame
bottleneck <- cbind(bottleneck_out[!exc_cols], sep_cols)
# bottleneck <- bottleneck[, !(str_detect(names(bottleneck), "TPM99"))]

# just take datasets in largest cluster cluster
# ids <- c("antarctic_fur_seal", "galapagos_fur_seal", "stellers_sea_lion_cl_2",
#     "grey_seal_orkneys", "harbour_seal_waddensee_cl_2", "galapagos_sea_lion",
#     "south_american_fur_seal_cl_2", "hooded_seal", "mediterranean_monk_seal",
#     "hawaiian_monk_seal", "bearded_seal_cl_2", "crabeater_seal",
#     "leopard_seal", "arctic_ringed_seal", "ross_seal",
#     "weddell_seal_cl_1", "northern_fur_seal_cl_1", "atlantic_walrus_cl_1",
#     "nes_cl_1", "ses_cl_1", "california_sea_lion", "south_american_sea_lion",
#     "new_zealand_sea_lion", "saimaa_ringed_seal_cl_2", "lagoda_ringed_seal",
#     "baltic_ringed_seal", "new_zealand_fur_seal", "australian_fur_seal")

# correct naming of species from the bottleneck program
#bottleneck$id[11] <- "crabeater_seal" # crabeater_seal_recoded to crabeater_seal
bottleneck$id[7] <- "bearded_seal_cl_2" # bearded seal is now cluster 2 --> still has to be changed in the bottleneck tests
# bottleneck$id[33] <- "arctic_ringed_seal" # ringed seal to arctic ringed seal
# check if all names written correctly


###### take full datasets // here is the step where one could sort out the clusters too
ids <- bottleneck$id[!str_detect(bottleneck$id, "cl")]
#ids <- ids[!str_detect(ids, "HW")]
#ids <- ids[!str_detect(ids, "pop")]    

ids %in%  bottleneck$id 

# some name changes due to renaming , also bearded seal is now cl_2 (has to be changed in bottleneck)
unique_seals <- bottleneck[bottleneck$id %in% ids, ]
bottleneck <- unique_seals

# calculate ratio of het exc
prop_het_exc <- function(mut_mod){
    het_exc <- paste0(mut_mod, "_het_exc")
    het_def <- paste0(mut_mod, "_het_def")
    ratio <- bottleneck[[het_exc]] / (bottleneck[[het_exc]] + bottleneck[[het_def]])
}

all_ratios <- data.frame(do.call(cbind, lapply(c("IAM", "TPM70","TPM80", "TPM90",  "SMM"), prop_het_exc)))
names(all_ratios) <- c("IAM_ratio", "TPM70_ratio","TPM80_ratio", "TPM90_ratio",  "SMM_ratio")

bottleneck <- cbind(bottleneck, all_ratios)

wilc_tests <- str_detect(names(bottleneck), "Wilc_Exc")
het_exc_ratios <- str_detect(names(bottleneck), "_ratio")

bottle_tests <- bottleneck[, wilc_tests]
bottle_tests_ratio <- bottleneck[, het_exc_ratios]
bottle_tests$id <- bottleneck$id
bottle_tests_ratio$id <- bottleneck$id


bottleneck_results <- merge(bottle_tests, bottle_tests_ratio, by = "id")
names(bottleneck_results)
bottleneck_final <- bottleneck_results[c(1,7,2,8,3,9,5,10,6,11,4)]

write_delim(bottleneck_final, path = "data/processed/bottleneck_results.txt")


