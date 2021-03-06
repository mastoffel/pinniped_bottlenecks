# Second processing step of BOTTLENECK output
# calculation of prop. loci in het-exc etc.

# Bottleneck output files 
# full data (incl non-HW loci)

# files needed for this script:
# (a) out_bottleneck_stats_30.xls containing the bottleneck output for full data and largest clusters 
# this file is provided, as the BOTTLENECK program was used to calculate all heterozygosity-excess tests.
# Then the script bottleneck_results.R has was run to transform the BOTTLENECK results into 
# a tidy data frame.

# file output from this script:
# (a) bottleneck_results_30.txt formatted bottleneck results for full data sets
# (b) bottleneck_results_30_cl.txt formatted bottleneck results for largest genetic clusters
# (c) if selected: bottleneck_results_HW_30.txt formatted bottleneck results for full data sets containing 
# only loci in HW equilibrium


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
library(readr)
library(readxl)

# load data with original population names ---------------------------------------------------------

# this is the output of the bottleneck tests, merged into one data.frame

# calculate everything on HW filtered bottleneck stats?
HW <- TRUE

if (HW) {
    name_short <- "_HW"
} else {
    name_short <- ""
}

bottleneck_out <- read_excel(paste0("data/processed/out_bottleneck_stats", name_short, "_30.xls"))
names(bottleneck_out)[1] <- "id"

# extract pure names
bottleneck_out$id <- sapply(strsplit(bottleneck_out$id, "_genepop"), `[[`, 1)
bottleneck_out$id <- str_replace_all(bottleneck_out$id, "_HW", "")

# check which columns should be numeric
charcols <- str_detect(names(bottleneck_out), "Def.Exc") | str_detect(names(bottleneck_out), "id") | str_detect(names(bottleneck_out), "Mode_Shift")
bottleneck_out[!charcols] <- lapply(bottleneck_out[!charcols], as.numeric)

# here, we calculate the proportion of loci in heterozygosity-excess by
# splitting up the Def.Exc column with number of loci in het excess
str(bottleneck_out)
exc_cols <- str_detect(names(bottleneck_out), "Def.Exc")

split_up <- function(x){
    df <- apply(data.frame(str_split_fixed(x, "vs", 2)), 2, as.numeric)
    out <- as.data.frame(df)
}
sep_cols <- do.call(cbind, lapply(bottleneck_out[exc_cols], split_up))

# rename splitted columns
names(sep_cols) <- str_replace(names(sep_cols), "X1", "het_def")
names(sep_cols) <- str_replace(names(sep_cols), "X2", "het_exc")
names(sep_cols) <- str_replace(names(sep_cols), "Def.Exc.", "")

# add the het.exc data to original data.frame
bottleneck <- cbind(bottleneck_out[!exc_cols], sep_cols)

# correct naming of species from the bottleneck program
# this has to be double checked in case the clusters are taken
#bottleneck$id[7] <- "bearded_seal_cl_2" # bearded seal is now cluster 2 
#if (HW) bottleneck$id[7] <- "bearded_seal_cl_2_HW"

# (a) just take the full datasets (alternatively, one can get the largest cluster datasets too here)
ids <- bottleneck$id[!str_detect(bottleneck$id, "cl")]

# (b) take full datasets for unclustered data and clustered datasets for clustered data

ids_cl_temp <- bottleneck$id[str_detect(bottleneck$id, "cl")]
ids_cl_fullname <- str_replace(ids_cl_temp, "_cl_[1-9]", "")
# delete 
ids_cl <- ids[!(ids %in% ids_cl_fullname)]
# put together 
ids_cl <- c(ids_cl, ids_cl_temp)



format_bottleneck <- function(ids){ 

    # some name changes due to renaming , also bearded seal is now cl_2 (has to be changed in bottleneck)
    unique_seals <- bottleneck[bottleneck$id %in% ids, ]
    bottleneck <- unique_seals
    
    # calculate proportion of loci in heterozygosity-excess
    prop_het_exc <- function(mut_mod){
        het_exc <- paste0(mut_mod, "_het_exc")
        het_def <- paste0(mut_mod, "_het_def")
        ratio <- bottleneck[[het_exc]] / (bottleneck[[het_exc]] + bottleneck[[het_def]])
    }
    
    all_ratios <- data.frame(do.call(cbind, lapply(c("IAM", "TPM70","TPM80", "TPM90",  "SMM"), prop_het_exc)))
    names(all_ratios) <- c("IAM_ratio", "TPM70_ratio","TPM80_ratio", "TPM90_ratio",  "SMM_ratio")
    
    # bind ratios to data.frame
    bottleneck <- cbind(bottleneck, all_ratios)
    # extract wilcinson tests and ratio
    wilc_tests <- str_detect(names(bottleneck), "Wilc_Exc")
    het_exc_ratios <- str_detect(names(bottleneck), "_ratio")
    
    bottle_tests <- bottleneck[, wilc_tests]
    bottle_tests_ratio <- bottleneck[, het_exc_ratios]
    bottle_tests$id <- bottleneck$id
    bottle_tests_ratio$id <- bottleneck$id
    
    # merge all results
    bottleneck_results <- merge(bottle_tests, bottle_tests_ratio, by = "id")
    names(bottleneck_results)
    bottleneck_final <- bottleneck_results[c(1,7,2,8,3,9,5,10,6,11,4)]
}

# format bottleneck results
out <- apply(cbind(ids, ids_cl), 2, format_bottleneck)
bottleneck_final_full <- out[[1]]
bottleneck_final_cl <- out[[2]]
# for some reason, the bearded seal is cl 1 not 2. 
#bottleneck_final_cl[which(bottleneck_final_cl$id == "bearded_seal_cl_2"), "id"] <- "bearded_seal_cl_1"

# check that bottleneck_final_XX are sorted in the same way
unlist(sapply(bottleneck_final_cl$id, function(x) which(bottleneck_final_full$id %in% str_replace(x, "_cl_[1-9]", ""))))
# double check
data.frame(bottleneck_final_cl$id, bottleneck_final_full$id)

# save
write_delim(bottleneck_final_full, path = paste0("data/processed/bottleneck_results", name_short, "_30.txt"))
write_delim(bottleneck_final_cl, path = "data/processed/bottleneck_results_30_cl.txt")

