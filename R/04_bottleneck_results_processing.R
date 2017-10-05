# old bottleneck results processing, has to be looked over

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
# sheet numbers to load
bottleneck_out <- read_excel("data/processed/out_bottleneck_stats.xls")
bottleneck_out_hw <- read_excel("data/processed/out_bottleneck_stats_HW.xls")

bottleneck_out <- rbind(bottleneck_out, bottleneck_out_hw)
names(bottleneck_out)[1] <- "id"

# extract pure names
bottleneck_out$id <- sapply(strsplit(bottleneck_out$id, "_genepop"), `[[`, 1)

# numeric
charcols <- str_detect(names(bottleneck_out), "Def.Exc") | str_detect(names(bottleneck_out), "id")
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
bottleneck <- bottleneck[, !(str_detect(names(bottleneck), "TPM99"))]

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
bottleneck$id[11] <- "crabeater_seal" # crabeater_seal_recoded to crabeater_seal
bottleneck$id[7] <- "bearded_seal_cl_2" # bearded seal is now cluster 2 --> still has to be changed in the bottleneck tests
bottleneck$id[33] <- "arctic_ringed_seal" # ringed seal to arctic ringed seal
# check if all names written correctly


###### take full datasets
ids <- bottleneck$id[!str_detect(bottleneck$id, "cl")]
ids <- ids[!str_detect(ids, "HW")]
ids <- ids[!str_detect(ids, "pop")]    

ids %in%  bottleneck$id 

# some name changes due to renaming , also bearded seal is now cl_2 (has to be changed in bottleneck)
unique_seals <- bottleneck[bottleneck$id %in% ids, ]
bottleneck <- unique_seals

# change names to unique names

# dont run for full datasets
# add factor for full / pop / cl
# bottleneck$dataset[str_detect(bottleneck$id, "_pop")] <- "pop"
# bottleneck$dataset[str_detect(bottleneck$id, "_cl")] <- "cl"
# bottleneck$dataset[is.na(bottleneck$dataset)] <- "full"

prop_het_exc <- function(mut_mod){
    het_exc <- paste0(mut_mod, "_het_exc")
    het_def <- paste0(mut_mod, "_het_def")
    ratio <- bottleneck[[het_exc]] / (bottleneck[[het_exc]] + bottleneck[[het_def]])
}

all_ratios <- data.frame(do.call(cbind, lapply(c("IAM", "TPM70","TPM90", "TPM95",  "SMM"), prop_het_exc)))
names(all_ratios) <- c("IAM_ratio", "TPM70_ratio","TPM90_ratio", "TPM95_ratio",  "SMM_ratio")

bottleneck <- cbind(bottleneck, all_ratios)

wilc_tests <- str_detect(names(bottleneck), "Wilc_Exc")
het_exc_ratios <- str_detect(names(bottleneck), "_ratio")

bottle_tests <- bottleneck[, wilc_tests]
bottle_tests_ratio <- bottleneck[, het_exc_ratios]
bottle_tests$id <- bottleneck$id
bottle_tests_ratio$id <- bottleneck$id


bottleneck_results <- merge(bottle_tests, bottle_tests_ratio, by = "id")
write.xlsx(bottleneck_results, file = "data/processed/bottleneck_results.xlsx", row.names = FALSE)
# save(bottle_tests, bottle_tests_ratio, file = "bottleneck_results.RData")

bot <- melt(bottle_tests, id.vars = "id")
bot$value <- as.numeric(bot$value)

names(bot) <- c("id", "mutation_model", "value")

ggplot(bot, aes(x= mutation_model, y = id, fill = value)) + 
            #facet_grid(.~dataset) + 
            geom_tile(color = "white", size = 0.1) +
           # scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#1d91c0", "#225ea8", "#0c2c84", "#081d58")) +
            #scale_fill_gradientn(colours=c("#fff7fb", "#d0d1e6", "#67a9cf", "#02818a", "#014636")) +
            #scale_fill_gradientn(colours=c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#cb181d", "#a50f15", "#67000d")) +
            scale_fill_gradientn(colours=c("#ffffd9", "#edf8b1", "#c7e9b4", "#41b6c4", "#225ea8", "#253494", "#081d58")) +
            # scale_fill_gradient(name = "p-value", label = comma,  breaks = c(0.05, 0.5)) +
            theme_tufte(base_family="Helvetica") +
            coord_equal() +
            theme(plot.title=element_text(hjust=0),
                  axis.ticks=element_blank(),
                  axis.text=element_text(size=10),
                  legend.title=element_text(size=10),
                  legend.text=element_text(size=9)) 
           # scale_fill_gradientn(colours=c("darkred", "red", "orange", "yellow"),                      
        #    values  = rescale(c(min(bot$value), 0.05, 0.10, max(bot$value))))
            # coord_equal()
             # scale_fill_gradientn(colors =  breaks=c(0.01, 0.05, 0.1, 1))
   
            #coord_equal()


