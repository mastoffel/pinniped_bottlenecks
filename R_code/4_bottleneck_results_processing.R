# plotting bottleneck

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

# just take datasets in HW and largest cluster
ids <- c("weddell_seal_cl_2",
         "stellers_sea_lion_cl_2_HW",
         "south_american_sea_lion",
         "south_american_fur_seal_cl_2",
         "ses_cl_1_HW",
         "saimaa_ringed_seal_cl_2_HW",
         "ross_seal",
         "ringed_seal",
         "northern_fur_seal_cl_1",
         "new_zealand_sea_lion_cl_1_HW",
         "new_zealand_fur_seal_HW",
         "nes_cl_2_HW",
         "mediterranean_monk_seal_HW",
         "leopard_seal_HW",
         "lagoda_ringed_seal_HW",
         "hooded_seal",
         "hawaiian_monk_seal_HW",
         "harbour_seal_waddensee_cl_1",
         "grey_seal_orkneys_HW",
         "galapagos_sea_lion",
         "galapagos_fur_seal",
         "crabeater_seal_recoded",
         "california_sea_lion_HW",
         "bearded_seal_cl_1",
         "baltic_ringed_seal_HW",
         "australian_fur_seal",
         "atlantic_walrus_schafer_cl_1_HW",
         "antarctic_fur_seal")


# check if all names written correctly
ids %in% bottleneck$id

unique_seals <- bottleneck[bottleneck$id %in% ids, ]
bottleneck <- unique_seals

# change names to unique names

simple_names <- function(seals) {
    seals$id
}

# add factor for full / pop / cl
bottleneck$dataset[str_detect(bottleneck$id, "_pop")] <- "pop"
bottleneck$dataset[str_detect(bottleneck$id, "_cl")] <- "cl"
bottleneck$dataset[is.na(bottleneck$dataset)] <- "full"

prop_het_exc <- function(mut_mod){
    het_exc <- paste0(mut_mod, "_het_exc")
    het_def <- paste0(mut_mod, "_het_def")
    ratio <- bottleneck[[het_exc]] / (bottleneck[[het_exc]] + bottleneck[[het_def]])
}

all_ratios <- data.frame(do.call(cbind, lapply(c("IAM", "TPM70","TPM90", "TPM95",  "SMM"), prop_het_exc)))
names(all_ratios) <- c("IAM_ratio", "TPM70_ratio","TPM90_ratio", "TPM95_ratio",  "SMM_ratio")

bottleneck <- cbind(bottleneck, all_ratios)

# bottle_tests <- bottleneck %>%
#                     mutate(IAM_het_exc_ratio = IAM_Heq / IAM_het_def) %>%
#                     mutate(TPM70_het_exc_ratio = TPM70_Heq / TPM70_het_def) %>%
#                     mutate(TPM95_het_exc_ratio = TPM95_Heq / TPM95_het_def) %>%
#                     mutate(TPM99_het_exc_ratio = TPM99_Heq / TPM99_het_def) %>%
#                     mutate(SMM_het_exc_ratio = SMM_Heq / SMM_het_def) 

# extract sign tests for all models

wilc_tests <- str_detect(names(bottleneck), "Wilc_Exc")
het_exc_ratios <- str_detect(names(bottleneck), "_ratio")

bottle_tests <- bottleneck[, wilc_tests]
bottle_tests_ratio <- bottleneck[, het_exc_ratios]
bottle_tests$id <- bottleneck$id
bottle_tests_ratio$id <- bottleneck$id

save(bottle_tests, bottle_tests_ratio, file = "bottleneck_results.RData")

library(reshape2)
library(ggplot2)
library(dplyr)
bot <- melt(bottle_tests, id.vars = "id")
bot$value <- as.numeric(bot$value)

names(bot) <- c("id", "mutation_model", "value")

library(ggthemes)
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


