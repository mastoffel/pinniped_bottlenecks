# This script cleans all ice breeding and ringed seal datasets by exchanging
# the 900 + genotypes. For the species from Davis (2008), it will remove the leading
# 9. For the species from Nyman (2014), i.e. the 3 ringed seal subspecies, it willö
# remove the last digit.


# load all seal data
all_seals <- sealABC::read_excel_sheets("data/processed/seal_data_largest_clust_and_pop_all_hw.xlsx")

lapply(all_seals, function(x) range(unlist(x[4:ncol(x)]), na.rm = TRUE))


## Discard the first 9 for all of Davis´s datasets
change_900_alleles_disc1st <- function(x, keep1, keep2){ # which digits to keep
    if (any(x > 899, na.rm = TRUE)) {
        x <- as.character(x)
        x[(nchar(x, keepNA = FALSE) == 3) & (substr(x, 1, 1) == "9")] <- substr(x[(nchar(x, keepNA = FALSE) == 3) & (substr(x, 1, 1) == "9")], keep1, keep2)
        x <- as.numeric(x)
    }
    x
}

change_900_alleles_df_disc1st <- function(df,  keep1, keep2){
    out <- as.data.frame(do.call(cbind, lapply(df[4:ncol(df)], change_900_alleles_disc1st, keep1, keep2)))
    out <- cbind(df[1:3], out)
    out
}

clean_DAVISandNyman_data <- function(x){
    out <- all_seals[[x]]
    if (str_detect(names(all_seals)[x], paste(c("bearded", "crabeater", "leopard", "arctic_ringed", "ross", "weddell"), collapse = "|"))){
        out <- change_900_alleles_df_disc1st(all_seals[[x]], 2, 3)
    }
    if (str_detect(names(all_seals)[x], paste(c("lagoda", "saimaa", "baltic"), collapse = "|"))){ ### care here as lagoda might be correctly renamed to ladoga at some point
        out <- change_900_alleles_df_disc1st(all_seals[[x]], 1, 2)
    }
    out
}

all_seals_new <- lapply(c(1:length(all_seals)), clean_DAVISandNyman_data)

lapply(all_seals_new, function(x) range(unlist(x[4:ncol(x)]), na.rm = TRUE))
names(all_seals_new) <- names(all_seals)


# write excel file with each dataset 
library(WriteXLS)
envir <- environment()
list_to_df <- function(species, dfs, envir){
    assign(species, dfs[[species]], envir)
}
lapply(names(all_seals_new), list_to_df, all_seals_new, envir)
WriteXLS(names(all_seals_new), ExcelFileName = "seal_data_largest_clust_and_pop_all_hw.xls")

