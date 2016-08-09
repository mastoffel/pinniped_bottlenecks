# script to calculate HW for all seal microsatellite datasets and select HW loci and put these
# into a new excel file as additional datasets

# Hardy Weinberg
library(pegas)
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library(stringr)
library(dplyr)
library(sealABC)

# load cleaned seal data
library(readxl)

seal_data <- "data/processed/seal_data_largest_clust_and_pop.xlsx"
all_seals <- sealABC::read_excel_sheets(seal_data)


# create genind files where genotypes are stored as "134/140" instead of two column format and save
# them in folder ../data/genind_formatted as tab seperated text files, including id and population--------

# where do the genotypes start?
gene_start <- 4
# create file with genotypes stored as "134/140" from two-column per locus format
change_geno_format <- function(seal_data_df) {
    loci_names <- names(seal_data_df[gene_start:ncol(seal_data_df)])
    # filter loci names and remove the allel add-on
    loci_names <- str_replace_all(loci_names[seq(from = 1, to = length(loci_names), by = 2)], "_a", "")
    short_geno <- data.frame(matrix(nrow = nrow(seal_data_df), ncol = length(loci_names)))
    names(short_geno) <- loci_names

    genotypes <- seal_data_df[, gene_start:ncol(seal_data_df)]
    length_data <- ncol(genotypes)

    col_num <- 1
    for (i in seq(from = 1, to = length_data, by = 2)) {
        short_geno[col_num] <- paste(genotypes[[i]], genotypes[[i+1]], sep = "/")
        col_num <- col_num + 1
    }

    short_geno[short_geno == "NA/NA"] <- NA
    geno_out <- cbind(seal_data_df[c(1,2)], short_geno)
    # format pop and ID to be one string
    geno_out$pop <- str_replace_all(as.character(geno_out$pop), ",", "")
    geno_out$pop <- str_replace_all(as.character(geno_out$pop), " ", "_")
    geno_out$id <- str_replace_all(as.character(geno_out$id), " ", "_")
    geno_out
}

seal_data_locus_format <- lapply(all_seals, change_geno_format)



check_na_genotypes <- function(df) {
                                df[3:ncol(df)] <- lapply(df[gene_start:ncol(df)], function(x) {
                                    if (sum(grepl("NA", x))>0) {
                                        x[which(grepl("NA", x))] <- NA
                                    }
                                    x
                                })
                        df
}

# MAIN FORMAT to convert to genind and loci classes
seal_data_locus_format <- lapply(seal_data_locus_format, check_na_genotypes)


# write to txt files
for (i in 1:length(seal_data_locus_format)) {
        write.table(seal_data_locus_format[[i]], file = paste("output/genind_formatted/", names(seal_data_locus_format[i]), ".txt", sep = ""
                                                        ), sep = " ", quote = FALSE, row.names = FALSE)
}
# load it for pegas
seal_data_pegas <- list()
library(pegas)
for (i in 1:length(seal_data_locus_format)) {
seal_data_pegas[[names(seal_data_locus_format)[i]]] <- read.loci(paste("output/genind_formatted/", names(seal_data_locus_format)[i], ".txt", sep = ""),
               header = TRUE, loci.sep = " ", allele.sep = "/", col.pop = 2, col.loci = c(3:ncol(seal_data_locus_format[[i]])))
}

# save(seal_data_pegas, file = "seal_data_pegas.RData")

# create data.frame with all descriptive variables
sample_size <- as.numeric(lapply(seal_data_pegas, function(x) nrow(x)))
loci_full <- as.numeric(as.numeric(lapply(seal_data_pegas, function(x) ncol(x)-2)))
total_genotypes <- sample_size * loci_full

# Exact test of Hardy-Weinberg equilibrium by Markov chain Monte Carlo----------------
# load("seal_data_pegas.RData")

#### hardy weinberg tests on everything
# all_hw <- lapply(seal_data_pegas, hw.test, B = 10000)
# save(all_hw, file = "data/processed/all_hw_10000iter.RData")
load("data/processed/all_hw_10000iter.RData")
all_non_hw <- lapply(all_hw, function(x) {
                                    x <- as.data.frame(x)
                                    x[which(x["Pr.exact"] < 0.05), ]
                                    })
all_non_hw_number <- unlist(lapply(all_non_hw, nrow))


# Chisqu test of Hardy-Weinberg equilibrium by Markov chain Monte Carlo----------------
all_non_hw_chi <- lapply(all_hw, function(x) {
    x <- as.data.frame(x)
    x[which(x["Pr(chi^2 >)"] < 0.05), ]
})
all_non_hw_number_chi <- unlist(lapply(all_non_hw_chi, nrow))

# number of populations
num_pop <- unlist(lapply(seal_data_locus_format, function(x) out <- length(levels(as.factor(x$pop)))))


### with bonferroni correction
bonf_corr <- function(x) {
    x <- as.data.frame(x)
    x$p_exac_bonf <- p.adjust(x$Pr.exact, method = "bonferroni")
    x$chi_bonf <- p.adjust(x[[3]], method = "bonferroni")
    x
}
all_hw_bonf <- lapply(all_hw, bonf_corr)

## get loci in HW according to exact test
bonf_corr_loci <- function(x) {
    x <- as.data.frame(x)
    x$p_exac_bonf <- p.adjust(x$Pr.exact, method = "bonferroni")
    if (any(x$p_exac_bonf < 0.05)){
        loci_non_hw <- which(x$p_exac_bonf < 0.05)
        return(row.names(x)[-loci_non_hw])
    } else {
    return(row.names(x))
    }
}
all_hw_bonf_loci <- lapply(all_hw, bonf_corr_loci)


# how many loci not in hw per test?
non_hw_p_bonf <- unlist(lapply(all_hw_bonf, function(x) sum(x$p_exac_bonf < 0.05)))
non_hw_chi_bonf <- unlist(lapply(all_hw_bonf, function(x) sum(x$chi_bonf < 0.05)))

# overlap both tests
non_hw_both_bonf <- unlist(lapply(all_hw_bonf, function(x) sum((x$p_exac_bonf < 0.05) & (x$chi_bonf < 0.05))))

seal_data_descr <- data.frame("names" = names(all_seals),
                               "sample_size" = sample_size,
                              "loci_full" = loci_full ,
                              "total_genotypes" = total_genotypes,
                              "number_populations" = num_pop,
                              "non_hw_exact_test" = all_non_hw_number,
                              "non_hw_chisqu_test" = all_non_hw_number_chi,
                              "non_hw_exact_bonf" = non_hw_p_bonf,
                            "non_hw_chisqu_bonf" = non_hw_chi_bonf,
                              "non_hw_both_tests" = non_hw_both_bonf)
                              
library(WriteXLS)
WriteXLS(seal_data_descr, "data/processed/seal_data_descriptives.xls")


## get loci in HW according to exact test
bonf_corr_loci <- function(x) {
    x <- as.data.frame(x)
    x$p_exac_bonf <- p.adjust(x$Pr.exact, method = "bonferroni")
    if (any(x$p_exac_bonf < 0.05)){
        loci_non_hw <- which(x$p_exac_bonf < 0.05)
        return(row.names(x)[-loci_non_hw])
    } else {
        return(row.names(x))
    }
}
all_hw_bonf_loci <- lapply(all_hw, bonf_corr_loci)


library(stringr)
extract_loci_in_hw <- function(species, all_hw_bonf_loci, seal_data){
    loci_hw <- all_hw_bonf_loci[[species]]
    if (length(loci_hw) == 0) return(seal_data[[species]][c(1:3)])
    loci_names <- names(seal_data[[species]])
    # problem of multiple matching
    which_match <- function(loci_in_hw) out <- which(!is.na(str_match(loci_names, loci_in_hw)))
    loci_a <- paste0(loci_hw, "_a")
    loci_b <- paste0(loci_hw, "_b")
    
    col_in_hw_a <- unique(as.numeric(unlist(sapply(loci_a, which_match))))
    col_in_hw_b <- unique(as.numeric(unlist(sapply(loci_b, which_match))))
    
    col_in_hw <-  sort(c(col_in_hw_a, col_in_hw_b))
    
    seals <- seal_data[[species]]
    out <- cbind(seals[, c(1:3)], seals[, col_in_hw])
}

species <- names(all_seals)

all_seals_in_hw <- lapply(species , extract_loci_in_hw, all_hw_bonf_loci,  all_seals)
names(all_seals_in_hw) <- species
unlist(lapply(all_seals_in_hw, function(x) ncol(x) - 3)) / 2
# extract get data.frames just with loci in HW

# write excel file with each dataset 
write_dflist_to_xls(all_seals_in_hw, "seal_data_largest_clust_and_pop_all_hw.xls")


####### LAST step ######
# save all xls files as xlsx for postprocessing 
