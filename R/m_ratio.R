# 
library(readxl)
library(strataG)

# sheet numbers to load
dataset_names <- excel_sheets("data/processed/seal_data_largest_clust_and_pop_all_hw.xlsx")

load_dataset <- function(dataset_names) {
    read_excel("data/processed/seal_data_largest_clust_and_pop_all_hw.xlsx", sheet = dataset_names)
}
# load all datasets
all_seals <- lapply(dataset_names, load_dataset)
names(all_seals) <- dataset_names

# delete genotypes with weird ranges (allele sizes over 900)
lapply(all_seals, function(x) range(x[, 4:ncol(x)], na.rm = TRUE))

exclude_weird_genotypes <- function(x) {
    x
    new_geno <- do.call(cbind, lapply(x[4:ncol(x)], function(y){
                        y[y > 899] <- NA
                        y
                        }))
    out <- cbind(x[, 1:3], new_geno)
}

all_seals <- lapply(all_seals, exclude_weird_genotypes)

# mratio

mean_mRatio <- function(genotypes){
    msats <- new("gtypes", gen.data = genotypes[4:ncol(genotypes)], ploidy = 2)
    out <- mean(mRatio_new(msats, prop_corr_rep_length = 0.3), na.rm = TRUE)
}

all_mratios3 <- lapply(all_seals, mean_mRatio)
all_mratios2

msats <- new("gtypes", gen.data = afs[4:ncol(afs)], ind.names = afs$id, ploidy = 2)
mRatio(msats)



# test whether strataG mratio is correct
afs <- all_seals$antarctic_fur_seal
afs <- afs[-c(1:3)]
head(afs)


# calculate mRatios myself
calc_mratio_dataset <- function(genotypes){
    
    genotypes <- genotypes[4:ncol(genotypes)]
    
    loci <- seq(from = 1, to = ncol(genotypes), by = 2)
    
    calc_mratio_locus <- function(locus_ind, genotypes) {
        range_allele_size <- range(genotypes[c(locus_ind, locus_ind+1)], na.rm = TRUE)
        diff_as <- range_allele_size[2] - range_allele_size[1]
        number_of_alleles <- ncol(table(genotypes[c(locus_ind, locus_ind+1)]))
        out <- number_of_alleles / diff_as
    }
    
    all_mratios <- (unlist(lapply(loci, calc_mratio_locus, genotypes)))
    all_mratios
}

all_ratios <- lapply(all_seals, calc_mratio_dataset)


hist(unlist(lapply(all_ratios, mean)))



# Start writing to an output file
out <- system("echo '246 21 100 3.5 0.3' | R/Critical_M", intern = TRUE)



sink() 
sink(type="message")
cat(readLines("test.log"), sep="\n")

system("./Critical_M")
system("echo Y | 5")
system("30")


data(dolph.msats)
data(dolph.strata)
strata.schemes <- dolph.strata[, c("broad", "fine")]
rownames(strata.schemes) <- dolph.strata$id


msats.g <- new("gtypes", gen.data = dolph.msats[, -1], ploidy = 2,
    ind.names = dolph.msats[, 1], schemes = strata.schemes)

library(strataG)
browseVignettes("strataG")

mRatio()