
library(readxl)
# sheet numbers to load
data <- "data/processed/seal_data_largest_clust_and_pop.xlsx"
dataset_names <- excel_sheets(data)

load_dataset <- function(dataset_names) {
    read_excel(data, sheet = dataset_names)
}
# load all datasets
seal_data <- lapply(dataset_names, load_dataset)
names(seal_data) <- dataset_names
# afs
genotypes <- seal_data[[1]]
genotypes <- genotypes[4:ncol(genotypes)]


library(strataG)
#### calculate empirical summary statistics
g_types_geno <- new("gtypes", genotypes, ploidy = 2)
# summary statistics
# according to DIYabc
# number of alleles across loci
num_alleles <- strataG::numAlleles(g_types_geno)
num_alleles_mean <- mean(num_alleles, na.rm = TRUE)
num_alleles_sd <- sd(num_alleles, na.rm = TRUE)

# allele size variance (actually, sd) across loci
allele_size_sd <- unlist(lapply(seq(from = 1, to = ncol(genotypes), by = 2), 
    function(x) sd(unlist(genotypes[x:(x+1)]), na.rm = TRUE)))

mean_allele_size_sd <- mean(allele_size_sd, na.rm = TRUE)
sd_allele_size_sd <- sd(allele_size_sd, na.rm = TRUE)

# expected heterozygosity
exp_het <- exptdHet(g_types_geno)
exp_het_mean <- mean(exp_het, na.rm = TRUE)
exp_het_sd<- sd(exp_het, na.rm = TRUE)
# observed heterozygosity
obs_het <- obsvdHet(g_types_geno)
obs_het_mean <- mean(obs_het, na.rm = TRUE)
obs_het_sd <- sd(obs_het, na.rm = TRUE)

obs_stats <- data.frame(num_alleles_mean, num_alleles_sd, mean_allele_size_sd, sd_allele_size_sd,
    exp_het_mean, exp_het_sd, obs_het_mean, obs_het_sd)




## calculate same statistics with diversity
library(diveRsity)
genepop_afs <- "data/processed/genepop_files/antarctic_fur_seal_genepop.txt"
test <- read.table(genepop_afs)


out <- divBasic(infile = genepop_afs)
out <- basicStats(genepop_afs, rarefaction = FALSE)

out <- readGenepop(genepop_afs)

library("diveRsity")
data(Test_data)

# divide stats and parameters
sims_stats <- sims[6:(ncol(sims)-1)] # last column is model factor
sims_param <- sims[1:5]
models <- sims[["model"]]

library(abc)

# abc
cv.modsel <- cv4postpr(models, sims_stats, nval=5, tol=.25, method="rejection")
s <- summary(cv.modsel)
plot(cv.modsel, names.arg=c("bot", "neut", "decl"))

mod_prob <- postpr(obs_stats, models, sims_stats, tol = 0.3, method = "mnlogistic")
summary(mod_prob)


res.gfit.bott <- gfit(target=obs_stats, sumstat= sims_stats[models == "bot",], statistic=mean, nb.replicate=100)
plot(res.gfit.bott, main="Histogram under H0")

summary(res.gfit.bott)

# or PCA

gfitpca(target=obs_stats, sumstat=sims_stats, index=models, cprob=.1)

## 
