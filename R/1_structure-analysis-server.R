# This script runs from a server where STRUCTURE is installed.
# It uses the parallelStructure package to interface with STRUCTURE.

library(readxl)
library(ParallelStructure)
library(sealABC)

## load a list of genetic data frames and runs STRUCTURE in parallel on all of them
seal_data <- "data/processed/seal_data_largest_clust_and_pop.xlsx"
all_seals <- sealABC::read_excel_sheets(seal_data)

# just get the 28 single species
all_seals <- all_seals[1:28]

# prepare data
all_seals_struc <- lapply(all_seals, function(x){
                             x[is.na(x)] <- -9
                             x
                    })

# format all seal data to just one population
seals_no_pop <- lapply(all_seals_struc, function(x) {
                            x$pop <- 1
                            x
                    })
seals_no_clus <- lapply(seals_no_pop, function(x){
                            x$cluster <- NULL
                            x
} )

options(scipen=999)
# format

# construct job matrix and write to job file
nrep <- 5
burnin <- 10000
niter <- 100000
up_to_k <- 10

# define variables for job matrix
k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

# put matrix together --> all populations are seen as one!
pop <- 1
seal_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k), 
                    rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)
# write job file
write(t(seal_jobs), ncol=length(seal_jobs[1,]), file='seal_jobs.txt')

# file path to structure
STR_path='/usr/local/bin/'

# prepare runs
system("mkdir structure_results")
# transform into char
seals <- lapply(seals_no_clus, function(x) as.matrix(apply(x, 2, as.character)))


for (i in 1:length(seals)){
    seal_struc <- seals[[i]]
    # create outfile
    system(paste0("mkdir structure_results/", names(all_seals)[i]))
    outpath <- paste0("structure_results/", names(all_seals)[i], "/")
    filename <- paste0(names(all_seals)[i], ".txt")
    write(t(seal_struc), ncol=length(seal_struc[1,]), file= filename)
    
    ParallelStructure::parallel_structure(structure_path=STR_path, joblist='seal_jobs.txt', n_cpu=20, infile=filename, outpath= outpath, 
                        numinds=nrow(seal_struc), numloci=(ncol(seal_struc)-2)/2, printqhat=1, plot_output = 1, onerowperind=1)
}

