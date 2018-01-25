# checking all microsatellite datasets and output a basic excel file
# with genotypes, ids and population column for 28 species

library(devtools)
# install_github("mastoffel/sealABC")
library(sealABC)

# load all seal data
all_seals <- sealABC::read_excel_sheets("data/processed/seal_data_largest_clust_and_pop.xlsx")
names(all_seals)

# get all full datasets
all_seals <- all_seals[1:28]

# delete the cluster column
seals <- lapply(all_seals, function(x) x <- x[-3])

# check all ranges
all_ranges <- lapply(seals, function(x) lapply(x[3:ncol(x)], range, na.rm = TRUE))
names(all_ranges)
str(all_ranges)


write_dflist_to_xls(seals, "data/processed/seal_genotypes_basic_28.xls")
