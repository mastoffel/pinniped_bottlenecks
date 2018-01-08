library(readr)
library(readxl)
library(dplyr)
library(knitr)
library(kableExtra)
# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_29_modeling.csv"))

# number of genotypes
sum(all_stats$nind * all_stats$nloc)

# Supplementary Table 1: IUCN, life history data and sample size/locus information
all_stats_origin <- read_xlsx("data/raw/table_data.xlsx")
all_stats_table <- all_stats %>% 
    left_join(all_stats_origin, by = c("species", "common", "tip_label", "latin")) %>% 
    dplyr::mutate(Genotypes = nloc * nind) %>% 
    dplyr::select(common, latin, IUCN_rating, Abundance, BreedingType, 
        SSD, harem_size, nloc, nind, Genotypes, origin, published_short) %>% 
    dplyr::rename(`Common name` = common,
        `Scientific` = latin,
        `IUCN status` = IUCN_rating,
        `Breeding habitat` = BreedingType,
        `Harem size` = harem_size,
        `Loci` = nloc,
        `Individuals` = nind,
        `Sampling location` = origin,
        `Published` = published_short) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(Scientific = cell_spec(Scientific, italic = TRUE))

options(knitr.table.format = "latex")

kable(all_stats_table , format = "latex",  escape = F, 
      booktabs = TRUE, align = "c", digits = 3, linesep = "") %>% 
    kable_styling(latex_options = c( "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("my_latex_table", keep_pdf = TRUE)

# Supplementary Table 2: genetic, bottleneck data
# Table 2: Genetic information and outcome of the ABC
transform_ss <- function(PE, CIlow, CIhigh) {
    out <- paste0(PE, " (", CIlow," ,", CIhigh, ")")
    out
}

all_stats_table2 <- all_stats %>% 
    dplyr::select(common, latin,
        num_alleles_mean, num_alleles_mean.CIlow, num_alleles_mean.CIhigh,
        obs_het_mean, obs_het_mean.CIlow, obs_het_mean.CIhigh,
        exp_het_mean, exp_het_mean.CIlow, exp_het_mean.CIhigh,
        prop_low_afs_mean, prop_low_afs_mean.CIlow, prop_low_afs_mean.CIhigh,
        mean_allele_range, mean_allele_range.CIlow, mean_allele_range.CIhigh,
        mratio_mean, mratio_mean.CIlow, mratio_mean.CIhigh
    ) %>% 
    mutate_if(is.numeric, funs(round(., 2))) 

# bin ci and means together
all_stats_table2_short <- data.frame(do.call(cbind, lapply(seq(from = 3, to = ncol(all_stats_table2 ), by = 3), 
    function(x) transform_ss(all_stats_table2[[x]], all_stats_table2[[x + 1]], 
        all_stats_table2[[x + 2]])))) %>% 
    dplyr::rename(`Allelic richness` = X1,
        `Obs. heterozygosity` = X2,
        `Exp. heterozygosity` = X3,
        `Prop. low frequency alleles` = X4,
        `Allelic range` = X5,
        `M-ratio` = X6) %>% 
    bind_cols(all_stats_table2[1:2], .) %>% 
    dplyr::arrange(-row_number())




all_stats_table <- all_stats %>% 
    left_join(all_stats_origin, by = c("species", "common", "tip_label", "latin")) %>% 
    dplyr::select(common, latin) %>% 
    dplyr::rename(`Common name` = common,
        `Scientific` = latin,
        `IUCN status` = IUCN_rating,
        `Breeding habitat` = BreedingType,
        `Harem size` = harem_size,
        `Loci` = nloc,
        `Individuals` = nind,
        `Sampling location` = origin,
        `Published` = published_short) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(Scientific = cell_spec(Scientific, italic = TRUE))

