library(readr)
library(readxl)
library(dplyr)
library(knitr)
# create tables for supplementary material
library(kableExtra)
# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_29_modeling.csv"))
all_stats$latin[7] <- "Otaria byronia / flavescens"
# number of genotypes
sum(all_stats$nind * all_stats$nloc)

# Supplementary Table 1: IUCN, life history data and sample size/locus information
all_stats_origin <- read_xlsx("data/raw/table_data.xlsx")
all_stats_table <- all_stats %>% 
    left_join(all_stats_origin, by = c("species", "common", "tip_label", "latin")) %>% 
    dplyr::mutate(Genotypes = nloc * nind) %>% 
    dplyr::select(common, latin, IUCN_rating, Abundance, BreedingType, Generation_time,
        SSD, harem_size, nloc, nind, Genotypes, origin, published_short) %>% 
    dplyr::rename(`Common name` = common,
        `Scientific` = latin,
        `IUCN status` = IUCN_rating,
        `Breeding habitat` = BreedingType,
        `Generation time` = Generation_time,
        `Harem size` = harem_size,
        `Loci` = nloc,
        `Individuals` = nind,
        `Sampling location` = origin,
        `Published` = published_short) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(Scientific = cell_spec(Scientific, italic = TRUE))

options(knitr.table.format = "latex")

align_tab1 <- c("l", "l", rep(x = "c", ncol(all_stats_table) - 2))

kable(all_stats_table , format = "latex",  escape = F, 
      booktabs = TRUE, align = align_tab1, digits = 3, linesep = "") %>% 
    kable_styling(latex_options = c( "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("my_latex_table", keep_pdf = TRUE)


# Supplementary Table 2: genetic, bottleneck data
# Table 2: Genetic summary statistics per 10 ind.
transform_ss <- function(PE, CIlow, CIhigh) {
    out <- paste0(PE, " (", CIlow,", ", CIhigh, ")")
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
all_stats_table2_short <- data.frame(
    do.call(cbind, lapply(seq(from = 3, to = ncol(all_stats_table2 ), by = 3), 
    function(x) transform_ss(all_stats_table2[[x]], all_stats_table2[[x + 1]], 
        all_stats_table2[[x + 2]])))) %>% 
    bind_cols(all_stats_table2[1:2], .) %>% 
    dplyr::rename(
        `Common name` = common,
        `Scientific` = latin,
        `Allelic richness` = X1,
        `Obs. heterozygosity` = X2,
        `Exp. heterozygosity` = X3,
        `Prop. low frequency alleles` = X4,
        `Allelic range` = X5,
        `M-ratio` = X6) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(Scientific = cell_spec(Scientific, italic = TRUE))

align_tab2 <- c("l", "l", rep(x = "c", ncol(all_stats_table2_short ) - 2))

kable(all_stats_table2_short, format = "latex",  escape = F, 
    booktabs = TRUE, align = align_tab2, digits = 3, linesep = "") %>% 
    kable_styling(latex_options = c( "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("Sup_tab_2", keep_pdf = TRUE)

# supplementary table 3: Bottleneck stats

all_stats_table_sub3 <- all_stats %>% 
    left_join(all_stats_origin, by = c( "common",  "latin")) %>% 
    dplyr::select(4, 3, 78, 80, 82, 84, 89, 90) %>% 
    dplyr::rename(`Common name` = common,
        `Scientific` = latin,
        `TPM 70` =  TPM70_ratio,
        `TPM 80` = TPM80_ratio,
        `TPM 90` = TPM90_ratio,
        `SMM` = SMM_ratio,
        `$\\boldsymbol{p_{bot}}$` = bot,
        `$\\boldsymbol{p_{neut}}$` = neut) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(Scientific = cell_spec(Scientific, italic = TRUE))

align_tab3 <- c("l", "l", rep(x = "c", ncol(all_stats_table_sub3 ) - 2))

kable(all_stats_table_sub3, format = "latex", escape = F, 
    booktabs = TRUE, align = align_tab3, digits = 3, linesep = "") %>% 
    add_header_above(c(" "," ", "Heterozygosity-excess ($prop_{het-exc}$)" = 4, "ABC" = 2), escape = F) %>%
    kable_styling(latex_options =  "scale_down") %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("Sup_tab_3", keep_pdf = TRUE)


# Supplementary table 3: Model fit ABC
model_fit <- read_delim("output/model_evaluation/check3_modeval/sims_10000k_p_vals_fit.txt", 
                        delim = " ") 

all_stats_table_sub4_temp <- all_stats %>% 
    left_join(all_stats_origin, by = c( "common",  "latin", "species")) %>% 
    full_join(model_fit, by = c("species")) 
 
all_stats_table_sub4 <-all_stats_table_sub4_temp %>% 
    dplyr::select("common", "latin", "bot") %>% 
    dplyr::mutate(model = ifelse(bot > 0.5, "bot", "neut")) %>% 
    dplyr::mutate(p_val = ifelse(bot > 0.5, all_stats_table_sub4_temp$bot_p, 
                                             all_stats_table_sub4_temp$neut_p)) %>% 
    dplyr::select(-bot) %>% 
    dplyr::rename(`Common name` = common,
        `Scientific` = latin,
        `Selected model` =  model,
         `p value` = p_val) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(Scientific = cell_spec(Scientific, italic = TRUE))

align_tab4 <- c("l", "l", rep(x = "c", ncol(all_stats_table_sub4 ) - 2))

kable(all_stats_table_sub4, format = "latex", escape = F, 
    booktabs = TRUE, align = align_tab4 , digits = 3, linesep = "") %>% 
    kable_styling(latex_options =  "scale_down") %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("Sup_tab_4", keep_pdf = TRUE)
