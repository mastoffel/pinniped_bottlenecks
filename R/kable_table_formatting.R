library(readr)
library(readxl)
library(dplyr)
library(knitr)
# create tables for supplementary material
library(kableExtra)
options(knitr.table.format = "latex")
# all_stats for modeling
all_stats <- as.data.frame(read_csv("data/processed/all_stats_30_modeling.csv"))
all_stats$latin[7] <- "Otaria byronia / flavescens"
# number of genotypes
sum(all_stats$nind * all_stats$nloc)

# Supplementary Table 1: IUCN, life history data and sample size/locus information
all_stats_origin <- read_xlsx("data/raw/table_data.xlsx")
all_stats_table <- all_stats %>% 
    left_join(all_stats_origin, by = c("species", "common", "tip_label", "latin")) %>% 
    dplyr::mutate(Genotypes = nloc * nind) %>% 
    dplyr::select(common_abbr, latin, IUCN_rating, Abundance, BreedingType, Generation_time,
        SSD, harem_size, nind, nloc, Genotypes,  origin, published_short) %>% 
    dplyr::mutate(Abundance = prettyNum(Abundance, big.mark = ",", scientific = FALSE),
                  Genotypes = prettyNum(Genotypes, big.mark = ",", scientific = FALSE)) %>% 
    dplyr::rename(`Common name` = common_abbr,
        `Scientific name` = latin,
        `IUCN status` = IUCN_rating,
        `Breeding habitat` = BreedingType,
        `Generation time (days)` = Generation_time,
        `Harem size` = harem_size,
        `Loci` = nloc,
        `Individuals` = nind,
        `Sampling location(s)` = origin,
        `Publication` = published_short) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(`Scientific name` = cell_spec(`Scientific name`, format = "latex", italic = TRUE))

#options(knitr.table.format = "latex")

align_tab1 <- c("l", "l", rep(x = "c", ncol(all_stats_table) - 2))

kable(all_stats_table , format = "latex",  escape = F, 
      booktabs = TRUE, align = align_tab1, digits = 3, linesep = "") %>% 
    kable_styling(latex_options = c( "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    add_header_above(c(" " = 2, "Conservation, demography, ecology and life-history data" = 6,
                     "Genetic data" = 5), italic = TRUE) %>% 
    kable_as_image("other_stuff/tables/SupTab1_updated12052018", keep_pdf = TRUE)


# Supplementary Table 2: genetic, bottleneck data
# Table 2: Genetic summary statistics per 10 ind.
transform_ss <- function(PE, CIlow, CIhigh) {
    out <- paste0(PE, " (", CIlow,", ", CIhigh, ")")
    out
}

all_stats_table2 <- all_stats %>% 
    dplyr::select(common_abbr, latin,
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
        `Common name` = common_abbr,
        `Scientific name` = latin,
        `Ar` = X1,
        `Ho` = X2,
        `He` = X3,
        `Proportion of low frequency alleles` = X4,
        `Allelic range` = X5,
        `M-ratio` = X6) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate( `Scientific name` = cell_spec( `Scientific name`, italic = TRUE))

align_tab2 <- c("l", "l", rep(x = "c", ncol(all_stats_table2_short ) - 2))

kable(all_stats_table2_short, format = "latex",  escape = F, 
    booktabs = TRUE, align = align_tab2, digits = 3, linesep = "", #) %>% 
    col.names = c("Common name", "Scientific name", "Ar", "Ho", "He",
                   "Proportion of low frequency alleles",
                   "Allelic range", "M-ratio")) %>% 
    kable_styling(latex_options = c( "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("SupTab2", keep_pdf = TRUE, file_format = "jpeg")

# kable(all_stats_table , format = "latex",  escape = F, 
#     booktabs = TRUE, align = align_tab1, digits = 3, linesep = "") %>% 
#     kable_styling(latex_options = c( "scale_down")) %>% 
#     row_spec(0, bold = TRUE) %>% 
#     add_header_above(c(" " = 2, "Conservation, demography, ecology and life-history data" = 6,
#         "Genetic data" = 5), italic = TRUE) %>% 
#     kable_as_image("SupTab1_test", keep_pdf = TRUE)
# supplementary table 3: Bottleneck stats

all_stats_table_sub3 <- all_stats %>% 
    left_join(all_stats_origin, by = c( "common",  "latin")) %>% 
    dplyr::select(91, 3, 78, 80, 82, 84, 89, 90) %>% 
    dplyr::rename(`$\\textbf{Common name}` = common_abbr,
        `$\\textbf{Scientific name}` = latin,
        `$\\textbf{TPM 70}` =  TPM70_ratio,
        `$\\textbf{TPM 80}` = TPM80_ratio,
        `$\\textbf{TPM 90}` = TPM90_ratio,
        `$\\textbf{SMM}$` = SMM_ratio,
        #`$p_{bot}$` = bot,
        #`$p_{neut}$` = neut) %>% 
         `\\boldsymbol{p_{bot}}` = bot,
         `\\boldsymbol{p_{neut}}` = neut) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(`$\\textbf{Scientific name}` = cell_spec(`$\\textbf{Scientific name}`, italic = TRUE))

align_tab3 <- c("l", "l", rep(x = "c", ncol(all_stats_table_sub3 ) - 2))

kable(all_stats_table_sub3, format = "latex", escape = F, 
    booktabs = TRUE, align = align_tab3, digits = 3, linesep = "") %>% 
    kable_styling(latex_options =  "scale_down") %>% 
   # row_spec(0, bold = TRUE) %>% 
    add_header_above(c(" "," ", "Heterozygosity-excess (prop_{het-exc})" = 4, "ABC" = 2), italic = F, escape = F) %>%
    kable_as_image("SupTab3", keep_pdf = TRUE)


# Supplementary table 3: Model fit ABC
model_fit <- read_delim("output/model_evaluation/check3_modeval/sims_10000k_p_vals_fit_30.txt", 
                        delim = " ") 

all_stats_table_sub4_temp <- all_stats %>% 
    left_join(all_stats_origin, by = c( "common",  "latin", "species")) %>% 
    full_join(model_fit, by = c("species")) 
 
all_stats_table_sub4 <-all_stats_table_sub4_temp %>% 
    dplyr::select("common_abbr", "latin", "bot") %>% 
    dplyr::mutate(model = ifelse(bot > 0.5, "bot", "neut")) %>% 
    dplyr::mutate(p_val = ifelse(bot > 0.5, all_stats_table_sub4_temp$bot_p, 
                                             all_stats_table_sub4_temp$neut_p)) %>% 
    dplyr::mutate(model = ifelse(model == "neut", "neutral", "bottleneck")) %>% 
    dplyr::select(-bot) %>% 
    dplyr::rename(`Common name` = common_abbr,
        `Scientific name` = latin,
        `Selected model` =  model,
         `\\textit{p} value` = p_val) %>% 
    dplyr::arrange(-row_number()) %>% 
    mutate(`Scientific name` = cell_spec(`Scientific name`, italic = TRUE))

align_tab4 <- c("l", "l", rep(x = "c", ncol(all_stats_table_sub4 ) - 2))

kable(all_stats_table_sub4, format = "latex", escape = F, 
    booktabs = TRUE, align = align_tab4 , digits = 3, linesep = "") %>% 
    kable_styling(latex_options =  "scale_down") %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("SupTab5", keep_pdf = TRUE)



################## Supplementary table 4 for ABC estimates ##############

## bottleneck model
library(tidyr)
library(dplyr)
library(MCMCglmm)
# load abc posterior data
load("data/processed/abc_estimates/abc_10000k_bot_complete_30.RData")
# unnest the abc posterior values for plotting
abc_bot <- unnest(abc_complete)

# load abc posterior data
load("data/processed/abc_estimates/abc_10000k_neut_complete.RData")
# unnest the abc posterior values for plotting
abc_neut <- unnest(abc_complete)

estimate_mode <- function(s) {
    d <- density(s, adjust = 2.5)
    d$x[which.max(d$y)]
}

names(abc_bot)
options(scipen=0)
abc_table_nbot <- abc_bot %>% 
    filter(pars %in% c("nbot")) %>% #"mut_rate"
    group_by(species, pars) %>% 
    summarise(mean = mean(adj_vals),
              median = median(adj_vals),
              mode = estimate_mode(adj_vals),
              `HPD lower`= HPDinterval(mcmc(adj_vals), 0.95)[1],
              `HPD upper` = HPDinterval(mcmc(adj_vals), 0.95)[2]) %>% 
             # HPD_low_50 = HPDinterval(mcmc(adj_vals), 0.50)[1],
             # HPD_high_50 = HPDinterval(mcmc(adj_vals), 0.50)[2]) %>% 
    left_join(all_stats_origin[c(2,4)]) %>% 
    #mutate_if(is.numeric, funs(formatC(., format = "e", digits = 2)))
    mutate_if(is.numeric, funs(round(., 0))) %>% 
    dplyr::select(common, everything()) %>% 
    dplyr::select(-species, -pars) %>% 
    mutate_if(is.numeric, funs(as.character(.)))


abc_table_mut <- abc_bot %>% 
    filter(pars %in% c("mut_rate")) %>% #"mut_rate"
    group_by(species, pars) %>% 
    summarise(mean = mean(adj_vals),
        median = median(adj_vals),
        mode = estimate_mode(adj_vals),
        `HPD lower`= HPDinterval(mcmc(adj_vals), 0.95)[1],
        `HPD upper` = HPDinterval(mcmc(adj_vals), 0.95)[2]) %>% 
    # HPD_low_50 = HPDinterval(mcmc(adj_vals), 0.50)[1],
    # HPD_high_50 = HPDinterval(mcmc(adj_vals), 0.50)[2]) %>% 
    left_join(all_stats_origin[c(2,4)]) %>% 
    mutate_if(is.numeric, funs(formatC(., format = "e", digits = 2))) %>% 
    #mutate_if(is.numeric, funs(round(., 0))) %>% 
    dplyr::select(common, everything()) %>% 
    dplyr::select(-species, -pars)
abc_table_mut <- abc_table_mut[-c(1,2)]

abc_bot_full <- bind_cols(abc_table_nbot, abc_table_mut) %>% 
    left_join(all_stats[c("species", "common_abbr")]) %>% 
    dplyr::select(species, common_abbr, mean:`HPD upper1`)
abc_bot_full <- abc_bot_full[-1]
names(abc_bot_full)[1] <- "Common name"
align_tab5 <- c("l", "l", rep(x = "c", ncol(abc_bot_full ) - 2))
#$N_{\\mathrm{e}}Bot$

kable(abc_bot_full,  format = "latex", escape = FALSE, 
    booktabs = TRUE, align = align_tab5, digits = 3, linesep = "",
    col.names = c("Common name", rep(names(abc_table_mut), 2))) %>% #c("Common name", rep(names(abc_table_mut), 2)))
    #kable_styling(latex_options =  "scale_down") %>% 
    add_header_above(c(" " = 1, "N_{e}bot" = 5, "Mu" = 5), escape = FALSE, bold = TRUE) %>%  #bold = T
    #add_header_above(c(" " = 1, "N_{e}Bot" = 5, "Mu" = 5), escape = FALSE, bold = TRUE) %>% 
    add_header_above(c(" " = 1, "Bottleneck model" = 10), bold = TRUE) %>% 
    row_spec(0, bold = TRUE) %>% 
    # group_rows("N_{bot}", 1,10, escape = F) %>% 
    # group_rows("Mutation rate mu", 11,20, escape = F) %>% 
    kable_as_image("Sup_tab_4a", file_format = "jpeg", keep_pdf = TRUE, density = 600)

# prediction errors:
# bottleneck: Nebot = 0.56, mu = 0.74
# neutral: mu = 0.68, GSM = 0.84

## Neutral model

abc_table_GSM_neut <- abc_neut %>% 
    filter(pars %in% c("gsm_param")) %>% #"mut_rate"
    group_by(species, pars) %>% 
    summarise(mean = mean(adj_vals),
        median = median(adj_vals),
        mode = estimate_mode(adj_vals),
        `HPD lower`= HPDinterval(mcmc(adj_vals), 0.95)[1],
        `HPD upper` = HPDinterval(mcmc(adj_vals), 0.95)[2]) %>% 
    # HPD_low_50 = HPDinterval(mcmc(adj_vals), 0.50)[1],
    # HPD_high_50 = HPDinterval(mcmc(adj_vals), 0.50)[2]) %>% 
    left_join(all_stats_origin[c(2,4)]) %>% 
    #mutate_if(is.numeric, funs(formatC(., format = "e", digits = 2)))
    mutate_if(is.numeric, funs(round(., 3))) %>% 
    #mutate_if(is.numeric, funs(formatC(., format = "e", digits = 2))) %>% 
    select(common, everything()) %>% 
    select(-species, -pars) %>% 
    mutate_if(is.numeric, funs(as.character(.)))


abc_table_mut_neut <- abc_neut %>% 
    filter(pars %in% c("mut_rate")) %>% #"mut_rate"
    group_by(species, pars) %>% 
    summarise(mean = mean(adj_vals),
        median = median(adj_vals),
        mode = estimate_mode(adj_vals),
        `HPD lower`= HPDinterval(mcmc(adj_vals), 0.95)[1],
        `HPD upper` = HPDinterval(mcmc(adj_vals), 0.95)[2]) %>% 
    # HPD_low_50 = HPDinterval(mcmc(adj_vals), 0.50)[1],
    # HPD_high_50 = HPDinterval(mcmc(adj_vals), 0.50)[2]) %>% 
    left_join(all_stats_origin[c(2,4)]) %>% 
    mutate_if(is.numeric, funs(formatC(., format = "e", digits = 2))) %>% 
    #mutate_if(is.numeric, funs(round(., 0))) %>% 
    select(common, everything()) %>% 
    select(-species, -pars)
abc_table_mut_neut <- abc_table_mut_neut[-c(1,2)]


abc_neut_full <- bind_cols(abc_table_GSM_neut, abc_table_mut_neut) %>% 
    left_join(all_stats[c("species", "common_abbr")]) %>% 
    select(species, common_abbr, mean:`HPD upper1`)
abc_neut_full <- abc_neut_full[-1]
align_tab7 <- c("l", "l", rep(x = "c", ncol(abc_neut_full) - 2))

kable(abc_neut_full,  format = "latex", escape = F, 
    booktabs = TRUE, align = align_tab5, digits = 3, linesep = "",
    col.names = c("Common name", rep(names(abc_table_mut_neut), 2))) %>% 
    kable_styling(latex_options =  "scale_down") %>% 
    add_header_above(c(" " = 1, "GSM_{par}" = 5, "mu" = 5), escape = F, bold = TRUE) %>% 
    add_header_above(c(" " = 1, "Neutral model" = 10), bold = TRUE) %>% 
    row_spec(0, bold = TRUE) %>% 
    # group_rows("N_{bot}", 1,10, escape = F) %>% 
    # group_rows("Mutation rate mu", 11,20, escape = F) %>% 
    kable_as_image("Sup_tab_4b", file_format = "jpeg", keep_pdf = TRUE, density = 600)


#### correlation matrix of the five predictor variables 
library(purrr)
pred_vars <- all_stats %>% 
    dplyr::select(logAbundance, SSD, TPM80_ratio, bot,  BreedingType)
all_combs <- expand.grid(names(pred_vars), names(pred_vars)) %>% 
                filter(Var1 != "BreedingType") %>% 
                filter(Var1 != Var2)

calc_r2 <- function(x,y){
    if(!(x==y)){
    mod <- lm(pred_vars[[x]]~pred_vars[[y]])
    summary(mod)$r.squared
    } else {
        1
    }
}
r2s <- unlist(map2(all_combs[[1]], all_combs[[2]], calc_r2))
# add to df
all_combs$r2s <- r2s
all_combs[1:2] <- lapply(all_combs[1:2], as.character)
# create matrix:
cor_mat <- data.frame(matrix(nrow = 5, ncol = 5), row.names = names(pred_vars))
names(cor_mat) <- names(pred_vars)
# fill in
for (x in names(cor_mat)){
    for (y in rownames(cor_mat)){
        if (!(x==y) & !(y=="BreedingType")){
            whichrow <- which((all_combs$Var1 == y) & (all_combs$Var2 == x))
            cor_mat[y, x] <- all_combs[whichrow, "r2s"]
        }
    }
}
#fill
cor_mat["BreedingType", 1:4] <- cor_mat[1:4, "BreedingType"]
#cor_mat[is.na(cor_mat)] <- 1
cor_mat[lower.tri(cor_mat)] <- NA
options(knitr.kable.NA = '')

cor_mat <- cor_mat %>% 
    rename("prop_{het-exc}" = "TPM80_ratio",
        "p_{bot}" = "bot",
        "Abundance" = "logAbundance",
        "Breeding Habitat" = "BreedingType")

rownames(cor_mat) <- c("Abundance", "SSD", "prop_{het-exc}",  "p_{bot}",
                       "Breeding Habitat")
str(cor_mat)

options(digits=2)
library(formattable)
cor_mat <- as.data.frame(cor_mat)
cor_mat[1] <- formattable(x = cor_mat[1], digits = 3, format = "f")
kable(cor_mat,  format = "latex", escape = F, 
    booktabs = TRUE, digits = 3, linesep = "") %>% 
   # kable_styling(latex_options =  "scale_down") %>% 
    #row_spec(0, bold = TRUE) %>% 
    # group_rows("N_{bot}", 1,10, escape = F) %>% 
    # group_rows("Mutation rate mu", 11,20, escape = F) %>% 
    kable_as_image("SupTab6", file_format = "jpeg", keep_pdf = TRUE, density = 600)


#### model output table

# bot ~ lh
mod_beta <- read_delim("output/mcmcmodels/hetexc_vs_lh_SSD_beta.txt", delim = " ")
mod_R2 <- read_delim("output/mcmcmodels/hetexc_vs_lh_SSD_R2.txt", delim = " ")
mod_SC <- read_delim("output/mcmcmodels/hetexc_vs_lh_SSD_SC.txt", delim = " ")
mod_beta2 <- read_delim("output/mcmcmodels/bot_vs_lh_SSD_beta.txt", delim = " ")
mod_R22 <- read_delim("output/mcmcmodels/bot_vs_lh_SSD_R2.txt", delim = " ")
mod_SC2 <- read_delim("output/mcmcmodels/bot_vs_lh_SSD_SC.txt", delim = " ")

options(digits=2)
betas <- rbind(mod_beta[2:3, ], mod_beta2[2:3, ])[c("post_median", "lower", "upper")]
betas <- paste0(round(betas[[1]], 2), " (", round(betas[[2]], 2),", ", round(betas[[3]], 2), ")")
r2s <- rbind(mod_R2[2:3, 2:4], mod_R22[2:3, 2:4])
r2s <- paste0(round(r2s[[1]], 2), " (", round(r2s[[2]], 2),", ", round(r2s[[3]], 2), ")")
scs <- rbind(mod_SC[3:5], mod_SC2[3:5])
scs  <- paste0(round(scs[[1]], 2), " (", round(scs[[2]], 2),", ", round(scs[[3]], 2), ")")
mod1_df <- data.frame(variables = rep(mod_SC$pred, 2), betas, r2s, scs)

# paste("Plot of ", alpha^beta, " versus ", hat(mu)[0]))
kable(mod1_df, format = "latex", booktabs = T, align = "c",
    col.names = c("Model", "$\\beta$", expression(R^2), "$r(\\hat{Y},x)$"), escape = F) %>%
    group_rows("prop_{het-exc}", 1, 2, latex_gap_space = "1em", escape = F) %>% 
    group_rows("p$_{bot}$", 3, 4, latex_gap_space = "1em",escape = F) %>% 
    kable_as_image("SupTab8", file_format = "jpeg", keep_pdf = TRUE, density = 600)


# Ar
# $prop_{\\mathrm{het-exc}}$
mod_beta_ar <- read_delim("output/mcmcmodels/gendiv_vs_lh_plus_bot_beta.txt", delim = " ")
mod_R2_ar <- read_delim("output/mcmcmodels/gendiv_vs_lh_plus_bot_R2.txt", delim = " ")
mod_SC_ar <- read_delim("output/mcmcmodels/gendiv_vs_lh_plus_bot_SC.txt", delim = " ")

betas_ar <- mod_beta_ar[2:6, c("post_median", "lower", "upper")]
betas_ar <- paste0(round(betas_ar[[1]], 2), " (", round(betas_ar[[2]], 2),", ", round(betas_ar[[3]], 2), ")")
r2s_ar <- mod_R2_ar[2:6, 2:4]
r2s_ar <- paste0(round(r2s_ar[[1]], 2), " (", round(r2s_ar[[2]], 2),", ", round(r2s_ar[[3]], 2), ")")
scs_ar <- mod_SC_ar[3:5]
scs_ar   <- paste0(round(scs_ar[[1]], 2), " (", round(scs_ar[[2]], 2),", ", round(scs_ar[[3]], 2), ")")
mod2_df <- data.frame(variables = mod_SC_ar$pred, betas_ar, r2s_ar, scs_ar)
mod2_df$variables <- c("SSD", "Breeding Habitat", "Abundance", "p$_{bot}$", "prop$_{het-exc}$")

kable(mod2_df, format = "latex", booktabs = T, align = "c",
    col.names = c("Model", "$\\beta$", expression(R^2), "$r(\\hat{Y},x)$"), escape = F) %>%
    group_rows("A_{r}", 1, 2, latex_gap_space = "1em", escape = F) %>% 
    kable_as_image("SupTab9", file_format = "jpeg", keep_pdf = TRUE, density = 600)
names(mod2_df) <- names(mod1_df)

# conservation
gendiv_beta <- read_delim("output/mcmcmodels/IUCN_gendiv_beta.txt", delim = " ")
gendiv_R2 <- read_delim("output/mcmcmodels/IUCN_gendiv_R2.txt", delim = " ")
hetexc_beta <- read_delim("output/mcmcmodels/IUCN_hetexc_beta.txt", delim = " ")
hetexc_R2 <- read_delim("output/mcmcmodels/IUCN_hetexc_R2.txt", delim = " ")
bot_beta <- read_delim("output/mcmcmodels/IUCN_bot_beta.txt", delim = " ")
bot_R2 <- read_delim("output/mcmcmodels/IUCN_bot_R2.txt", delim = " ")


iucn_mod <- data.frame(variables = rep("IUCN status", 3),
                        betas = c(
                            paste0(round(gendiv_beta[2, 3], 2), " (", round(gendiv_beta[2,5], 2), ", ", round(gendiv_beta[2,6], 2), ")" ),
                            paste0(round(hetexc_beta[2, 3], 2), " (", round(hetexc_beta[2,5], 2), ", ", round(hetexc_beta[2,6], 2), ")" ),
                            paste0(round(bot_beta[2, 3], 2), " (", round(bot_beta[2,5], 2), ", ", round(bot_beta[2,6], 2), ")" )),
                        r2s = rep(NA, 3),
                        scs = rep(NA, 3))
# bind all together
mod_full <- bind_rows(mod1_df, mod2_df, iucn_mod)

options(knitr.kable.NA = "")

mod_full$r2marginal <- c(
                        paste0(round(mod_R2[1,2], 2), " (", round(mod_R2[1,3], 2), ", ", round(mod_R2[1,4], 2), ")"),
                        NA,
                        paste0(round(mod_R22[1,2], 2), " (", round(mod_R22[1,3], 2), ", ", round(mod_R22[1,4], 2), ")"),
                        NA,
                        paste0(round(mod_R2_ar[1,2], 2), " (", round(mod_R2_ar[1,3], 2), ", ", round(mod_R2_ar[1,4], 2), ")"),
                        rep(NA,4),
                        paste0(round(gendiv_R2[1,2], 2), " (", round(gendiv_R2[1,3], 2), ", ", round(gendiv_R2[1,4], 2), ")"),
                        paste0(round(hetexc_R2[1,2], 2), " (", round(hetexc_R2[1,3], 2), ", ", round(hetexc_R2[1,4], 2), ")"),
                        paste0(round(bot_R2[1,2], 2), " (", round(bot_R2[1,3], 2), ", ", round(bot_R2[1,4], 2), ")")
    )

# expression(R^2~unique)
kable(mod_full, format = "latex", booktabs = T, align = "c",
    col.names = c("Model", "stand. $\\beta$", "$R^2_{unique}$", "$r(\\hat{Y},x)$", "$R^2_{marginal}$"), escape = F) %>%
    group_rows("prop_{het-exc}", 1, 2, latex_gap_space = "1em", escape = F) %>% 
    group_rows("p$_{bot}$", 3, 4, latex_gap_space = "1em",escape = F) %>% 
    group_rows("A_{r}", 5, 9, latex_gap_space = "1em", escape = F) %>% 
    group_rows("A_{r}", 10,10, latex_gap_space = "1em", escape = F) %>% 
    group_rows("prop_{het-exc}", 11,11, latex_gap_space = "1em", escape = F) %>% 
    group_rows("p$_{bot}$", 12,12, latex_gap_space = "1em", escape = F) %>% 
    kable_as_image("SupTab8_full", file_format = "jpeg", keep_pdf = TRUE, density = 600)
