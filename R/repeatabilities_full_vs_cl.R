# check whether clustered and not clustered results are similar using repeatabilities

library(tidyr)
library(rptR)
library(dplyr)
library(readr)

seals <- read_csv("data/processed/all_stats_tree_30.csv")
seals_cl <- read_csv("data/processed/all_stats_tree_30_cl.csv")
all_seals <- rbind(seals, seals_cl)

all_stats_long <- all_seals %>% 
    dplyr::select(species, num_alleles_mean, obs_het_mean, TPM70_ratio, TPM80_ratio, TPM90_ratio, SMM_ratio, bot) 

calc_rpts <- function(var_in_all_stats){
    formula_temp <- as.formula(paste0(var_in_all_stats, "~ (1|species)"))
    rpt_temp <- rpt(formula = formula_temp, data = all_stats_long, datatype = "Gaussian", 
                    grname = "species",
                    nboot = 10, npermut = 0)
    out <- data.frame(variable = var_in_all_stats, R = rpt_temp$R, CI = rpt_temp$CI_emp, 
            row.names = NULL) %>% 
            mutate_if(is.numeric, funs(round(., 3))) %>% 
            rename(R = species, CI25 = CI.2.5., CI975 = CI.97.5.) %>% 
            mutate(CI = paste0("[",CI25," ", CI975, "]")) %>% 
            dplyr::select(variable, R, CI)
    out
}

# calculate all repeatabilities
all_rpts <- lapply(names(all_stats_long[2:ncol(all_stats_long)]), calc_rpts)
all_rpts_df <- do.call(rbind, all_rpts)
library(knitr)
library(kableExtra)

var_names <- c("A_r", "H_o", "TPM 70", "TPM 80", "TPM 90", "SMM", "p_{bot}")
all_rpts_df$variable <- var_names

kable(all_rpts_df, format = "latex",  escape = F, 
    booktabs = TRUE, align = "c", digits = 3, linesep = "") %>% 
    kable_styling(latex_options = c( "scale_down")) %>% 
    row_spec(0, bold = TRUE) %>% 
    kable_as_image("SupTab5", keep_pdf = TRUE)


