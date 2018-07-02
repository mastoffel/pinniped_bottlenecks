# search terms for microsat data on web of science

library(readxl)
library(dplyr)
library(stringr)
all_species <- read_xlsx("data/raw/search_terms.xlsx",col_names = FALSE) %>% 
                   # mutate(search_term = paste0('("', X__1, '" OR "', X__2, '") AND microsat*'))
             mutate(search_term = paste0('("', X__1, '" OR "', X__2, '" OR "', X__3, '" OR "', X__4, '") AND microsat*')) 
            
all_species$search_term <- str_remove_all(all_species$search_term, paste0('" OR "', NA))

WriteXLS::WriteXLS(all_species,ExcelFileName = "data/raw/search_terms_formatted.xls", col.names = FALSE)

write_excel_csv(all_species, path = "data/raw/search_terms_formatted.csv", col_names = FALSE)
