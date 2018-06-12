# parses output of bottleneck
# out1: TPM 70, IAM and SMM
# out2: TPM 95
# out3: TPM 99
# out24: TPM 90

# run these bash commands on all files in bottleneck/out to remove header information
# make sure it works

system("for i in output/bottleneck_out/*.txt; do cat $i | awk '/SIGN/,0' > output/bottleneck_out/$i.edit.txt; done")

system("for i in output/bottleneck_out/hw/*.txt; do cat $i | awk '/SIGN/,0' > output/bottleneck_out/hw/$i.edit.txt; done")

# for i in *.txt; do cat $i | awk '/SIGN/,0' > $i.edit.txt; done

library(stringr)
library(dplyr)
library(plyr)

# paths to full datasets
in_1 <- paste("output/bottleneck_out/", list.files(path = "output/bottleneck_out", pattern="*1.txt.edit.txt"), sep = "")
in_2 <- paste("output/bottleneck_out/", list.files(path = "output/bottleneck_out", pattern="*2.txt.edit.txt"), sep = "")
in_3 <- paste("output/bottleneck_out/", list.files(path = "output/bottleneck_out", pattern="*3.txt.edit.txt"), sep = "")
#in_4 <- paste("output/bottleneck_out/", list.files(path = "output/bottleneck_out", pattern="*4.txt.edit.txt"), sep = "")

# paths to datasets in hwe
in1_hw <- paste("output/bottleneck_out/hw/", list.files(path = "output/bottleneck_out/hw", pattern="*1.txt.edit.txt"), sep = "")
in2_hw <- paste("output/bottleneck_out/hw/", list.files(path = "output/bottleneck_out/hw", pattern="*2.txt.edit.txt"), sep = "")
in3_hw <- paste("output/bottleneck_out/hw/", list.files(path = "output/bottleneck_out/hw", pattern="*3.txt.edit.txt"), sep = "")
#in4_hw <- paste("output/bottleneck_out/hw/", list.files(path = "output/bottleneck_out/hw", pattern="*4.txt.edit.txt"), sep = "")

import_files <- function(file){
        file <- readLines(file)
}

bottleneck_out1 <- lapply(in_1, import_files)
bottleneck_out2 <- lapply(in_2, import_files)
bottleneck_out3 <- lapply(in_3, import_files)
#bottleneck_out4 <- lapply(in_4, import_files)

bottleneck_out1_hw <- lapply(in1_hw, import_files)
bottleneck_out2_hw <- lapply(in2_hw, import_files)
bottleneck_out3_hw <- lapply(in3_hw, import_files)
#bottleneck_out4_hw <- lapply(in4_hw, import_files)

# x <- bottleneck_out1#[[2]]
# y <- bottleneck_out2#[[2]]
# z <- bottleneck_out3#[[2]]

#------------------------------------------------------------------------------
# function to get the main values from bottleneck output file 1, TPM 70, IAM, SMM

get_vals <- function(x) {
# remove caution messages if present
if (grepl("Caution", x[19]) == TRUE) {
        run_1 <- x[-c(19,20)]} else {
        run_1 <- x
        }
        
# get vals from first file
iam <- data.frame(run_1[c(3,4,5,20,31,32,33)])
colnames(iam) <- "out"
tpm70 <- data.frame(run_1[c(8,9,10,23,36,37,38)])
colnames(tpm70) <- "out"
smm <- data.frame(run_1[c(13,14,15,26,41,42,43)])
colnames(smm) <- "out"
modeshift <- data.frame(run_1[49])
colnames(modeshift) <- "out"
out <- rbind(iam, tpm70, smm, modeshift)
}

out1 <- lapply(bottleneck_out1, get_vals)
out1_hw <- lapply(bottleneck_out1_hw, get_vals)

#------------------------------------------------------------------------------
# function to get the main values from bottleneck output file 2, TPM 80

get_vals2 <- function(y) {
# remove caution messages if present
if (grepl("Caution", y[9]) == TRUE) {
        run_2 <- y[-c(9,10)]} else {
        run_2 <- y
        }
# get vals from second file
tpm80 <- data.frame(run_2[c(3,4,5,10,15,16,17)])
colnames(tpm80) <- "out"
tpm80 <- tpm80
}

out2 <- lapply(bottleneck_out2, get_vals2)
out2_hw <- lapply(bottleneck_out2_hw, get_vals2)

#------------------------------------------------------------------------------
# function to get the main values from bottleneck output file 3, TPM 90

get_vals3 <- function(y) {
        # remove caution messages if present
        if (grepl("Caution", y[9]) == TRUE) {
                run_3 <- y[-c(9,10)]} else {
                        run_3 <- y
                }
        # get vals from second file
        tpm90 <- data.frame(run_3[c(3,4,5,10,15,16,17)])
        colnames(tpm90) <- "out"
        tpm90 <- tpm90
}

out3 <- lapply(bottleneck_out3, get_vals3)
out3_hw <- lapply(bottleneck_out3_hw, get_vals3)

#------------------------------------------------------------------------------
# function to get the main values from bottleneck output file 4, TPM X

get_vals4 <- function(y) {
        # remove caution messages if present
        if (grepl("Caution", y[9]) == TRUE) {
                run_4 <- y[-c(9,10)]} else {
                        run_4 <- y
                }
        # get vals from second file
        tpm90 <- data.frame(run_4[c(3,4,5,10,15,16,17)])
        colnames(tpm90) <- "out"
        tpm90 <- tpm90
}

#out4 <- lapply(bottleneck_out4, get_vals4)
#out4_hw <- lapply(bottleneck_out4_hw, get_vals4)

#------------------------------------------------------------------------------
# rbind all vals together

data <- mapply(rbind, out1, out2, out3, SIMPLIFY=F)
data_hw <- mapply(rbind, out1_hw, out2_hw, out3_hw, SIMPLIFY=F)

# grep everything out that we don't need

get_stats <- function(x){
        library(dplyr)
x <- data.frame(x) %>%
        mutate(out = gsub("\\s.+:", "", out)) %>%
        mutate(out = gsub(":", "", out)) %>%
        mutate(out = gsub("Probability", "", out)) %>%
        mutate(out = gsub("Expected", "", out)) %>%
        mutate(out = gsub("T2", "", out)) %>%
        mutate(out = gsub("loci with heterozygosity deficiency", "", out)) %>%
        mutate(out = gsub("\\sloci with heterozygosity excess.", "", out)) %>%
        mutate(out = gsub("and", "vs", out))

out <- t(x)

}

out <- lapply(data, get_stats)
out_hw <- lapply(data_hw, get_stats)

#------------------------------------------------------------------------------
# make dataframe and add name elements

get_colnames <- function(x) {
        x <- data.frame(x)
        colnames(x) <- c("IAM_Heq", "IAM_Def/Exc", "IAM_Sign", "IAM_Stand_Diff", "IAM_Wilc_Def", "IAM_Wilc_Exc", "IAM_Wilc_2Tail", 
                         "TPM70_Heq", "TPM70_Def/Exc", "TPM70_Sign", "TPM70_Stand_Diff", "TPM70_Wilc_Def", "TPM70_Wilc_Exc", "TPM70_Wilc_2Tail",
                         "SMM_Heq", "SMM_Def/Exc", "SMM_Sign", "SMM_Stand_Diff", "SMM_Wilc_Def", "SMM_Wilc_Exc", "SMM_Wilc_2Tail",
                         "Mode_Shift",
                         "TPM80_Heq", "TPM80_Def/Exc", "TPM80_Sign", "TPM80_Stand_Diff", "TPM80_Wilc_Def", "TPM80_Wilc_Exc", "TPM80_Wilc_2Tail",
                         "TPM90_Heq", "TPM90_Def/Exc", "TPM90_Sign", "TPM90_Stand_Diff", "TPM90_Wilc_Def", "TPM90_Wilc_Exc", "TPM90_Wilc_2Tail")
        x <- x[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36)] # CHECK THIS
        #x <- x[c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,15,16,17,18,19,20,21,22)] # CHECK THIS

}


out <- lapply(out, get_colnames)

out_hw <- lapply(out_hw, get_colnames)

#------------------------------------------------------------------------------
# specify dataset names
library(plyr)
dataset_names <-list.files(path = "output/bottleneck_out", pattern="*2.txt.edit.txt")
dataset_names <- gsub("_genepop_out2.txt.edit.txt", "", dataset_names)
names(out) <- dataset_names
out <- ldply(out, data.frame)

dataset_names_hw <-list.files(path = "output/bottleneck_out/hw", pattern="*2.txt.edit.txt")
dataset_names_hw <- gsub("_genepop_out2.txt.edit.txt", "", dataset_names_hw)
names(out_hw) <- dataset_names_hw
out_hw <- ldply(out_hw, data.frame)

# write out files
library(WriteXLS)
WriteXLS(out, ExcelFileName = "output/out_bottleneck_stats_30.xls", na = "", row.names = F, col.names = T)
WriteXLS(out_hw, ExcelFileName = "output/out_bottleneck_stats_HW_30_Jun2018.xls", na = "", row.names = F, col.names = T)

write.csv(out, file = "output/out_bottleneck_stats_30.csv", na = "", row.names = F, col.names = T, quote = F)
write.csv(out_hw, file = "output/out_bottleneck_stats_HW_30_Jun2018.csv", na = "", row.names = F, col.names = T, quote = F)
