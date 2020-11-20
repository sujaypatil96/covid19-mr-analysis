setwd(dir = "C:\\Users/sujay/Desktop/USC Assignments and Material/PM570/Final Project/covid19-colocalization-analysis/scripts/")

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("snpStats")
# install("coloc")

# install.packages("devtools")
# devtools::install_github("boxiangliu/locuscomparer")

# install.packages("tidyverse")
# install.packages("metafor")
# install.packages("qqman")

library(tidyverse)
library(coloc)
library(locuscomparer)
library(metafor)
library(qqman)

# read in eQTL data
eqtl <- read.table(file = "../data/eQTL_spQTL_results/BLUEPRINT_eQTL/T-cell.txt",
                   header = T,
                   as.is = T)

head(eqtl)

# read in GWAS data, `nrows` argument added for circumventing memory constraints
gwas <- read.table(file = "../data/COVID19_HGI_B2_ALL_eur_leave_23andme_20201020.txt",
                   header = T,
                   as.is = T)

head(gwas)

# strict p-value filter used in the paper = 5*10^(-7)
sorted_filtered_gwas <- gwas %>% filter(all_inv_var_meta_p < 5*10^(-2)) %>% arrange(all_inv_var_meta_p)

# function to create 1Mbp locus/window for colocalization analysis
filter_snp <- function(gwas) {
    snps_passed = data.frame() %>% tibble()
    snps_filtered_empty = FALSE
    
    while(snps_filtered_empty != TRUE) {
        snps_filtered <- data.frame()
        
        snp_to_test <- gwas[1,]
        snps_left <- gwas[-1,]
        snps_filtered <- snps_left[abs(snp_to_test$POS - snps_left$POS) > 5e5,]
        
        snps_filtered_empty <- all(is.na(snps_filtered))
        print(snps_filtered_empty)
        
        gwas <- snps_filtered
        snps_passed <- rbind(snps_passed, snp_to_test)
    }
    
    snps_passed <- snps_passed %>% arrange(POS) %>% mutate(BP_LAG = POS - lag(POS))
    return(snps_passed)
}

# apply filter_snp() function to each chromosome
filtered_gwas_list <- lapply(sorted_filtered_gwas$CHR %>% unique() %>% sort(),
                             FUN = function(i) { filter_snp(sorted_filtered_gwas %>% filter(CHR==i)) })

filtered_gwas <- filtered_gwas_list %>% bind_rows()

# create column in "chromsome:position" format to match eQTL and GWAS datasets on
gwas["sid"] <- paste(gwas$CHR, gwas$POS, sep = ":")

# add MAF column to gwas dataset
gwas$all_meta_MAF <- (1 - gwas$all_meta_AF)

# match eQTL and GWAS datasets
input <- merge(x = eqtl,
               y = gwas,
               by = "sid",
               all = F,
               suffixes = c("_eqtl", "_gwas"))

head(input)

MR_df <- input %>% mutate(MR.est = all_inv_var_meta_beta / slope,
                          MR.se = all_inv_var_meta_sebeta / slope,
                          MR.z = MR.est / MR.se,
                          MR.p = 2 * pnorm(-abs(MR.z)))

forest(
    MR_df$MR.est,
    sei = MR_df$MR.se,
    slab = sprintf("    %s", MR_df$pid),
    xlab = "Causal Effect",
    annotate = FALSE,
    ilab = data.frame(
        sprintf("%.2f", exp(MR_df$MR.est)),
        sprintf("(%.2f, %.2f)", exp(MR_df$MR.est - MR_df$MR.z * MR_df$MR.se), exp(MR_df$MR.est + MR_df$MR.z * MR_df$MR.se)),
        sprintf("%.2e", MR_df$MR.p)),
    ilab.xpos = c(3, 3.5, 4.25),
    pch = 16,
    atransf = exp,
    at = log(c(0.5, 1, 2, 4, 8)),
    xlim = c(-2, 4.5)
)

data_manhattan <- data.frame(SNP=MR_df$pid,
                             CHR=MR_df$CHR,
                             BP=MR_df$POS,
                             P=MR_df$MR.p,
                             zscore=MR_df$MR.z)

manhattan(data_manhattan, main="Manhattan Plot - T-cells", 
          col = c("blue4", "orange3"), 
          suggestiveline = F, genomewideline = -log10(0.05))
