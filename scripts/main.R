setwd(dir = "C:\\Users/sujay/Desktop/USC Assignments and Material/PM570/Final Project/covid19-colocalization-analysis/scripts/")

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("snpStats")
# install("coloc")

# install.packages("tidyverse")

library(tidyverse)
library(coloc)

# read in eQTL data
eqtl <- read.table(file = "../data/eQTL_spQTL_results/BLUEPRINT_eQTL/Monocyte.txt",
                   header = T,
                   as.is = T)

head(eqtl)

# read in GWAS data, `nrows` argument added for circumventing memory constraints
gwas <- read.table(file = "../data/COVID19_HGI_B2_ALL_eur_leave_23andme_20201020.txt",
                   header = T,
                   nrows = 1000000,
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

# match eQTL and GWAS datasets
input <- merge(x = eqtl,
               y = gwas,
               by = "sid",
               all = F,
               suffixes = c("_eqtl", "_gwas"))

head(input)

# proportion = no. of cases/no. of controls
s_gwas <- 6404/902088

result <- coloc.abf(dataset1 = list(pvalues=input$all_inv_var_meta_p, type="cc", s=s_gwas, N=input$all_meta_sample_N),
                    dataset2 = list(pvalues=input$npval, type="quant", N=nrow(eqtl)), MAF=input$all_meta_AF)

result$summary
