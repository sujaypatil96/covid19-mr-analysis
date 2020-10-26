library(tidyverse)
library(coloc)
setwd("~/Desktop/pm570_project/covid19-colocalization-analysis/")
#https://hanruizhang.github.io/GWAS-eQTL-Colocalization/



gwas <- read_table2(file = "../covid19-colocalization-analysis/data/hosp_vs_pop_23andMe_10K.gz")

eqtl <- read.table(file = "../covid19-colocalization-analysis/data/eQTL_spQTL_results/BLUEPRINT_eQTL/Monocyte.txt.gz",
                   header = T,
                   as.is = T)

eqtl <- eqtl %>% tibble()


# kallisto expression data
# psi: splicing
# the different data types are BLUEPRINT, DICE, DGN, GEUVADIS
#gwas <- gwas %>% mutate(P=all_inv_var_meta_p) 
gwas <- rename(gwas, P = all_inv_var_meta_p, CHR='#CHR')

sorted_filtered_gwas <- gwas %>% filter(P < .05) %>% arrange(P)


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


filter_snp(chr_gwas <- sorted_filtered_gwas %>% filter(CHR == 5))

filtered_gwas_list <- lapply(sorted_filtered_gwas$CHR %>% unique() %>% sort(),
                        FUN = function(i) {filter_snp(sorted_filtered_gwas 
                                                      %>% filter(CHR==i))})

filtered_gwas <- filtered_gwas_list %>% bind_rows()




filtered_gwas["sid"] <- paste(filtered_gwas$CHR, filtered_gwas$POS, sep = ":")

input <- merge(x = eqtl,
               y = filtered_gwas,
               by = "sid",
               all = F,
               suffixes = c("_eqtl", "gwas"))


# next step is colocalizing the eqtl




result <- coloc.abf(dataset1 = list(pvalues=input$P, type="cc", s=0.005, N=nrow(filtered_gwas)),
                    dataset2 = list(pvalues=input$npval, type="quant", N=nrow(eqtl)), MAF=input$all_meta_AF)
result$summary




