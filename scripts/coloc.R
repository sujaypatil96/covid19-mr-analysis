library(tidyverse)
library(coloc)
setwd("~/Desktop/pm570_project/covid19-colocalization-analysis/")
#https://hanruizhang.github.io/GWAS-eQTL-Colocalization/



gwas <- read_table2(file = "../covid19-colocalization-analysis/data/hosp_vs_pop_GWAS")

eqtl <- read.table(file = "../covid19-colocalization-analysis/data/eQTL_spQTL_results/BLUEPRINT_eQTL/Monocyte.txt.gz",
                   header = T,
                   as.is = T)

eqtl <- eqtl %>% tibble()

gwas <- rename(gwas, P = all_inv_var_meta_p, CHR='#CHR')

# compute MAF
gwas$all_meta_MAF <- (1 - gwas$all_meta_AF)

# label snps
gwas["sid"] <- paste(gwas$CHR, gwas$POS, sep = ":")


input <- merge(x = eqtl,
               y = gwas,
               by = "sid",
               all = F,
               suffixes = c("_eqtl", "gwas"))




result <- coloc.abf(dataset1 = list(pvalues=input$P, type="cc", s=0.005, N=nrow(gwas)),
                    dataset2 = list(pvalues=input$npval, type="quant", N=nrow(eqtl)), MAF=input$all_meta_AF)

result <- coloc.abf(dataset1 = list(pvalues=input$P, type="cc", s=0.005, N=nrow(gwas)),
                    dataset2 = list(pvalues=input$npval, type="quant", N=nrow(eqtl)), MAF=input$all_meta_MAF)


result$summary
