setwd(dir = "C:\\Users/sujay/Desktop/USC Assignments and Material/PM570/Final Project/covid19-colocalization-analysis/scripts/")

# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("snpStats")
# install("coloc")

library(coloc)

eqtl <- read.table(file = "../data/eQTL_spQTL_results/BLUEPRINT_eQTL/Monocyte.txt",
                   header = T,
                   as.is = T)

eqtl_snps <- eqtl[,"sid"]

write.table(x = eqtl_snps, 
            file = "../data/snps.txt",
            row.names = F,
            quote = F)

head(eqtl)

gwas <- read.table(file = "../data/COVID19_HGI_B2_ALL_eur_leave_23andme_20201020.txt",
                   header = T,
                   nrows = 10000000,
                   as.is = T)

gwas["sid"] <- paste(gwas$CHR, gwas$POS, sep = ":")

input <- merge(x = eqtl,
               y = gwas,
               by = "sid",
               all = F,
               suffixes = c("_eqtl", "gwas"))

head(input)

result <- coloc.abf(dataset1 = list(pvalues=input$all_inv_var_meta_p, type="cc", s=0.005, N=nrow(gwas)),
                    dataset2 = list(pvalues=input$npval, type="quant", N=nrow(eqtl)), MAF=input$all_meta_AF)

result$summary
