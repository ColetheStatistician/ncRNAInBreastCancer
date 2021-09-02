library(dplyr)
library(magrittr)
# Goal of script is to take the results from the optimized PC and NC eQTL models, find the SNP intersection, and output them in a conveniant format for mediation analysis
# Read in results from optimized models
optPC <- readRDS("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/pc_final_results.rds")
optNC <- readRDS("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/nc_final_results.rds")

# Get intersection
snp_int <- intersect(optPC$eQTLs$trans$eqtls$snps, optNC$eQTLs$cis$eqtls$snps)

# Format as triples
## Get genes corresponding to SNPs
trans_pc <- optPC$eQTLs$trans$eqtls[which(optPC$eQTLs$trans$eqtls$snps %in% snp_int),]
cis_nc <- optNC$eQTLs$cis$eqtls[which(optNC$eQTLs$cis$eqtls$snps %in% snp_int),]
## Match trans-pc and cis-ns by SNP
both <- dplyr::inner_join(trans_pc, cis_nc, by = "snps")
colnames(both) <- c("snps", "pc_trans_genes", "pc_trans_statistic", "pc_trans_pvalue", "pc_trans_FDR", "pc_trans_beta",
                    "nc_cis_genes", "nc_cis_statistic", "nc_cis_pvalue", "nc_cis_FDR", "nc_cis_beta")
# Arrange by significance
both <- both %>%
  arrange("pc_trans_pvalue", "pc_trans_FDR", "nc_cis_statistic", "nc_cis_FDR")

## SNP pc-gene pairs have multiple nc-RNAs
multiple <- both %>%
  group_by(snps, pc_trans_genes) %>%
  mutate(num = n() ) %>%
  arrange(desc(num)) %>%
  filter(num>1)
which(both$snps 

# Univariate mediation analysis situation
singles <- both %>%
  group_by(snps, pc_trans_genes) %>%
  mutate(num = n() ) %>%
  arrange(desc(num)) %>%
  filter(num == 1)

# Save files
write.csv(both, "triplets.csv", row.names = FALSE)
write.csv(multiple, "multiple_nc_RNAs.csv", row.names = FALSE)
write.csv(singles, "singles.csv", row.names = FALSE)
