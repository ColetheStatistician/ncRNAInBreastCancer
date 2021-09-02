library(magicfor)
source("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/scripts/MediAnalysisFunctions.R")
# Test Data
pc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_Covariates.csv")[,-1]
nc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_Covariates.csv")[,-1]
pc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv"))
nc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv"))
genotypes <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/22.Breast_Carcinoma_genotype.csv")) # CHANGE ME
# SNPS to get from CHR22
optPC <- readRDS("pc_final_trans_results.rds") ## DELETE ME FOR TOTAL GENOME
optNC <- readRDS("nc_final_results.rds") ## DELETE ME FOR TOTAL GENOME
## Get list of unique SNPs in intersection between the pc cis-eQTLs and nc trans-eQTLs ## DELETE ME FOR TOTAL GENOME
uniq_snps <- intersect(optPC$eQTLs$trans$eqtls$snps, optNC$eQTLs$cis$eqtls$snps) ## DELETE ME FOR TOTAL GENOME
## Table of top p value trans qtls with overlapping cis ## DELETE ME FOR TOTAL GENOME
trans_w_overlapping_cis <- optPC$eQTLs$trans$eqtls[which(optPC$eQTLs$trans$eqtls$snps %in% uniq_snps),] ## DELETE ME FOR TOTAL GENOME
CHR22snps <- trans_w_overlapping_cis[which(trans_w_overlapping_cis$snps %in% genotypes$id),] ## DELETE ME FOR TOTAL GENOME
# Read in triplets
triplets <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/triplets.csv")
triplets <- triplets[which(triplets$snps %in% CHR22snps$snps),] ## DELETE ME FOR TOTAL GENOME

MAresults <- list()

# Mediation Analysis
magic_for(silent = TRUE)
MAresults <- for (i in 1:dim(triplets)[1]) {
  tempMAresults <- medAn(desired_SNP=triplets$snps[i], 
                         genotypes=genotypes, 
                         pc_covariates=pc_covariates, 
                         nc_covariates=nc_covariates, 
                         pc_desired_ENSEMBL=triplets$pc_trans_genes[i], 
                         nc_desired_ENSEMBL=triplets$nc_cis_genes[i], 
                         total_scenarios=dim(triplets)[1], 
                         scenario=i)
  put(tempMAresults$test.stat, tempMAresults$p.value)
}
MAresults <- magic_result_as_dataframe()[,-1]

# Join results with triplet information
final <- cbind(triplets, MAresults)
colnames(final)[(dim(final)[2]-1):dim(final)[2]] <- c("MedA.test.stat", "MedA.p.value")

# write results
write.csv(final, "MAsubsetresults.csv", row.names=FALSE)
