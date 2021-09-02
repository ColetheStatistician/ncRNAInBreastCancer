source("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/scripts/MediAnalysisFunctions.R")
output_file <- "multiple_nc_MAresults.csv" # CHANGE ME

# Test Data
pc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_Covariates.csv")[,-1]
nc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_Covariates.csv")[,-1]
pc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv"))
nc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv"))
genotypes <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/tot.Breast_Carcinoma_genotype.csv")) # CHANGE ME
# Read in triplets
triplets <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/multiple_nc_RNAs.csv") # CHANGE ME
# Group by SNP and pc-gene (excluding repeats)
snps_and_pcs <- triplets %>% dplyr::group_by(snps, pc_trans_genes) %>% dplyr::arrange(nc_cis_pvalue) %>% dplyr::slice(1) %>% dplyr::select(snps, pc_trans_genes)
# snps_and_pcs <- snps_and_pcs[1:10,] # DELETE ME

# Mediation Analysis
MAresults <- list()
for(i in 1:dim(snps_and_pcs)[1]) {
  MAresults[[i]] <- Mediation(desired_SNP=snps_and_pcs$snps[i], 
                             genotypes=genotypes, 
                             expressionMatrix = pc_exp, 
                             mediators_matrix = nc_exp,
                             pc_covariates=pc_covariates, 
                             nc_covariates=nc_covariates, 
                             pc_desired_ENSEMBL=snps_and_pcs$pc_trans_genes[i], 
                             triplets_df=triplets, 
                             total_scenarios=dim(snps_and_pcs)[1], 
                             scenario=i)
  # Get nc RNAs
  ncRNAs <- triplets %>% filter(snps == snps_and_pcs$snps[i], pc_trans_genes == snps_and_pcs$pc_trans_genes[i]) %>% dplyr::arrange(nc_cis_pvalue) %>% select(nc_cis_genes) # TODO: Select column containing corresponding ncRNA IDs
  # Combine results with SNP, pc-gene, and nc-genes information in list
  MAresults[[i]] <- list(
    "snp" = snps_and_pcs$snps[i],
    "pcgene" = snps_and_pcs$pc_trans_genes[i],
    "ncRNAs" = ncRNAs,
    "results" = MAresults[[i]]
  )
}

# write results
write.csv(MAresults, output_file, row.names=FALSE)

