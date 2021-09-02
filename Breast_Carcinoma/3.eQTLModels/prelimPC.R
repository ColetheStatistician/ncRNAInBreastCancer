# Load libraries and functions
source("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/scripts/IterateMatrixeQTLFunctions.R")

# Read in formatted data
genotype <- "tot.Breast_Carcinoma_genotype.csv"
snploc <- "tot.Breast_Carcinoma_snploc.csv"
geneloc <- "formatted_gene_loc.csv"
exp_pc <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv"

# Read in expression data
covariates <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formattedCovariates.csv")
Y_pc <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv")

# Prep data
formattedHCP <- formatHCPs(covariates = covariates, expression_matrix = Y_pc)
Z <- formattedHCP$Z
Y <- formattedHCP$Y

# Fix parameters
GENOTYPE_file_name = genotype; EXP_file_name = exp_pc; gene_location_file_name = geneloc; snps_location_file_name =snploc
k=10; lambda1=10; lambda2=2; lambda3=2; first_two_plus_k_pcs=2

# Get eQTLs
eQTLs_pc <- MeQTLandHcp(Z=Z, Y=Y, GENOTYPE_file_name=GENOTYPE_file_name, gene_location_file_name=gene_location_file_name, snps_location_file_name=snps_location_file_name, EXP_file_name=EXP_file_name,
			covariates=covariates, first_two_plus_k_pcs=first_two_plus_k_pcs, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, k=k)
                        
# Organize results
parameters <- data.frame("k"=10, 'lambda1'=10, 'lambda2'=2, 'lambda3'=2, 'first_two_plus_k_pcs'=2)

end_statistics <- list("parameters" = parameters,
		       "eQTLs_pc" = eQTLs_pc)

saveRDS(end_statistics, "pc_prelim_results.rds") 
