# Load libraries and functions
source("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/scripts/IterateMatrixeQTLFunctions.R")
output_filename <- "pc_final_results.rds" # CHANGE ME

## Read in formatted data
genotype <- "tot.Breast_Carcinoma_genotype.csv"
snploc <- "tot.Breast_Carcinoma_snploc.csv"
geneloc <- "formatted_gene_loc.csv"
exp <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv" # CHANGE ME
### Read in data for HCPs
covariates <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formattedCovariates.csv")
Y <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv") # CHANGE ME

# Fix parameters
GENOTYPE_file_name = genotype; EXP_file_name = exp; gene_location_file_name = geneloc; snps_location_file_name =snploc
k=8; lambda1=2; lambda2=1; lambda3=1; first_two_plus_k_pcs=1

# Get eQTLs
eQTLs <- MeQTLandHcp(Y=Y, GENOTYPE_file_name=GENOTYPE_file_name, gene_location_file_name=gene_location_file_name, 
		     snps_location_file_name=snps_location_file_name, EXP_file_name=EXP_file_name,
		     covariates=covariates, first_two_plus_k_pcs=first_two_plus_k_pcs, lambda1=lambda1, 
		     lambda2=lambda2, lambda3=lambda3, k=k)

# Organize results
parameters <- data.frame("k"=k, 'lambda1'=lambda1, 'lambda2'=lambda2, 'lambda3'=lambda3, 'PCs'= 2+first_two_plus_k_pcs)

end_statistics <- list("parameters" = parameters,
		       	"eQTLs" = eQTLs)

saveRDS(end_statistics, output_filename) 
