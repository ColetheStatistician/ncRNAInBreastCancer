# Load libraries and functions
source("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/scripts/IterateMatrixeQTLFunctions.R")
output_filename <- "nc_iterate_PC_results.rds" # CHANGE ME

## Read in formatted data
genotype <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/tot.Breast_Carcinoma_genotype.csv"
snploc <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/tot.Breast_Carcinoma_snploc.csv"
geneloc <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_gene_loc.csv"
exp <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv" # CHANGE ME
### Read in data for HCPs
covariates <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formattedCovariates.csv")
Y <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv") # CHANGE ME

## Set variables
GENOTYPE_file_name = genotype; EXP_file_name = exp; gene_location_file_name = geneloc; snps_location_file_name =snploc
first_two_plus_k_pcs=6; k=1; lambda1=1; lambda2=1; lambda3=1; covariates=covariates # k and lambdas obtained from other script
toOpt = c("first_two_plus_k_pcs")

# Iterate
eQTLs <- iterateMeQTL(first_two_plus_k_pcs=first_two_plus_k_pcs, # how many of first principal components to use ( after the first 5)
                      k=k, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, # hcp variables to optimize
                      Y=Y,
                      covariates=covariates, # matrix of covariates
                      GENOTYPE_file_name=GENOTYPE_file_name, 
                      EXP_file_name=EXP_file_name, 
                      gene_location_file_name=gene_location_file_name, 
                      snps_location_file_name=snps_location_file_name,
                      toOpt = toOpt)

# Find which scenario is optimized
optResults <- getOptimized(iterateMeQTL_results = eQTLs,
                           first_two_plus_k_pcs=first_two_plus_k_pcs, # how many of first principal components to use
                           k=k, lambda1=lambda1, lambda2=lambda2, lambda3=lambda3, # hcp variables to optimize
                           toOpt = toOpt)
  
saveRDS(optResults, output_filename)
