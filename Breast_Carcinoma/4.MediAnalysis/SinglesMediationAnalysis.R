library(magicfor)
source("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/scripts/MediAnalysisFunctions.R")
output_file <- "singleMAresults.csv" # CHANGE ME

# Test Data
pc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_Covariates.csv")[,-1]
nc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_Covariates.csv")[,-1]
pc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv"))
nc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv"))
genotypes <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/tot.Breast_Carcinoma_genotype.csv")) # CHANGE ME
# Read in triplets
triplets <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/singles.csv") # CHANGE ME: Change from singles to all
# Get subset of triplets
#triplets <- triplets[1:200,] # DELETE ME

# Mediation Analysis
magic_for(silent = TRUE)
MAresults <- for (i in 1:dim(triplets)[1]) {
  tempMAresults <- Mediation(desired_SNP=triplets$snps[i], 
                             genotypes=genotypes, 
                             expressionMatrix = pc_exp, 
                             mediators_matrix = nc_exp,
                             pc_covariates=pc_covariates, 
                             nc_covariates=nc_covariates, 
                             pc_desired_ENSEMBL=triplets$pc_trans_genes[i], 
                             triplets_df=triplets, 
                             total_scenarios=dim(triplets)[1], 
                             scenario=i)
  put(tempMAresults$test.stat, tempMAresults$p.value)
}
MAresults <- magic_result_as_dataframe()[,-1]

# Join results with triplet information
final <- cbind(triplets, MAresults)
colnames(final)[(dim(final)[2]-1):dim(final)[2]] <- c("MedA.test.stat", "MedA.p.value")

# write results
write.csv(final, output_file, row.names=FALSE)
