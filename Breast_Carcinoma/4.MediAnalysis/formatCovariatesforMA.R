# Load libraries and functions
library(Rhcpp)
source("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/scripts/IterateMatrixeQTLFunctions.R")

# Read in necessary data
covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formattedCovariates.csv")
raw_expr <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv"))
raw_mi <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv"))

# Set output file names
pc_output_filename <- "formatted_pc_Covariates.csv"
nc_output_filename <- "formatted_nc_Covariates.csv"

# incorporate HCPs
## Prep data for HCPs
pc_formattedHCP <- formatHCPs(covariates = covariates, expression_matrix = raw_expr)
nc_formattedHCP <- formatHCPs(covariates = covariates, expression_matrix = raw_mi)

pc_Z <- pc_formattedHCP$Z
nc_Z <- nc_formattedHCP$Z

pc_Y <- pc_formattedHCP$Y
nc_Y <- nc_formattedHCP$Y

# Subset to desired PCs
pc_chosenPCs <- make.names(rep("PC", 4+3), unique = TRUE)[-1] # the `+3` ensures there is at least 2 PCs as covariates
pc_variablesToGet <- c(c(covariates[,1])[(dim(covariates)[1] - 4):dim(covariates)[1]], pc_chosenPCs)
pc_covariates <- covariates[which(covariates$id %in% pc_variablesToGet),]

nc_chosenPCs <- make.names(rep("PC", 2+3), unique = TRUE)[-1] # the `+3` ensures there is at least 2 PCs as covariates
nc_variablesToGet <- c(c(covariates[,1])[(dim(covariates)[1] - 4):dim(covariates)[1]], nc_chosenPCs)
nc_covariates <- covariates[which(covariates$id %in% nc_variablesToGet),]

## Calculate HCPs
### Values obtained from optimized models
pc_HCPs <- getHCPs(first_two_plus_k_pcs = 4, k = 7, lambda1 = 1, lambda2 = 2, lambda3 = 1, Z = pc_Z, Y = pc_Y)
nc_HCPs <- getHCPs(first_two_plus_k_pcs = 2, k = 1, lambda1 = 1, lambda2 = 1, lambda3 = 1, Z = nc_Z, Y = nc_Y)

# Format HCPS and make covariates
pc_covariates <- formatHCPstoCov(hcps = pc_HCPs, covariates = pc_covariates)
nc_covariates <- formatHCPstoCov(hcps = nc_HCPs, covariates = nc_covariates)

# Exclude categorical variables
drop <- c("patient.breast_carcinoma_estrogen_receptor_status", "patient.gender", "patient.menopause_status", "patient.stage_event.pathologic_stage")
pc_covariates <- pc_covariates[which(!(pc_covariates$id %in% drop)),]
nc_covariates <- nc_covariates[which(!(nc_covariates$id %in% drop)),]

# Get final format
## Rows correspond to probes and columns to samples. All numeric
pc_covariates <- pc_covariates[,-1]
colnames(pc_covariates) <- NULL
rownames(pc_covariates) <- NULL
pc_covariates <- apply(pc_covariates, 2, as.numeric)

nc_covariates <- nc_covariates[,-1]
colnames(nc_covariates) <- NULL
rownames(nc_covariates) <- NULL
nc_covariates <- apply(nc_covariates, 2, as.numeric)

# Write files
write.csv(pc_covariates, pc_output_filename)
write.csv(nc_covariates, nc_output_filename)
