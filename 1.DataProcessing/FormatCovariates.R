# Goal is to produce covariate file that matches this format: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt
cat("Formatting covariates\n")
# Load libraries
library(data.table)

# Load data
clinical <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/Clinical/Breast_Carcinoma_Covariates.txt")
# Set covariate ids to the side
temp <- clinical[,1]

#######################################
# Subset to European Descent patients #
#######################################
# Get ED participants
EDpatients <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/uniquePatients.csv")

# Remove non-European ancestry individuals
colstoget <- which(clinical[1358,] == "white")
clinical <- clinical[,..colstoget]

# get patients
p_clinical <- toupper(clinical[20,])

# Cut to first 3 TCGA identification sections
p_clinical_sub <- substr(p_clinical, 1, 12)

## subset to patients in all datasets
clinical <- cbind(temp, data.frame(clinical)[,which(p_clinical_sub %in% as.matrix(EDpatients))])

##############################
# Get covariates of interest #
##############################
# Covariates:
# my descriptor - data descriptor - row index
# age - patient.age_at_initial_pathologic_diagnosis - 12
# patient tcga barcode - patient.bcr_patient_barcode - 20
# estrogen - patient.breast_carcinoma_estrogen_receptor_status - 23
# gender - patient.gender - 1138
# menopause - patient.menopause_status - 1163
# race - patient.race_list.race - 1358
# pathologic stage - patient.stage_event.pathologic_stage - 1460

# Select these covariates
covariates <- clinical[c(12, 23, 1138, 1163, 1460),]
# Make patient barcodes the column names
colnames(covariates) <- c("id", toupper(as.character(clinical[20,-1])))
# Set the covariate ids to the side minus race
temp <- covariates[-5,1]

################################################
# Code values to be compatible with MatrixeQTL #
################################################
## estrogen
colstoget <- which(covariates[2,] == "negative")
covariates[2, colstoget] <- 0
colstoget <- which(covariates[2,] == "positive")
covariates[2, colstoget] <- 1
colstoget <- which(covariates[2,] == "indeterminate")
covariates[2, colstoget] <- 2
colstoget <- which(is.na(covariates[2,]))
covariates[2, colstoget] <- 3
## gender
colstoget <- which(covariates[3,] == "female")
covariates[3, colstoget] <- 0
colstoget <- which(covariates[3,] == "male")
covariates[3, colstoget] <- 1
## menopause
colstoget <- which(covariates[4,] == "pre (<6 months since lmp and no prior bilateral ovariectomy and not on estrogen replacement)")
covariates[4, colstoget] <- 0
colstoget <- which(covariates[4,] == "peri (6-12 months since last menstrual period)")
covariates[4, colstoget] <- 1
colstoget <- which(covariates[4,] == "post (prior bilateral ovariectomy or >12 mo since lmp with no prior hysterectomy)")
covariates[4, colstoget] <- 2
colstoget <- which(covariates[4,] == "indeterminate (neither pre or postmenopausal)")
covariates[4, colstoget] <- 3
colstoget <- which(is.na(covariates[4,]))
covariates[4, colstoget] <- 4
## pathologic stage
colstoget <- which(covariates[5,] == "stage i")
covariates[5, colstoget] <- 0
colstoget <- which(covariates[5,] == "stage ia")
covariates[5, colstoget] <- 1
colstoget <- which(covariates[5,] == "stage ib")
covariates[5, colstoget] <- 2
colstoget <- which(covariates[5,] == "stage ii")
covariates[5, colstoget] <- 3
colstoget <- which(covariates[5,] == "stage iia")
covariates[5, colstoget] <- 4
colstoget <- which(covariates[5,] == "stage iib")
covariates[5, colstoget] <- 5
colstoget <- which(covariates[5,] == "stage iii")
covariates[5, colstoget] <- 6
colstoget <- which(covariates[5,] == "stage iiia")
covariates[5, colstoget] <- 7
colstoget <- which(covariates[5,] == "stage iiib")
covariates[5, colstoget] <- 8
colstoget <- which(covariates[5,] == "stage iiic")
covariates[5, colstoget] <- 9
colstoget <- which(covariates[5,] == "stage iv")
covariates[5, colstoget] <- 10
colstoget <- which(covariates[5,] == "stage x")
covariates[5, colstoget] <- 11
colstoget <- which(is.na(covariates[5,]))
covariates[5, colstoget] <- 12

####################################################
# Subset PCs to be patients with European ancestry #
####################################################
# Read in PCs
allPCs <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/plink.eigenvec")

# Subset to intersection
## format PCs TCGA barcodes
allPCs$V2 <- substr(allPCs$V2, 1, 12)
intersection <- RVenn::overlap(RVenn::Venn(list(allPCs$V2, EDpatients$intersection)))
PCs <- allPCs[which(allPCs$V2 %in%  intersection),]

#######################
# Make PCs covariates #
#######################
PCs <- PCs[,-1]
PCs <- t(PCs)

colnames(PCs) <- PCs[1,]
PCs <- PCs[-1,]

rownames(PCs) <- NULL
PCs <- cbind("id" = make.names(names = rep("PC", dim(PCs)[1]+1), unique = TRUE)[-1], PCs)

common_column_names <- intersect(names(covariates), colnames(PCs))
covariates <- merge(covariates, PCs, by = common_column_names, all = TRUE)

# Write covariates file
data.table::fwrite(covariates, "formattedCovariates.csv")

cat("Done formatting covariates\n")
# Clear space
rm(list=ls())
