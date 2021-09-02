# Get intersection of all patients with race covariates and expression data
# libraries
library(data.table)
library(ggplot2)

# Read in data
clinical <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/Clinical/Breast_Carcinoma_Covariates.txt") 
allPCs <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/plink.eigenvec")

## Remove NA ethinicity patients
colsToGet <- which(!is.na(clinical[1358,]))
clinical <- clinical[,..colsToGet]

# Get intersection of PCA patients and covariate patients
## Make patient ID the column name
p_clinical <- toupper(clinical[20,])
p_pcs <- as.matrix(allPCs[,2])

## Cut to first 3 TCGA identification sections
p_clinical_sub <- substr(p_clinical, 1, 12)
p_pcs_sub <- substr(p_pcs, 1, 12)

## Get intersection
intersection <- intersect(p_clinical_sub, p_pcs_sub)

## Subset PCs to patients in intersection
PCs <- allPCs[which(p_pcs_sub %in% intersection),]
ancestries <- t(data.frame(clinical)[1358, which(p_clinical_sub %in% intersection)])

## Make dataframe of PC1, PC2, ancestry
df <- data.frame(
    PCs[,3],
    PCs[,4],
    ancestries
)
colnames(df) <- c("PC1", "PC2", "Ancestry")
df[,1] <- as.numeric(df[,1])
df[,2] <- as.numeric(df[,2])
df[,3] <- as.factor(df[,3])

# Plot PC1 vs PC2 colored by ancestry
pdf("PC1vsPC2.pdf")
ggplot(data=df, aes(x=PC1, y=PC2, color=Ancestry)) + geom_point()
dev.off()
