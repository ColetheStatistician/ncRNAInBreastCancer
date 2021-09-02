# Get intersection of all European patients from clinical covariates, geneotype, and expression data
cat("Getting patient intersection\n")

# libraries
library(data.table)
library(magrittr)
library(bigsnpr)
library(RVenn)

# Read in data
clinical <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/Clinical/Breast_Carcinoma_Covariates.txt")
genotype <- bigsnpr::snp_attach(bigsnpr::snp_readBed2("/u/scratch/n/nolanc/Breast_Carcinoma/genotypes/22.Breast_Carcinoma.bed"))
data <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/mRNA/tumor_Breast_Carcinoma_Expression.txt")
miRNA <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/miRNA/Breast_Carcinoma_Expression.txt")
#n_expression <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/mRNA/normal_Breast_Carcinoma_Expression.txt")
#t_expression <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/mRNA/tumor_Breast_Carcinoma_Expression.txt")

# Subset clinical covariates to only contain European Ancestry individuals
# Remove non-European ancestry individuals
colstoget <- which(clinical[1358,] == "white")
clinical <- clinical[,..colstoget]

# Get intersect of miRNA and RNA patients before getting intersection with genotype and clinical patients
## clean
### first row has no necessary information
data <- data[-1,]

## Change column names to be easier to match samples between datasets
colnames(data)[1] <- "gene_id"
colnames(miRNA)[1] <- "gene_id"

data_part <- colnames(data)
mi_part <- colnames(miRNA)

data_samples <- c(data_part[1], substr(data_part[-1], 1, 16))
mi_samples <- c(mi_part[1], substr(mi_part[-1], 1, 16))

colnames(data) <- data_samples
colnames(miRNA) <- mi_samples

## subset columns to patients that are in both miRNA and data
expression <- data %>%
        dplyr::select(gene_id, intersect(data_samples, mi_samples))
mi_intersection_col <- miRNA %>%
        dplyr::select(gene_id, intersect(data_samples, mi_samples))

# Get patient ids
## Make all upper case to prevent any issues with matching
p_clinical <- toupper(clinical[20,])
p_genotype <- toupper(genotype$fam$sample.ID)
p_expression <- toupper(colnames(expression))

# Cut to first 3 TCGA identification sections
p_clinical_sub <- substr(p_clinical, 1, 12)
p_genotype_sub <- substr(p_genotype, 1, 12)
p_expression_sub <- substr(p_expression, 1, 12)

# Get intersection of all patients
intersection <- RVenn::overlap(RVenn::Venn(list(p_clinical_sub, p_genotype_sub, p_expression_sub)))

# Output
data.table::fwrite(data.frame(intersection), "uniquePatients.csv")

# Delete .rds/.bk files
file.remove(list.files()[stringr::str_detect(list.files(), ".rds")])
file.remove(list.files()[stringr::str_detect(list.files(), ".bk")])

cat("Done getting patient intersection\n")
# Clear space
rm(list=ls())
