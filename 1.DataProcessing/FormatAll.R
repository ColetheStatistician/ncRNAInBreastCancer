# Format the data in one script in the following sequence:
# 1. Get intersection of European Ancestry patients in all data sets
# 2. Format covariate file that matches this format: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt
# 3. Format genotypes
# 4. Formatting Expression and Gene Loc
# 5. Format SNP locations

#########################################################################################################
# 1. Get intersection of all European patients from clinical covariates, geneotype, and expression data #
#########################################################################################################
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

###############################################################################################################################################
# 2. Format covariate file that matches this format: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt #
###############################################################################################################################################
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

#######################
# 3. Format genotypes #
#######################
cat("Formatting genotypes\n")
# load libraries
require(bigsnpr)
require(data.table)
require(stringr)
require(tictoc)

#bedfile <- "22.Breast_Carcinoma.bed"
#snp <- bigsnpr::snp_attach(bigsnpr::snp_readBed2(bedfile))
# Get European Descent (ED) participants
EDpatients <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/uniquePatients.csv")

# read in data
bedfile <- c("/u/scratch/n/nolanc/Breast_Carcinoma/genotypes/22.Breast_Carcinoma.bed",
             "/u/scratch/n/nolanc/Breast_Carcinoma/genotypes/tot.Breast_Carcinoma.bed") # list.files()[stringr::str_detect(list.files(), ".bed")]

for (x in 1:length(bedfile)) {
  file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".bk")) # Deleting files here
  file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".rds"))
  tictoc::tic()
  print(paste0("Processing ", bedfile[x]))

  # Read in file
  snp <- bigsnpr::snp_attach(bigsnpr::snp_readBed2(bedfile[x]))

  # Subset to patients in all datasets
  ## Get genotype patients
  gPatients <- snp$fam$sample.ID
  ## Subset whole gPatients ID to match EDpatient IDs
  gPatients_sub <- substr(gPatients, 1, 12)
  ## subset our bigsnpr object to patients in all datasets
  backupfile <- bigsnpr::snp_subset(x = snp, ind.row = which(gPatients_sub %in% as.matrix(EDpatients)), backingfile=NULL)
  ## Read in subset file
  subsnp <- readRDS(backupfile)

  # Format
  ## Take transpose of snp$genotypes
  subsnp$genotypes_transpose <- big_transpose(big_copy(subsnp$genotypes))

  ## Convert snp$genotypes_tranpose from an S4 object to a matrix
  subsnp$genotypes_transpose <- subsnp$genotypes_transpose[]

  ## Add column names (sample ID) to the matrix snp$genotypes_tranpose
  colnames(subsnp$genotypes_transpose) <- subsnp$fam$sample.ID

  ## Convert to data frame and add first column with SNP ids
  snp_data <- data.frame("id" = snp$map$marker.ID, subsnp$genotypes_transpose)
  
  # Identify non-unique genes and make annotations unique
  repeatedGenes <- snp_data$id[duplicated(snp_data$id)]
  toChange <- which(snp_data$id %in% repeatedGenes)
  snp_data$id[toChange] <- make.names(snp_data$id[toChange], unique = TRUE)

  # Write data
  data.table::fwrite(snp_data, paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_genotype.csv"))
  # erase extra .bk and .rds files to conserve memory
  # snp .bk, .rds files
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".bk"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".bk"))
  }
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".rds"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".rds"))
  }
  # subsnp .rds, .bk files
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.bk"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.bk"))
  }
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.rds"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.rds"))
  }

  tictoc::toc()

  # prevent memory leak/overload
  rm(list = c("snp", "subsnp", "snp_data",
              "backupfile", "gPatients", "gPatients_sub"))
}

cat("Done formatting genotypes\n")
# Clear space
rm(list=ls())

#########################################
# 4. Formatting Expression and Gene Loc #
#########################################
cat("Formatting Expression and Gene Loc\n")
# Load libraries
library(data.table)
library(stringr)
library(biomaRt)

# Expression files to format
file <- "/u/scratch/n/nolanc/Breast_Carcinoma/mRNA/tumor_Breast_Carcinoma_Expression.txt"
rawfile <- "/u/scratch/n/nolanc/Breast_Carcinoma/mRNA/BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt"
miRNAfile <-"/u/scratch/n/nolanc/Breast_Carcinoma/miRNA//Breast_Carcinoma_Expression.txt"
# Load data
txtfile <- data.table::fread(file)
raw <- data.table::fread(rawfile)
miRNA <- data.table::fread(miRNAfile)

## first row has no necessary information
raw <- raw[-1,]

# Get European Descent (ED) patients that are in all other datasets
## Get ED patients
EDpatients <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/uniquePatients.csv")
## get expression patient IDs
p_exp <- stringr::str_replace_all(colnames(txtfile), "\\.", "-")
## there are 7 individuals with multiple tumor samples
## Exclude metastatic tumors ("06A")
p_exp <- p_exp[!stringr::str_detect(p_exp, "06A")]
## Cut to first 3 TCGA identification sections
p_exp_sub <- substr(p_exp, 1, 12)
## subset to patients in all datasets
t_patients <- c("id", p_exp[which(p_exp_sub %in% as.matrix(EDpatients))])

# Combine miRNA expression data with RNA expressiond data #
## Exclude rows comprised entriely of 0s
exclude <- -which(rowSums(miRNA[,-1]) < 1)
miRNA <- miRNA[exclude,]
## Normalize expression miRNA expression counts
### log2(expr + 1) transform then standardize
geneNormalize <- function(vector) {
 transformed <- log2(vector + 1)
 (transformed - mean(transformed)) / sd(transformed)
}
normalized_miRNA <- data.table(cbind(miRNA[,1], t(apply(miRNA[,-1], MARGIN=1, FUN=geneNormalize))))

# Change column names to be easier to match samples between datasets
colnames(raw)[1] <- "id"
colnames(normalized_miRNA)[1] <- "id"

raw_patients <- colnames(raw)
mi_patients <- colnames(normalized_miRNA)

raw_patients_sub <- c(raw_patients[1], substr(raw_patients[-1], 1, 16))
mi_patients_sub <- c(mi_patients[1], substr(mi_patients[-1], 1, 16))

# subset columns to patients that are in both miRNA and data
raw_colsToGet <- which(raw_patients_sub %in% intersect(raw_patients_sub, mi_patients_sub))
raw_sub <- raw[,..raw_colsToGet]
mi_colsToGet <- which(mi_patients_sub %in% intersect(raw_patients_sub, mi_patients_sub))
mi_sub <- normalized_miRNA[,..mi_colsToGet]
# raw_sub and mi_sub columns are in same order with regards to patients
colnames(mi_sub) <- colnames(raw_sub)

# Get raw, normalized expression using ED tumor patient ids
temp <- raw_sub[,1]
raw_colsToGet <- which(colnames(raw_sub) %in% t_patients)[-1] # Exclude first column
raw_sub <- cbind(temp, data.matrix(raw_sub[,..raw_colsToGet]))

mi_colsToGet <- which(colnames(mi_sub) %in% t_patients)
mi_sub <- mi_sub[,..mi_colsToGet]

## Remove genes with non unique counts
notunique <- -which(apply(raw_sub, 1, function(x) length(unique(x))) < 3)
raw_sub <- raw_sub[notunique,]

## Normalize and standardize non miRNA expression data
raw_sub <- data.table(cbind(raw_sub[,1], t(apply(raw_sub[,-1], MARGIN=1, FUN=geneNormalize))))

# Convert from entrez annotation to ENSEMBL annotations
## miRNA first
## Get hg19 genome (solution from https://www.biostars.org/p/136775/)
grch37 <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

## convert from hsa-mi format to ENSEMBL
mi_genes <- biomaRt::getBM(filters = "mirbase_id",
                           attributes = c("mirbase_id", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype"),
                           values = mi_sub$id,
                           mart = grch37)

## subset to gene_biotype Arjun recommended
mi_genes <- mi_genes[which(mi_genes$gene_biotype %in% c("lincRNA", "miRNA","snoRNA")),]
head(mi_genes)

## Map observed mirbase annotations to enesembl annotations
mi_sub <- mi_sub[match(mi_genes$mirbase_id, mi_sub$id),]
## Convert ensembl annotations
mi_sub$id <- mi_genes$ensembl_gene_id

# Now non miRNA
## Fix gene id column to be readable by biomaRt by getting the Entrez Gene ID (numbers after |)
raw_sub$id <- sub(".*\\|", "", raw_sub$id)

## Convert from entrez to ensembl annotations
raw_genes <- biomaRt::getBM(filters = "entrezgene_id",
                            attributes = c("entrezgene_id", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype"),
                            values = raw_sub$id,
                            mart = grch37)
head(raw_genes)

## Identify non-coding vs coding genes
### Gene_types in annotations
ncRNAs <- c("tRNA", "rRNA", "scRNA", "snRNA", "snoRNA", "miRNA", "misc_RNA", "lincRNA", # http://useast.ensembl.org/info/genome/genebuild/ncrna.html: tRNA, Mt-tRNA,rRNA, scRNA, snRNA, snoRNA, miRNA, misc_RNA, lincRNA
            "sRNA", "scaRNA") # https://www.gencodegenes.org/pages/biotypes.html:
pcRNAs <- c("protein_coding")
### Partition annotations to nc and pc genes
ncGenes <- raw_genes[which(raw_genes$gene_biotype %in% ncRNAs),]
pcGenes  <- raw_genes[which(raw_genes$gene_biotype %in% pcRNAs),]

## Map observed entrez annotations to enesembl annotations
nc_raw_sub <- raw_sub[match(ncGenes$entrezgene_id, raw_sub$id),]
pc_raw_sub <- raw_sub[match(pcGenes$entrezgene_id, raw_sub$id),]
## Convert ensembl annotations
nc_raw_sub$id <- ncGenes$ensembl_gene_id
pc_raw_sub$id <- pcGenes$ensembl_gene_id
## Combine miRNA with nc genes
nc_raw_sub <- rbind(nc_raw_sub, mi_sub)

# Make gene IDs unique
# Identify non-unique genes and make annotations unique
nc_repeatedGenes <- nc_raw_sub$id[duplicated(nc_raw_sub$id)]
nc_toChange <- which(nc_raw_sub$id %in% nc_repeatedGenes)
nc_raw_sub$id[nc_toChange] <- make.names(nc_raw_sub$id[nc_toChange], unique = TRUE)

pc_repeatedGenes <- pc_raw_sub$id[duplicated(pc_raw_sub$id)]
pc_toChange <- which(pc_raw_sub$id %in% pc_repeatedGenes)
pc_raw_sub$id[pc_toChange] <- make.names(pc_raw_sub$id[pc_toChange], unique = TRUE)

# Output formatted files
## expression files
data.table::fwrite(nc_raw_sub, "formatted_nc_expression.csv")
data.table::fwrite(pc_raw_sub, "formatted_pc_expression.csv")

## gene locations
gene_locations <- data.frame("id" = c(ncGenes$ensembl_gene_id, pcGenes$ensembl_gene_id, mi_genes$ensembl_gene_id),
                             "chr" = c(ncGenes$chromosome_name, pcGenes$chromosome_name, mi_genes$chromosome_name),
                             "s1" = c(ncGenes$start_position, pcGenes$start_position, mi_genes$start_position),
                             "s2" = c(ncGenes$end_position, pcGenes$end_position, mi_genes$end_position))
data.table::fwrite(gene_locations, "formatted_gene_loc.csv")

cat("Done formatting Expression and Gene Loc\n")

# Clear space
rm(list=ls())

###########################
# 5. Format SNP locations #
###########################                      
cat("Formatting SNP locations\n")
# load libraries
require(bigsnpr)
require(data.table)
require(stringr)
require(tictoc)

# Get European Descent (ED) participants
EDpatients <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/uniquePatients.csv")

# read in data
bedfile <- c("/u/scratch/n/nolanc/Breast_Carcinoma/genotypes/22.Breast_Carcinoma.bed",
             "/u/scratch/n/nolanc/Breast_Carcinoma/genotypes/tot.Breast_Carcinoma.bed") # list.files()[stringr::str_detect(list.files(), ".bed")]

for (x in 1:length(bedfile)) {
  tictoc::tic()
  print(paste0("Processing ", bedfile[x]))
  # Read in file
  snp <- bigsnpr::snp_attach(bigsnpr::snp_readBed2(bedfile[x]))

  # Subset to patients in all datasets
  ## Get genotype patients
  gPatients <- snp$fam$sample.ID
  ## Subset whole gPatients ID to match EDpatient IDs
  gPatients_sub <- substr(gPatients, 1, 12)
  ## subset our bigsnpr object to patients in all datasets
  backupfile <- bigsnpr::snp_subset(x = snp, ind.row = which(gPatients_sub %in% as.matrix(EDpatients)), backingfile=NULL)
  ## Read in subset file
  subsnp <- readRDS(backupfile)

  ## Convert to data.table with appropriate columns
  snp_loc <- data.table("snp" = subsnp$map$marker.ID, "chr" = subsnp$map$chromosome, "pos" = subsnp$map$physical.pos)
  
  # Identify non-unique genes and make annotations unique
  repeatedGenes <- snp_loc$snp[duplicated(snp_loc$snp)]
  toChange <- which(snp_loc$snp %in% repeatedGenes)
  snp_loc$snp[toChange] <- make.names(snp_loc$snp[toChange], unique = TRUE)

  # Write data
  data.table::fwrite(snp_loc, paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_snploc.csv"))
  # erase extra .bk and .rds files to conserve memory
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".bk"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".bk"))
  }
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".rds"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), ".rds"))
  }
  # subsnp .rds, .bk files
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.bk"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.bk"))
  }
  if(file.exists(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.rds"))) {
    file.remove(paste0(substr(bedfile[x], 1, nchar(bedfile[x])-4), "_sub1.rds"))
  }
  tictoc::toc()
  # prevent memory overload
  rm(list = c("snp", "subsnp", "snp_loc",
              "backupfile", "gPatients", "gPatients_sub"))
}
# Clear space
rm(list=ls())
cat("Done formatting SNP locations\n")
