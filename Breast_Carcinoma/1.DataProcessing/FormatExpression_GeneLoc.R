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
