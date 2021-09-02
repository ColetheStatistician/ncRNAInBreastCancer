# Script to generate 2 .csv files seperated by coding and noncoding RNAs (identified via http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz)
# Goal formatting from http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/geneloc.txt

# Load libraries
library(data.table)
library(dplyr)
library(biomaRt)
library(magrittr)

# Read in data
data <- data.table::fread("BRCA.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt")
annotations <- data.table::fread("gene_annotation.gencode.hg19.tsv")
miRNA <- data.table::fread("../miRNA/Breast_Carcinoma_Expression.txt")

# clean
## snip off transcript number at end (e.g. ENSEMBL######.transcript#)
annotations$ensembl_gene_id <- gsub("\\..*","", annotations$gene_id)
## first row has no necessary information
data <- data[-1,]

#######################
# miRNA incorporation #
#######################
# Normalize expression miRNA expression counts
## log2(expr + 1) transform then standardize
geneNormalize <- function(vector) {
 transformed <- log2(vector + 1)
 (transformed - mean(transformed)) / sd(transformed)
}
normalized_miRNA <- data.table(cbind(miRNA[,1], t(apply(miRNA[,-1], MARGIN=1, FUN=geneNormalize))))

# Change column names to be easier to match samples between datasets
colnames(data)[1] <- "gene_id"
colnames(normalized_miRNA)[1] <- "gene_id"

data_part <- colnames(data)
mi_part <- colnames(normalized_miRNA)

data_samples <- c(data_part[1], substr(data_part[-1], 1, 16))
mi_samples <- c(mi_part[1], substr(mi_part[-1], 1, 16))

colnames(data) <- data_samples
colnames(normalized_miRNA) <- mi_samples

# subset columns to patients that are in both miRNA and data
data_intersection_col <- data %>%
	dplyr::select(gene_id, intersect(data_samples, mi_samples))
mi_intersection_col <- normalized_miRNA %>%
	dplyr::select(gene_id, intersect(data_samples, mi_samples))

# make sure data_intersection_col has the same column order as mi_intersection_col
colnames(data_intersection_col) <- colnames(data_intersection_col)[match(colnames(mi_intersection_col), colnames(data_intersection_col))]

#######################################
# Subset to European Descent patients #
#######################################
# Get ED participants
EDpatients <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/genotypes/uniquePatients.csv")

# get patients
d_patients <- colnames(data_intersection_col)

# Subset whole patient ID to match EDpatient IDs
d_patients_sub <- substr(d_patients, 1, 12)

## subset to patients in all datasets
data_intersection_col <- data.frame("gene_id" = data_intersection_col[,1], data.frame(data_intersection_col)[,which(d_patients_sub %in% as.matrix(EDpatients))])
mi_intersection_col <- data.frame("gene_id" = mi_intersection_col[,1], data.frame(mi_intersection_col)[,which(d_patients_sub %in% as.matrix(EDpatients))])

##################################
# Convert to ENSEMBL annotations #
##################################
# Get hg19 genome (solution from https://www.biostars.org/p/136775/)
grch37 <- biomaRt::useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

# convert from hsa-mi format to ENSEMBL
mi_genes <- biomaRt::getBM(filters = "mirbase_id",
 			   attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype"),
                           values = mi_intersection_col$gene_id,
			   mart = grch37)
## subset to gene_biotype Arjun recommended
mi_genes <- mi_genes[which(mi_genes$gene_biotype %in% c("lincRNA", "miRNA","snoRNA")),]
head(mi_genes)

# Convert from Entrez gene annotations to ENSEMBL with hg19/gch37 genome
## Fix gene id column to be readable by biomaRt by getting the Entrez Gene ID (numbers after |)
data_intersection_col$entrezgene_id <- sub(".*\\|", "", data_intersection_col$gene_id)

## Convert 
genes <- biomaRt::getBM(filters = "entrezgene_id", 
                        attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype"), 
                        values = data_intersection_col$entrezgene_id, 
                        mart = grch37)    
head(genes)

## Identify non-coding vs coding genes
### Gene_types in annotations
ncRNAs <- c("tRNA", "Mt_tRNA", "rRNA", "scRNA", "snRNA", "snoRNA", "miRNA", "misc_RNA", "lincRNA", # http://useast.ensembl.org/info/genome/genebuild/ncrna.html: tRNA, Mt-tRNA,rRNA, scRNA, snRNA, snoRNA, miRNA, misc_RNA, lincRNA
            "Mt_rRNA", "ribozyme", "sRNA", "scaRNA") # https://www.gencodegenes.org/pages/biotypes.html: 
pcRNAs <- c("protein_coding")
### Partition annotations to nc and pc genes
ncGenes <- genes[which(genes$gene_biotype %in% ncRNAs),]
pcGenes  <- genes[which(genes$gene_biotype %in% pcRNAs),]

# Combine miRNA and ncGenes info
ncGenes <- rbind(ncGenes, mi_genes)

# Prep for writing
## Subset to necessary columns
ncGenes <- ncGenes[, -5]
colnames(ncGenes) <- c("geneid", "chr", "s1", "s2")
pcGenes <- pcGenes[, -5]
colnames(pcGenes) <- c("geneid", "chr", "s1", "s2")
## Order columns in descending order
ncGenes <- ncGenes %>% arrange(chr, s1)
pcGenes <- pcGenes %>% arrange(chr, s1)

head(ncGenes)
dim(ncGenes)
head(pcGenes)
dim(pcGenes)

# Write formatted gene location files
data.table::fwrite(ncGenes, "formatted_ncGenes_loc.csv")
data.table::fwrite(pcGenes, "formatted_pcGenes_loc.csv")

