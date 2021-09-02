# Code to parse through and get information from annotations
# data obtained from: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf.gz
require(data.table)
require(stringr)
require(readr)
# Read in data
genes <- data.table::fread("gencode.v19.chr_patch_hapl_scaff.annotation.gtf")
# Correct column names
setnames(genes, names(genes), c("chr","source","type","start","end","score","strand","phase","attributes"))
# Parse attributes
# Solution taken from: https://www.biostars.org/p/272889/
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}

genes$gene_id <- extract_attributes(gtf_attributes=genes$attributes, att_of_interest="gene_id")
genes$transcript_id <- extract_attributes(gtf_attributes=genes$attributes, att_of_interest="transcript_id")
genes$gene_type <- extract_attributes(gtf_attributes=genes$attributes, att_of_interest="gene_type")
genes$gene_name <- extract_attributes(gtf_attributes=genes$attributes, att_of_interest="gene_name")

# Drop other attribute information
genes <- genes[,-"attributes"]

# Output zipped file
fwrite(parsed_df,'gencode.v19.tsv.gz')
