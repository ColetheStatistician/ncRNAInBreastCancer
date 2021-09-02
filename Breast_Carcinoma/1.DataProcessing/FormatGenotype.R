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
