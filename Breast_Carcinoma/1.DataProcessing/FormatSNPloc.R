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
