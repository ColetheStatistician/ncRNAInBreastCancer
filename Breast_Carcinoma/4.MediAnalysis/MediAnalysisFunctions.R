# load libraries
library(boot)
library(magicfor)
library(tictoc)
library(magrittr)
library(dplyr)

# Test Data
# pc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_Covariates.csv")[,-1]
# nc_covariates <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_Covariates.csv")[,-1]
# pc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv"))
# nc_exp <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv"))
# genotypes <- data.frame(data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/tot.Breast_Carcinoma_genotype.csv")) # CHANGE ME
# SNPS to get from CHR22
# optPC <- readRDS("pc_final_trans_results.rds")
# optNC <- readRDS("nc_final_results.rds")
## Get list of unique SNPs in intersection between the pc cis-eQTLs and nc trans-eQTLs
# uniq_snps <- intersect(optPC$eQTLs$trans$eqtls$snps, optNC$eQTLs$cis$eqtls$snps)
## Table of top p value trans qtls with overlapping cis
# trans_w_overlapping_cis <- optPC$eQTLs$trans$eqtls[which(optPC$eQTLs$trans$eqtls$snps %in% uniq_snps),]
# CHR22snps <- trans_w_overlapping_cis[which(trans_w_overlapping_cis$snps %in% genotypes$id),]
# Read in triplets
# triplets <- read.csv("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/triplets.csv")
# CHR22triplets <- triplets[which(triplets$snps %in% CHR22snps$snps),]


# DEBUGGING/DEV 
## SNP pc-gene pairs have multiple nc-RNAs
#multiple <- triplets %>%
#  group_by(snps, pc_trans_genes) %>%
#  count() %>%
#  arrange(desc(n)) %>%
#  filter(n>1)
# Univariate mediation analysis situation
#singles <- triplets %>%
#  group_by(snps, pc_trans_genes) %>%
#  count() %>%
#  arrange(desc(n)) %>%
#  filter(n == 1)

# use singles as a simple case to set up analysis
#chr10119033471GA <- triplets %>% filter(snps == 'chr10:119033471:G:A', pc_trans_genes == 'ENSG00000070785')
# multiple corresponding NC example
#chr1955083602GA <- triplets %>% filter(snps == 'chr19:55083602:G:A', pc_trans_genes == 'ENSG00000168062')
#########################################################################################################################################################
# Function to extract SNP vector corresponding to genotype, regress continuous covariates our of expression matrix, and then permute mediation analysis #
#########################################################################################################################################################
# Argument values
#desired_SNP = 'chr19:55113783:T:C'; pc_desired_ENSEMBL ='ENSG00000007392'; triplets_df = triplets;
#genotypes = genotypes; expressionMatrix = pc_exp; mediators_matrix = nc_exp; pc_covariates = pc_covariates; nc_covariates = nc_covariates; 

Mediation <- function(desired_SNP, pc_desired_ENSEMBL, triplets_df = triplets,
                      genotypes, expressionMatrix = pc_exp, mediators_matrix = nc_exp,
                      pc_covariates, nc_covariates, scenario = 0, total_scenarios = 0) {
  tictoc::tic()
  
  # 1. Get SNP s genotypes
  snp_vec <- as.numeric(t(genotypes[which(genotypes$id == desired_SNP),])[-1])
  
  # 2. Get pc-gene G expression values (regress out covariates)
  ## Get expression values corresponding to desired_SNP s and pc-gene G (should only be one)
  matchingPCgenes <- which(expressionMatrix$id == pc_desired_ENSEMBL)
  if (length(matchingPCgenes) > 1) {
    warning("There are multiple matching annotations for the protein coding genes")
  } else {
    exp_vec <- as.numeric(t(expressionMatrix[which(expressionMatrix$id == pc_desired_ENSEMBL),])[-1])
  }
  ## Regress cont. covariates out of expression vector
  expdf <- data.frame(t(pc_covariates), "expression" = exp_vec)
  clean_exp_vec <- lm(expression ~ ., data=expdf)$residuals
  
  # 3. Get ALL nc-RNA expression values that are in a triplet with these other two (regress out covariates)
  ## Get (potentially) multiple expression values corresponding to desired_SNP and nc RNAs
  nc_desired_ENSEMBL_index <- which((triplets_df$snps == desired_SNP) & (triplets_df$pc_trans_genes == pc_desired_ENSEMBL))
  ## TODO: Is this the best solution?
  ## If num_mediators > 10, subset to the ten most significant
  if (length(nc_desired_ENSEMBL_index) > 10) {
    nc_desired_ENSEMBL_s <- triplets_df[nc_desired_ENSEMBL_index,] %>% dplyr::arrange(nc_cis_pvalue) %>% dplyr::slice(1:10) %>% dplyr::select(nc_cis_genes) %>% unlist()
  } else {
    nc_desired_ENSEMBL_s <- triplets_df$nc_cis_genes[nc_desired_ENSEMBL_index]
  }
  
  matchingNCgenesindex <- which(mediators_matrix$id %in% nc_desired_ENSEMBL_s)
  if (length(matchingNCgenesindex) > 1) {
    warning(paste0("Multiple nc-RNA for SNP ", desired_SNP, " and protein coding gene ", pc_desired_ENSEMBL))
    ## Get mediators expression matrix
    med_exp <- mediators_matrix[matchingNCgenesindex,-1]
    ## Get covariate names for linear model
    geneNames <- mediators_matrix[matchingNCgenesindex,1]
    rownames(med_exp) <- geneNames
    rownames(nc_covariates) <- make.names(rep("cov", dim(nc_covariates)[1]+1), unique = TRUE)[-1]
    covNames <- rownames(nc_covariates)
    ## Regress cont. covariates out of mediator vector
    meddf <- data.frame(t(med_exp), t(nc_covariates))
    storage <- list()
    for (j in 1:length(geneNames)) {
      ### Fit model and get residuals
      storage[[j]] <- lm(
        as.formula(paste(colnames(meddf)[j], "~",
          paste(colnames(meddf)[(dim(meddf)[2] - length(covNames)+1):dim(meddf)[2]], collapse = "+"),
          sep = ""
      )),
         data = meddf
       )$residuals
    }
    ## Get regressed out expression values
    clean_med <- Reduce(cbind, storage) 
    ## TODO: remove collinear components
    clean_med <- dropCollinearCov(mediator_matrix = clean_med)
  } else {
    med_vec <- as.numeric(t(mediators_matrix[matchingNCgenesindex,])[-1]) ## TODO: Change for multiple nc-RNAs
    ## Regress cont. covariates out of mediator vector
    meddf <- data.frame(t(nc_covariates), "mediator_exp" = med_vec)
    clean_med <- lm(mediator_exp ~ ., data=meddf)$residuals
  }
  
  # 4. Mediation analysis
  # Get mediation effect
  # snp = unlist(snp_vec); expression = unlist(clean_exp_vec); mediators = clean_med; covs = NULL; nperms = 1000; parallel = 'no'; nc = 1
  
  mediation_effect <- permuteTME(snp = unlist(snp_vec),
                                 expression = unlist(clean_exp_vec),
                                 mediators = clean_med,
                                 covs = NULL,
                                 nperms = 1000,
                                 parallel = 'no',
                                 nc = 1)
  
  cat("\n")
  cat(paste0("Scenario:", scenario, " / ", total_scenarios))
  cat("\n")
  tictoc::toc()
  cat("\n")
  
  print(list("test.stat" = mediation_effect$test.stat,
              "p.value" = mediation_effect$p.value,
              "ncRNAs_in_test" = nc_desired_ENSEMBL_s))
  
  return(list("test.stat" = mediation_effect$test.stat,
              "p.value" = mediation_effect$p.value,
              "ncRNAs_in_test" = nc_desired_ENSEMBL_s))
}

# Function to recursively drop collinear columns in matrix of numeric covariates
dropCollinearCov <- function(mediator_matrix) {
  
  if (is.null(dim(mediator_matrix)[2])) { # If is reduced to a 1-dim vector
    return(mediator_matrix)
  }
  
  if (sum(colSums(cor(mediator_matrix)==1) > 1) > 1) {
    collinear_index <- which(colSums(cor(mediator_matrix) == 1) > 1)
    # recursively drop collinear column
    return(dropCollinearCov(mediator_matrix[,-collinear_index[length(collinear_index)]]))
  } else {
    return(mediator_matrix)
  }
}

#' Compute total mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total mediation effect with boots structure
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#' @param indices blank, index for boot package
#'
#' @return estimate of TME
#'
#' @export
# Debug:
#snp = unlist(snp_vec); expression = (clean_exp_vec); mediators = (clean_med_vec); covs = NULL;

computeTME <- function(snp,
                       expression,
                       mediators,
                       covs,
                       indices){
  snp = c(snp)
  snp = snp[indices]

  if (ncol(data.frame(mediators)) == 0){return(0)}

  if (!is.null(covs)){
    direct = lm(expression ~ snp + mediators + covs)
    indirect = lm(mediators ~ snp + covs)
    }
  if (is.null(covs)){
    direct = lm(expression ~ snp + mediators)
    indirect = lm(mediators ~ snp)
  }
  b = coef(direct)[grepl('mediators',names(coef(direct)))]
  if (ncol(data.frame(mediators)) > 1){
    a = coef(indirect)[2,]
  } else {
    a = coef(indirect)[2]
  }

  TME = a %*% b

  return(as.numeric(TME))
}

#' Perform permutation test for the total mediation effect
#'
#' The function takes in a SNP-mediator set-gene triplet
#' and computes the total absolute mediation effect and permutation
#' test P-value
#'
#' @param snp vector, SNP of interest
#' @param expression vector, gene expression of interest
#' @param mediators data frame, mediators of interest
#' @param covs data frame, covariates
#' @param nperms integer, number of permutations for the null distribution
#' @param parallel character, boot parallel input
#' @param nc integer, number of cores
#'
#' @return estimate of TME and the permutation P-value
#'
#' @importFrom boot boot
#'
#' @export
permuteTME = function(snp,
                      expression,
                      mediators,
                      covs,
                      nperms = 1000,
                      parallel = 'no',
                      nc){

  a = boot::boot(data = snp,
                 statistic = computeTME,
                 R = nperms,
                 sim = 'permutation',
                 expression = expression,
                 mediators = mediators,
                 covs = covs,
                 parallel = parallel,
                 ncpus = nc)
  p = (sum(abs(a$t) >= abs(a$t0))+1)/(nperms + 1)
  return(list(test.stat = a$t0,
              p.value = p))

}
