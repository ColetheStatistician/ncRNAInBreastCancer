# Load libraries
library(data.table)
library(Rhcpp)
library(MatrixEQTL)
library(future.apply)

####################################
# Function to format data for hcps #
####################################
formatHCPs <- function(covariates, expression_matrix) {
  # ?hcp
  # Make Z Matrix: nxd where n=number of subjects and d=the number of known covariates
  Z <- data.frame(t(covariates))
  colnames(Z) <- Z[1,]
  Z <- Z[-1,]
  ## Standardize age
  Z$patient.age_at_initial_pathologic_diagnosis <- as.numeric(Z$patient.age_at_initial_pathologic_diagnosis)
  Z$patient.age_at_initial_pathologic_diagnosis <- (Z$patient.age_at_initial_pathologic_diagnosis - mean(Z$patient.age_at_initial_pathologic_diagnosis)) / sd(Z$patient.age_at_initial_pathologic_diagnosis)
  
  # Make Y Matrix: nxg where n=number of subjects and g are genes
  ## Make unique rownames 
  rows <- make.names(as.matrix(expression_matrix[,1]), unique = TRUE)
  expression_matrix <- data.frame(expression_matrix[,-1])
  rownames(expression_matrix) <- rows
  
  ## Make nxg
  expression_matrix <- t(expression_matrix)
  
  ## Remove constant genes
  if (length(which(apply(expression_matrix,2,sd) == 0)) > 0) {
    expression_matrix <- expression_matrix[,which(apply(expression_matrix,2,sd) != 0)]
  }
  
  return(list("Z" = Z,
              "Y" = expression_matrix))
}

########################
# Function to get HCPs #
########################

r_hcp <- function(Z, Y, k, lambda1, lambda2, lambda3, iter=100) {
  ## convergence criteria
  tol <- 1e-6
  
  A <- matrix(0, ncol(Z), k)
  W <- matrix(0, nrow(Z), k)
  diag(W) <- 1
  n <- k*ncol(Y)
  ##B <- matrix(runif(n), ncol(Z), ncol(Y))
  B <- matrix((1:n)/n, k, ncol(Y))
  
  n1 <- nrow(Z)
  d1 <- ncol(Z)
  
  n2 <- nrow(Y)
  d2 <- ncol(Y)
  
  if(n1 != n2)
    message('number of rows in F and Y must agree')
  
  if (k < 1 | lambda1 < 0 | lambda2 < 0 | lambda3 < 0 )
    message('lambda1, lambda2, lambda3 must be positive and/or k must be an integer')
  
  ##predefine for slight preformance improvement
  diagB <- diag(k)
  diagW <- diag(k)
  diagA <- diag(nrow(A))
  U1 <- lambda1*solve(lambda1*crossprod(Z) + diagA*lambda3)%*%t(Z)
  
  if(iter > 0) {
    o <- numeric(iter)
    for(ii in 1:iter) {
      ##o[ii] <- norm(Y-W%*%B, type="F")^2 + lambda1*norm(W-Z%*%A, type="F")^2 + lambda2*norm(B, type="F")^2 + lambda3*norm(A, type="F")^2
      o[ii] <- sum((Y-W%*%B)^2) + sum((W-Z%*%A)^2)*lambda1 + sum(B^2)*lambda2 + lambda3*sum(A^2)           
      W <- (tcrossprod(Y, B) + lambda1*Z%*%A) %*% solve(tcrossprod(B) + lambda1*diagB)
      B <- solve(crossprod(W) + lambda2*diagW, crossprod(W,Y))            
      ##A <- solve(lambda1*crossprod(Z) + diagA*lambda3), lambda1*crossprod(Z,W))
      A <- U1%*%W
      
      if(ii > 1)  {
        if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
          break
      }
      
    }
  }
  
  list(W=W, B=B, A=A, o=o, iter=ii)
}

# k=5; lambda1=5; lambda2=5; lambda3=5; Z=Z; Y=Y_nc; first_two_plus_k_pcs=4
getHCPs <- function(first_two_plus_k_pcs, # how many of principal components (after the first 5) to use
                    k, lambda1, lambda2, lambda3, Z, Y # hcp variables to optimize
) {
  # Get PCs
  chosenPCs <- make.names(rep("PC", first_two_plus_k_pcs+3), unique = TRUE)[-1] # the `+3` ensures there is at least 2 PCs as covariates
  colsToGet <- c(colnames(Z)[(length(colnames(Z))-4):length(colnames(Z))], chosenPCs)
  Z <- data.matrix(Z[,colsToGet])
  # Get HCPs
  HCPs <- r_hcp(Z = Z, Y = Y,# Z = scale(Z), Y = scale(Y),
                k = k,
                lambda1 = lambda1, 
                lambda2 = lambda2, 
                lambda3 = lambda3)$W
  
  return(HCPs)
}

#########################################
# Function to format HCPS as covariates #
#########################################
formatHCPstoCov <- function(hcps, covariates) {
  # correct dimensions/colnames
  hcps <- t(hcps)
  colnames(hcps) <- colnames(covariates)[-1]
  
  # Names in first column "id"
  hcp_names <- make.names(rep("HCP", dim(hcps)[1]+1), unique = TRUE)[-1]
  hcps <- cbind("id" = hcp_names, hcps)
  
  # Combine with covariates
  covariates <- rbind(covariates, hcps)
  
  return(covariates)
}

#######################################################
# Function to prep data in SlicedData class for MeQTL #
#######################################################
# Refer to: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/R.html#cis
# COV_file_name = "tempcovariates.csv"; GENOTYPE_file_name = genotype; EXP_file_name = exp_nc; gene_location_file_name = geneloc; snps_location_file_name =snploc
MeQTL<- function(
  COV_file_name = "tempcovariates.csv", 
  GENOTYPE_file_name, 
  EXP_file_name, 
  gene_location_file_name, 
  snps_location_file_name
) {
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  # Output file name
  output_file_name_cis = tempfile();
  output_file_name_tra = tempfile();
  
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = 1e-6;
  pvOutputThreshold_tra = 1e-6;
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6;
  
  ## Load genotype data
  snps = SlicedData$new();
  snps$fileDelimiter = ",";      # the , character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(GENOTYPE_file_name);
  
  ## Load gene expression data
  gene = SlicedData$new();
  gene$fileDelimiter = ",";      # the , character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(EXP_file_name);
  
  ## Load covariates
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = ",";      # the , character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(COV_file_name)>0) {
    cvrt$LoadFile(COV_file_name);
  }
  
  ## Run the analysis
  snpspos <- read.csv(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos <- read.csv(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  unlink(output_file_name_tra);
  unlink(output_file_name_cis);
  
  return(me)
  
}

########################################################
# Function that will incorporate hcps and pcs in MeQTL #
########################################################
# k=5; lambda1=5; lambda2=5; lambda3=5; Z=Z; Y=Y_nc; first_two_plus_k_pcs=4
MeQTLandHcp <- function(first_two_plus_k_pcs, # how many principal components (after first 5) to use
                        k, lambda1, lambda2, lambda3, # hcp variables to optimize
                        Y,
                        covariates,  
                        GENOTYPE_file_name, 
                        EXP_file_name, 
                        gene_location_file_name, 
                        snps_location_file_name,
                        scenario_num = 0,
                        total_scenarios = 0) {
  
  ## Prep data
  formattedHCP <- formatHCPs(covariates = covariates, expression_matrix = Y)
  Z <- formattedHCP$Z
  Y <- formattedHCP$Y
  
  # Get HCPs
  HCPs <- getHCPs(first_two_plus_k_pcs = first_two_plus_k_pcs, k = k, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, Z = Z, Y = Y)
  
  # Select desired PCs
  chosenPCs <- make.names(rep("PC", first_two_plus_k_pcs+3), unique = TRUE)[-1] # the `+3` ensures there is at least 2 PCs as covariates
  variablesToGet <- c(c(covariates[,1])$id[(dim(covariates)[1] - 4):dim(covariates)[1]], chosenPCs)
  covariates <- covariates[which(covariates$id %in% variablesToGet),]
  
  # Format HCPS and make covariates
  covariates <- formatHCPstoCov(hcps = HCPs, covariates = covariates)
  
  # save this new covariates matrix in a temporary, new file
  filename <- paste0(rnorm(1), "tempcovariates.csv") # make file name unique
  data.table::fwrite(covariates, filename)
  
  # Perform Matrix eQTL
  me <- MeQTL(COV_file_name = filename, 
              GENOTYPE_file_name = GENOTYPE_file_name, 
              EXP_file_name = EXP_file_name, 
              gene_location_file_name = gene_location_file_name, 
              snps_location_file_name = snps_location_file_name)
  
  cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
  cat("\n")
  cat(paste0("Finished Scenario:", scenario_num, "/", total_scenarios))
  cat("\n")
  
  # Delete temporary covariates file
  file.remove(filename)
  
  return(me)
}

#####################################################################################################
# Function that will use hcps, pcs, and mapply in order to maximize the number of significant eQTLs #
#####################################################################################################
# COV_file_name = "tempcovariates.csv"; GENOTYPE_file_name = genotype; EXP_file_name = exp_nc; gene_location_file_name = geneloc; snps_location_file_name =snploc
# first_two_plus_k_pcs=3; k=2; lambda1=2; lambda2=2; lambda3=2; Z=Z; Y=Y_nc; covariates=covariates
iterateMeQTL <- function(first_two_plus_k_pcs, # how many of first principal components to use
                         k, lambda1, lambda2, lambda3, # hcp variables to optimize
                         Y,
                         covariates, # matrix of covariates
                         GENOTYPE_file_name, 
                         EXP_file_name, 
                         gene_location_file_name, 
                         snps_location_file_name,
                         toOpt = c("k", "lambda1", "lambda2", "lambda3", "first_two_plus_k_pcs"), # Parameters to optimize. Those not listed will be held constant
                         fixPCs = FALSE
                        ) {
  
  # Get permutation of parameters
  if ("k" %in% toOpt) k <- 1:k
  if ("lambda1" %in% toOpt) lambda1 <- 1:lambda1
  if ("lambda2" %in% toOpt) lambda2 <- 1:lambda2
  if ("lambda3" %in% toOpt) lambda3 <- 1:lambda3
  if ("first_two_plus_k_pcs" %in% toOpt) first_two_plus_k_pcs <- 0:first_two_plus_k_pcs

  combos <- expand.grid(lambda1=lambda1,
                        lambda2=lambda2,
                        lambda3=lambda3,
                        k=k,
                        pcs = first_two_plus_k_pcs)
  
  ## If fixing num PCs, subset to rows w max num PCs
  if (fixPCs) combos <- combos[which(combos$pcs == max(combos$pcs)),]
  
  # iterate MeQTL
  results <-  mapply(FUN="MeQTLandHcp",
                     lambda1=combos$lambda1,
                     lambda2=combos$lambda2,
                     lambda3=combos$lambda3,
                     k=combos$k,
                     first_two_plus_k_pcs=combos$pcs,
                     scenario_num = 1:dim(combos)[1],
                     MoreArgs = list(Y=Y,
                                     covariates = covariates,
                                     GENOTYPE_file_name = GENOTYPE_file_name, 
                                     EXP_file_name = EXP_file_name, 
                                     gene_location_file_name = gene_location_file_name, 
                                     snps_location_file_name = snps_location_file_name,
                                     total_scenarios = dim(combos)[1]),
                     SIMPLIFY = FALSE)
  
  return(results)
}

#####################################
# Function to get optimized results #
#####################################
getOptimized <- function(iterateMeQTL_results,
                         first_two_plus_k_pcs, # how many of first principal components to use
                         k, lambda1, lambda2, lambda3, # hcp variables to optimize
                         toOpt = c("k", "lambda1", "lambda2", "lambda3", "first_two_plus_k_pcs"), # Parameters to optimize. Those not listed will be held constant
                         fixPCs = FALSE
                        ) {
  # Get permutation of parameters
  if ("k" %in% toOpt) k <- 1:k
  if ("lambda1" %in% toOpt) lambda1 <- 1:lambda1
  if ("lambda2" %in% toOpt) lambda2 <- 1:lambda2
  if ("lambda3" %in% toOpt) lambda3 <- 1:lambda3
  if ("first_two_plus_k_pcs" %in% toOpt) first_two_plus_k_pcs <- 0:first_two_plus_k_pcs
  
  combos <- expand.grid(lambda1=lambda1,
                        lambda2=lambda2,
                        lambda3=lambda3,
                        k=k,
                        pcs = first_two_plus_k_pcs)
  
  ## If fixing num PCs, subset to rows w max num PCs
  if (fixPCs) combos <- combos[which(combos$pcs == max(combos$pcs)),]
  
  ## Storage lists
  max_cis <- list()
  max_cis$index <- 0
  max_cis$n <- 0
  max_cis$eqtls <- data.frame()
  max_cis$param_comb <- data.frame()
  
  max_trans <- list()
  max_trans$index <- 0
  max_trans$n <- 0
  max_trans$eqtls <- data.frame()
  max_trans$param_comb <- data.frame()
  
  for (i in 1:length(iterateMeQTL_results)) {
    # cis
    if(dim(iterateMeQTL_results[[i]]$cis$eqtls)[1] > max_cis$n) {
      max_cis$index <- i
      max_cis$n <- dim(iterateMeQTL_results[[i]]$cis$eqtls)[1]
      max_cis$eqtls <- iterateMeQTL_results[[i]]$cis$eqtls
      max_cis$param_comb <- combos[i,]
    }
    
    # trans
    if(dim(iterateMeQTL_results[[i]]$trans$eqtls)[1] > max_trans$n) {
      max_trans$index <- i
      max_trans$n <- dim(iterateMeQTL_results[[i]]$trans$eqtls)[1]
      max_trans$eqtls <- iterateMeQTL_results[[i]]$trans$eqtls
      max_trans$param_comb <- combos[i,]
    }
  }
  if (max_cis$index == 0) {
    warning("No max_cis_index")
  }
  if (max_trans$index == 0) {
    warning("No max_trans_index")
  }
  
  # Convert from 2+k PCs to PCs
  max_cis$param_comb$pcs <- max_cis$param_comb$pcs + 2
  max_trans$param_comb$pcs <- max_trans$param_comb$pcs + 2
  
  # Output summary statistics
  end_statistics <- list("cis" = max_cis,
                         "trans" = max_trans)
  
  return(end_statistics)
}

##################################################################################################################
# Parallelized function that will use hcps, pcs, and mapply in order to maximize the number of significant eQTLs #
##################################################################################################################
# COV_file_name = "tempcovariates.csv"; GENOTYPE_file_name = genotype; EXP_file_name = exp_nc; gene_location_file_name = geneloc; snps_location_file_name =snploc
# first_two_plus_k_pcs=3; kRange=2; lambda1Range=2; lambda2Range=2; lambda3Range=2; Z=Z; Y=Y_nc; covariates=covariates
futureIterateMeQTL <- function(first_two_plus_k_pcs, # how many of first principal components to use
                               kRange, lambda1Range, lambda2Range, lambda3Range, # hcp variables to optimize
                               Z, Y,
                               covariates, # matrix of covariates
                               GENOTYPE_file_name, 
                               EXP_file_name, 
                               gene_location_file_name, 
                               snps_location_file_name) {
  # Get permutation of parameters
  combos <- expand.grid(lambda1=1:lambda1Range,
                        lambda2=1:lambda2Range,
                        lambda3=1:lambda3Range,
                        k=1:kRange,
                        pcs = 0:first_two_plus_k_pcs)
  plan(multisession, workers=3)
  # iterate MeQTL
  results <-  future.apply::future_mapply(FUN="MeQTLandHcp",
                                          lambda1=combos$lambda1,
                                          lambda2=combos$lambda2,
                                          lambda3=combos$lambda3,
                                          k=combos$k,
                                          first_two_plus_k_pcs=combos$pcs,
                                          MoreArgs = list(Z=Z,
                                                          Y=Y,
                                                          covariates = covariates,
                                                          GENOTYPE_file_name = GENOTYPE_file_name, 
                                                          EXP_file_name = EXP_file_name, 
                                                          gene_location_file_name = gene_location_file_name, 
                                                          snps_location_file_name = snps_location_file_name),
                                          SIMPLIFY = FALSE)
  
  return(results)
}
