# Load libraries
library(data.table)
library(Rhcpp)
library(MatrixEQTL)

##########################
# Read in formatted data #
##########################
genotype <- "22.Breast_Carcinoma_genotype.csv"
snploc <- "22.Breast_Carcinoma_snploc.csv"
geneloc <- "formatted_gene_loc.csv"
exp_nc <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv"
exp_pc <- "/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv"
# Read in expression data
covariates <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formattedCovariates.csv")
Y_nc <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_nc_expression.csv")
Y_pc <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/formatted_pc_expression.csv")

########################
# Format data for hcps #
########################
# ?hcp
# Make Z Matrix: nxd where n=number of subjects and d=the number of known covariates
Z <- data.frame(t(covariates))
colnames(Z) <- Z[1,]
Z <- Z[-1,]
Z[,1] <- as.numeric(Z[,1])
## Standardize age
Z$patient.age_at_initial_pathologic_diagnosis = as.numeric(Z$patient.age_at_initial_pathologic_diagnosis)
Z$patient.age_at_initial_pathologic_diagnosis <- (Z$patient.age_at_initial_pathologic_diagnosis - mean(Z$patient.age_at_initial_pathologic_diagnosis)) / sd(Z$patient.age_at_initial_pathologic_diagnosis)

# Make Y Matrix: nxg where n=number of subjects and g are genes
## Make unique rownames 
nc_rows <- make.names(as.matrix(Y_nc[,1]), unique = TRUE) 
pc_rows <- make.names(as.matrix(Y_pc[,1]), unique = TRUE) 

Y_nc <- data.frame(Y_nc[,-1])
Y_pc <- data.frame(Y_pc[,-1])

rownames(Y_nc) <- nc_rows
rownames(Y_pc) <- pc_rows

## Make nxg
Y_nc <- t(Y_nc)
Y_pc <- t(Y_pc)
Y_nc = Y_nc[,which(apply(Y_nc,2,sd) != 0)]

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

# k=5; lambda1=5; lambda2=5; lambda3=5; Z=Z; Y=Y_nc; first_k_pcs=4
getHCPs <- function(first_k_pcs, # how many of first principal components to use
                    k, lambda1, lambda2, lambda3, Z, Y # hcp variables to optimize
) {
    # Get PCs
    chosenPCs <- make.names(rep("PC", first_k_pcs+1), unique = TRUE)[-1]
    colsToGet <- c(colnames(Z)[(length(colnames(Z))-4):length(colnames(Z))], chosenPCs)
    Z <- data.matrix(Z[,colsToGet])
    # Get HCPs
    HCPs <- r_hcp(Z = scale(Z), Y = scale(Y),
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
    pvOutputThreshold_cis = 2e-2;
    pvOutputThreshold_tra = 1e-2;
    
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

#####################################################################################################
# Function that will use hcps, pcs, and mapply in order to maximize the number of significant eQTLs #
#####################################################################################################
# k=5; lambda1=5; lambda2=5; lambda3=5; Z=Z; Y=Y_nc; first_k_pcs=4
iterateMeQTL <- function(first_k_pcs, # how many of first principal components to use
                         k, lambda1, lambda2, lambda3, # hcp variables to optimize
                         Z, Y,
                         covariates,  
                         GENOTYPE_file_name, 
                         EXP_file_name, 
                         gene_location_file_name, 
                         snps_location_file_name) {
    
    # Get HCPs
    HCPs <- getHCPs(first_k_pcs = first_k_pcs, k = k, lambda1 = lambda1, lambda2 = lambda2, lambda3 = lambda3, Z = Z, Y = Y)
    
    # Format HCPS and make covariates
    covariates <- formatHCPstoCov(hcps = HCPs, covariates = covariates)
    # save this new covariates matrix in a new file
    filename <- "tempcovariates.csv"
    data.table::fwrite(covariates, filename)
    
    # Perform Matrix eQTL
    me <- MeQTL(COV_file_name = filename, 
                GENOTYPE_file_name = GENOTYPE_file_name, 
                EXP_file_name = EXP_file_name, 
                gene_location_file_name = gene_location_file_name, 
                snps_location_file_name = snps_location_file_name)
    
    ## Results:
    
    cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
    cat('Detected local eQTLs:', '\n');
    show(me$cis$eqtls)
    cat('Detected distant eQTLs:', '\n');
    show(me$trans$eqtls)
    
    ## Plot the Q-Q plot of local and distant p-values
    
    plot(me)
}
