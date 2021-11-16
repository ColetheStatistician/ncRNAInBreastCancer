require(data.table)
snps = fread(snp_dosage_file)
snplocs = fread(snp_location_File)

nc_exp = fread(nc_exp_file)
genelocs = fread(gene_location_file)
genelocs = subset(genelocs, geneid %in% nc_exp$gene_id)
genelocs = subset(genelocs, chr %in% c(1:22))
nc_exp = subset(nc_exp, gene_id %in% genelocs$geneid)

require(isoTWAS)
for (i in 1:nrow(genelocs)){
    
    print(paste0(i,' of ',nrow(genelocs)))
    gene = genelocs$geneid[i]
    if (!gene %in% a$gene_id){
        print(gene)
        pheno = as.numeric(nc_exp[nc_exp$gene_id == gene,-1])
        snp_cur = subset(snps,
                         SNP %in% subset(snplocs,
                                        chr == genelocs$chr[i] &
                                            pos >= genelocs$left[i] - 1e6 &
                                            pos <= genelocs$right[i] + 1e6)$id)
        snpMat = as.matrix(snp_cur[,-1])
        Y = as.matrix(pheno)
        Y.rep = as.matrix(pheno)
        R = 1
        id = colnames(snpMat)
        rownames(snpMat) = snp_cur$SNP
        
        seed = sample(1:1e6,1)
        if (nrow(snpMat) > 1 & length(unique(pheno)) >= 3){
            enet = univariate_elasticnet(X = t(snpMat),
                                         Y,
                                         Omega = 1,
                                         family = 'gaussian',
                                         scale = F,
                                         alpha = 0.5,
                                         nfolds = 5,
                                         verbose = F,
                                         par = F,
                                         n.cores = NULL,
                                         tx_names = c(gene),
                                         seed)[[1]]
            
            blup = univariate_blup(X = t(snpMat),
                                   Y,
                                   Omega = 1,
                                   scale = F,
                                   alpha = 0.5,
                                   nfolds = 5,
                                   verbose = F,
                                   par = F,
                                   n.cores = NULL,
                                   tx_names = c(gene),
                                   seed)[[1]]
            
            susie = univariate_susie(X = t(snpMat),
                                     Y,
                                     Omega = 1,
                                     scale = F,
                                     alpha = 0.5,
                                     nfolds = 5,
                                     verbose = F,
                                     par = F,
                                     n.cores = NULL,
                                     tx_names = c(gene),
                                     seed)[[1]]
            
            lasso = univariate_elasticnet(X = t(snpMat),
                                          Y,
                                          Omega = 1,
                                          family = 'gaussian',
                                          scale = F,
                                          alpha = 1,
                                          nfolds = 5,
                                          verbose = F,
                                          par = F,
                                          n.cores = NULL,
                                          tx_names = c(gene),
                                          seed)[[1]]
            
            r2vec = unlist(c(enet$R2,blup$R2,susie$R2,lasso$R2))
            mmm = which.max(r2vec)
            if (max(r2vec) >= 0.005){
            if (mmm == 1){
                row = c(gene,scale(c(enet$Pred)))
            } else if (mmm == 2){
                row = c(gene,scale(c(blup$Pred)))
            } else if (mmm == 3){
                row = c(gene,scale(c(susie$Pred)))
            } else {
                row = c(gene,scale(c(lasso$Pred)))
            }
            row = as.matrix(row)
            row = as.data.frame(t(row))
            colnames(row) = colnames(nc_exp)
            fwrite(row,
                   nc_grex_file,
                   append = T,
                   row.names=F,
                   sep='\t',
                   quote=F)}
        }
    }
    
}

library(MatrixEQTL)
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = nc_grex_file;
snps_location_file_name = gene_location_file;

# Gene expression file name
expression_file_name = pc_exp_file;
gene_location_file_name = gene_location_file;

# Covariates file name
# Set to character() for no covariates
covariates_file_name = covariate_file;

# Output file name
output_file_name_cis = gbat_cis_file;
output_file_name_tra = gbat_tra_file;

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1;

errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 3e6;

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
snpspos = snpspos[,1:3]
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

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


require(susieR)
setwd('/u/scratch/a/abtbhatt/ncRNA/gtex_breast')
rm(list=ls())
require(data.table)
snps = fread(snp_dosage_file)
snplocs = fread(snp_location_file)


nc_exp = fread(nc_exp_file)
genelocs = fread(gene_location_file)
pc_exp = fread(pc_exp_file)

cov = fread(covariates_file)
pc_exp = as.data.frame(pc_exp)
pc_exp[,-1] = limma::removeBatchEffect(as.matrix(pc_exp[,-1]),
                                       covariates = t(cov[,-1]))

trans = fread(gbat_tra_file)
trans$FDR = ifelse(p.adjust(trans$`p-value`,method='BH') <= trans$FDR,
                   p.adjust(trans$`p-value`,method='BH'),
                   trans$FDR)
trans = subset(trans,FDR < 0.05)

require(biomaRt)
ensembl <- useEnsembl(biomart = "genes",
                          dataset = "hsapiens_gene_ensembl")
bm = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
           filters = c('ensembl_gene_id'),
           values = unique(c(trans$SNP,trans$gene)),
           mart = ensembl)

for (i in 1:nrow(trans)){
    
    nc_gene = trans$SNP[i]
    pc_gene = trans$gene[i]
    
    pheno = as.numeric(pc_exp[pc_exp$gene_id == pc_gene,-1])
    if (length(unique(pheno)) >= 5 & !all(is.na(pheno))){
        
        snp_cur = subset(snps,
                         SNP %in% subset(snplocs,
                                        chr == genelocs$chr[genelocs$geneid == nc_gene] &
                                            pos >= genelocs$left[genelocs$geneid == nc_gene] - 1e6 &
                                            pos <= genelocs$right[genelocs$geneid == nc_gene] + 1e6)$id)
        X = as.matrix(snp_cur[,-1])
        if (nrow(X) > 2){
            fitted <- susie(Matrix::Matrix(t(X),sparse=T),
                            pheno,
                            estimate_residual_variance =
                                TRUE,
                            estimate_prior_variance =
                                FALSE,
                            scaled_prior_variance = 0.1,
                            verbose = T,
                            intercept = F)
            
            sets <- susie_get_cs(fitted,
                                 X = t(X),
                                 coverage = 0.9,
                                 min_abs_corr = 0.1)
            cs_snps = snp_cur$SNP[sets$cs[[which.max(sets$coverage)]]]
            if (length(cs_snps) != 0){
                pips = fitted$pip[sets$cs[[which.max(sets$coverage)]]]
                
                locs = subset(snplocs,id %in% cs_snps)
                locs = locs[match(cs_snps,locs$id),]
                
                
                df = data.frame(ncRNA = nc_gene,
                                ncRNA_hgnc = NA,
                                pcGene = pc_gene,
                                pcGene_hgnc = NA,
                                pcGene_Chr = genelocs$chr[genelocs$geneid == pc_gene][1],
                                pcGene_Start = genelocs$left[genelocs$geneid == pc_gene][1],
                                pcGene_End = genelocs$right[genelocs$geneid == pc_gene][1],
                                Effect = trans$beta[i],
                                `t-stat` = trans$`t-stat`[i],
                                P = trans$`p-value`[i],
                                FDR = trans$FDR[i],
                                CausalSNP = cs_snps,
                                PIP = pips,
                                SNP_Chr = genelocs$chr[genelocs$geneid == nc_gene][1],
                                SNP_Pos = locs$pos)
            }
        } else {
            cs_snps =  subset(snplocs,
                              chr == genelocs$chr[genelocs$geneid == nc_gene] &
                                  pos >= genelocs$left[genelocs$geneid == nc_gene] - 1e6 &
                                  pos <= genelocs$right[genelocs$geneid == nc_gene] + 1e6)$id
            locs = subset(snplocs,id %in% cs_snps)
            locs = locs[match(cs_snps,locs$id),]
            
            df = data.frame(ncRNA = nc_gene,
                            ncRNA_hgnc = NA,
                            pcGene = pc_gene,
                            pcGene_hgnc = NA,
                            pcGene_Chr = genelocs$chr[genelocs$geneid == pc_gene][1],
                            pcGene_Start = genelocs$left[genelocs$geneid == pc_gene][1],
                            pcGene_End = genelocs$right[genelocs$geneid == pc_gene][1],
                            Effect = trans$beta[i],
                            `t-stat` = trans$`t-stat`[i],
                            P = trans$`p-value`[i],
                            FDR = trans$FDR[i],
                            CausalSNP = locs$snp,
                            PIP = NA,
                            SNP_Chr = genelocs$chr[genelocs$geneid == nc_gene][1],
                            SNP_Pos = locs$pos)
        }
        
        fwrite(df,
               gbat_finemap_file,
               append = T,
               row.names=F,
               sep='\t',
               quote=F)
        
        
    }
    
}
