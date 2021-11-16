require(data.table)
require(glmnet)
require(susieR)

snps = fread(snp_tcga)
colnames(snps)[1] = 'id'
snplocs = fread(snploc_tcga)
colnames(snplocs)[1] = 'snp'
exp = fread(exp_tcga)
covs = fread(cov_tcga)
genelocs = fread(geneloc_tcga)
colnames(genelocs) = c('geneid','chr','left','right')

exp = as.data.frame(exp)
colnames(exp)[1] = 'Gene'
exp[,-1] = limma::removeBatchEffect(as.matrix(exp[,-1]),
                                    covariates = t(covs[,-1]))
genelocs = subset(genelocs,
                  geneid %in% exp$Gene)
exp = subset(exp,Gene %in% genelocs$geneid)
exp = subset(exp,Gene %in% rep_lincs)

require(doMC)
registerDoMC(cores = 5)

for (i in 1:nrow(exp)){
  print(i)
  w = which(genelocs$geneid == exp$Gene[i])
  snp.cur = subset(snps,id %in% subset(snplocs,
                                       chr == genelocs$chr[w] &
                                         pos <= genelocs$right[w] + 1e6 &
                                         pos >= genelocs$left[w] - 1e6)$snp)
  
  pheno = as.numeric(scale(as.numeric(exp[i,-1])))
  snpMat = as.matrix(snp.cur[,-1])
  
  if (nrow(snpMat) > 0){
    elastic.net = cv.glmnet(x = t(snpMat),
                            y = pheno,
                            nfolds = 5,
                            gamma = c(0,.25,.5,.75,1),
                            family = 'gaussian',
                            keep = T,
                            parallel = T,trace.it=1)
    
    r2 = cor(pheno,as.numeric(elastic.net$fit.preval[,which.min(elastic.net$cvm)]))^2
    
    if (r2 > .01){
      print(r2)
      Model = data.frame(Gene = exp$Gene[i],
                         SNP = snp.cur$id,
                         Chromosome = snplocs[match(snp.cur$id,
                                                    snplocs$snp),]$chr,
                         Position = snplocs[match(snp.cur$id,
                                                  snplocs$snp),]$pos,
                         Weight = coef(elastic.net,s = 'lambda.min')[-1])
      Model = subset(Model,Weight != 0)
      if (nrow(Model) > 0){
        
        print(exp$Gene[i])
        fwrite(Model,model_file,
               sep='\t',append=T,row.names=F)
        
      }
      
    }
  }
  
}
