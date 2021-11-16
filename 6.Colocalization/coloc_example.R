require(vroom)
require(MatrixEQTL)
require(data.table)
intersect = as.data.frame(fread(gwas_intersect_file,sep = '\t'))


gwas = vroom(gwas_sumstats_file)
gwas = gwas[,c('chr.iCOGs','Position.iCOGs',
               'Triple_Neg_log_or_meta',
               'Triple_Neg_se_meta',
               'EAFcases.iCOGs')]

colnames(gwas)[3:5] = c('Beta.meta',
                        'var.meta',
                        'MAF')


sdY.est <- function(vbeta, maf, n) {
  warning("estimating sdY from maf and varbeta, please directly supply sdY if known")
  oneover <- 1/vbeta
  nvx <- 2 * n * maf * (1-maf)
  m <- lm(nvx ~ oneover - 1)
  cf <- coef(m)[['oneover']]
  if(cf < 0)
    stop("estimated sdY is negative - this can happen with small datasets, or those with errors.  A reasonable estimate of sdY is required to continue.")
  return(sqrt(cf))
}

intersect.double = subset(intersect,
                          grepl(', ',ncRNA))
intersect = subset(intersect,
                   !grepl(', ',ncRNA))

for (i in 1:nrow(intersect.double)){
  
  nnn = unlist(strsplit(intersect.double$ncRNA[i],', '))
  df.now = intersect.double[i,]
  
  for (j in 1:length(nnn)){
    df.now$ncRNA = nnn[j]
    intersect = rbind(intersect,df.now)
  }
  
}

for (dataset in c('GTEx Breast','TCGA-BRCA')){
  print(dataset)
  intersect.cur = subset(intersect, Dataset == dataset)
  intersect.cur$Pair = paste(intersect.cur$ncRNA,
                             intersect.cur$pcGene,sep=':')
  intersect.cur = intersect.cur[!duplicated(intersect.cur$Pair),]
  
  if (dataset == 'GTEx Breast'){
    
    snpfile = snp_gtex
    expfile = exp_gtex
    covfile = cov_gtex
    snps_location_file_name = snploc_gtex
    gene_location_file_name = geneloc_gtex
    
  } else {
    snpfile = snp_tcga
    expfile = exp_tcga
    covfile = cov_tcga
    snps_location_file_name = snploc_tcga
    gene_location_file_name = geneloc_tcga
  }
  intersect.cur = intersect.cur[intersect.cur$TNBCRisk != '',]
  out.list = list()
  require(rlist)
  snp_full = data.table::fread(snpfile)
  for (i in 1:nrow(intersect.cur)){
    print(paste0(i,' out of ',nrow(intersect.cur)))
    ncRNA = unlist(strsplit(intersect.cur$ncRNA[i],', '))
    ncRNA = stringr::str_replace(ncRNA,' ','')
    pcGene = intersect.cur$pcGene[i]
    print(ncRNA)
    print(pcGene)
    sink()
    #### Extract SNPs and genes
    snppos = vroom(snps_location_file_name,progress=F,show_col_types=F)
    colnames(snppos) = c('id','chr','pos')
    genepos = vroom(gene_location_file_name,progress=F,show_col_types=F)
    colnames(genepos) = c('id','chr','left','right')
    genepos = subset(genepos, id %in% c(ncRNA,pcGene))
    snppos = subset(snppos, chr == genepos$chr[genepos$id == ncRNA] &
                      pos <= genepos$right[genepos$id == ncRNA] + 1e6 &
                      pos >= genepos$left[genepos$id == ncRNA] - 1e6)
    www = c()
    if (length(www) == 0){
      
      colnames(snp_full)[1] = 'SNP'
      snp_all = subset(snp_full,SNP %in% snppos$id)
      
      eee = fread(expfile)
      colnames(eee)[1] = 'id'
      eee = subset(eee,id %in% c(ncRNA,pcGene))
      cov = fread(covfile)
      
      snp_all = as.data.frame(snp_all)
      snp_all = snp_all[apply(snp_all[,-1],1,sd) != 0,]
      eee = as.data.frame(eee)
      eee = eee[match(c(ncRNA,pcGene),
                      eee$id),]
      
      Xc = t(as.matrix(cov[,-1]))
      y = t(eee[,-1])
      eee[,-1] = as.matrix(t(y - Xc %*% MASS::ginv(tcrossprod(t(Xc))) %*% t(Xc) %*% y))
      
      
      get_beta = function(x,y){
        
        reg = lm(y ~ x)
        return(as.numeric(coef(summary(reg))[2,1:2]))
        
      }
      
      one = apply(snp_all[,-1],1,get_beta,y = as.numeric(eee[1,-1]))
      two = apply(snp_all[,-1],1,get_beta,y = as.numeric(eee[2,-1]))
      
      
      gwas.cur = dplyr::filter(gwas, chr.iCOGs == genepos$chr[genepos$id == ncRNA] &
                                 Position.iCOGs <= genepos$right[genepos$id == ncRNA] + 1e6 &
                                 Position.iCOGs >= genepos$left[genepos$id == ncRNA] - 1e6)
      
      require(coloc)
      gwas.D1 = gwas.cur[,c('Beta.meta',
                            'var.meta',
                            'chr.iCOGs',
                            'Position.iCOGs')]
      gwas.D1$sdY = tryCatch(sdY.est(gwas.cur$var.meta^2,
                                     gwas.cur$MAF,
                                     2e5),error = function(e){2.5})
      gwas.D1$snp = paste(gwas.D1$chr.iCOGs,gwas.D1$Position.iCOGs,
                          sep=':')
      colnames(gwas.D1)[1:4] = c('beta',
                                 'varbeta',
                                 'type',
                                 'position')
      gwas.D1$varbeta = gwas.D1$varbeta^2
      gwas.D1$type = 'cc'
      
      snppos = subset(snppos, id %in% snp_all$SNP)
      snppos = snppos[match(snp_all$SNP,
                            snppos$id),]
      
      cis.eqtl = data.frame(beta = one[1,],
                            varbeta = one[2,]^2,
                            type = 'quant',
                            position = snppos$pos,
                            sdY = sd(eee[1,-1]),
                            snp = paste(snppos$chr,snppos$pos,sep=':'))
      
      tra.eqtl = data.frame(beta = two[1,],
                            varbeta = two[2,]^2,
                            type = 'quant',
                            position = snppos$pos,
                            sdY = sd(eee[2,-1]),
                            snp = paste(snppos$chr,snppos$pos,sep=':'))
      gwas.D1 = gwas.D1[!duplicated(gwas.D1$snp),]
      tra.eqtl = tra.eqtl[!duplicated(tra.eqtl$snp),]
      cis.eqtl = cis.eqtl[!duplicated(cis.eqtl$snp),]
        
        my.res.cis <- coloc::coloc.abf(dataset1=gwas.D1,
                                       dataset2=cis.eqtl)
        my.res.tra <- coloc::coloc.abf(dataset1=gwas.D1,
                                       dataset2=tra.eqtl) 
      
      out.df = list(Tissue = dataset,
                    Trait = 'TNBC BRCA risk',
                    ncRNA = ncRNA,
                    coloc.ncRNA = my.res.cis,
                    pcGene = pcGene,
                    coloc.pcGene = my.res.tra)
      out.list = list.append(out.list,out.df)
      
    }}
  saveRDS(out.list,paste0(ifelse(dataset == 'TCGA-BRCA',
                                 'TCGA','GTEx'),
                          '_CancerRisk_Colocalization'))
}
