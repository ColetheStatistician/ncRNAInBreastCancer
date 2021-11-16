require(bigsnpr)
require(data.table)

calc_twas = function(w,
                     z,
                     ld,
                     indices,
                     return = 'z'){
  w = w[indices]
  eff = (w %*% z)
  se = (w %*% ld %*% w)
  if (return == 'z'){
    return( eff/se )
  }
  if (return == 'all'){
    return(list(eff = eff,
                se = se))
  }
}

for (g in unique(models$Gene)){
  print(g)
  models.cur = subset(models,Gene == g)
  sumStats.cur = subset(sumStats,chr.iCOGs == models.cur$Chromosome[1] &
                          Position.iCOGs %in% models.cur$Position)
  models.cur = subset(models.cur,Position %in% sumStats.cur$Position.iCOGs)
  models.cur = models.cur[match(sumStats.cur$Position.iCOGs,
                                models.cur$Position),]
  ld.cur = subset(ld,id %in% models.cur$SNP)
  ld.cur = ld.cur[match(models.cur$SNP,
                        ld.cur$id),]
  ld.mat = tcrossprod(as.matrix(ld.cur[,-1]))/(ncol(ld.cur)-2)
  sumStats.cur$Z = sumStats.cur$Beta.meta/sumStats.cur$sdE.meta
  if (length(models.cur$Weight) > 0){
    Z.all = boot::boot(data = as.numeric(models.cur$Weight),
                       statistic = calc_twas,
                       R = 10000,sim='permutation',
                       z = as.numeric(sumStats.cur$Z),
                       ld = as.matrix(ld.mat))
    Z = as.numeric(Z.all$t0)
    P = 2*pnorm(-abs(Z))
    perm.P = (sum(abs(as.numeric(Z.all$t)) >= abs(as.numeric(Z.all$t0))) + 1)/10001
    data.table::fwrite(data.frame(Gene = g,
                                  Z = Z,
                                  P = P,
                                  Permutation.P = perm.P,
                                  Tissue = 'Breast',
                                  State = 'Tumor',
                                  Trait = 'Overall BRCA Risk'),
                       grex_results_file,
                       append=T,
                       row.names=F,
                       sep='\t')}
  
  
}
