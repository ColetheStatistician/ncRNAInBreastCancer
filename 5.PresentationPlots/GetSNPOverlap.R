# Load libraries
library(ggplot2)
library(ggbreak) 
library(knitr)

# read in results
optPC <- readRDS("pc_final_trans_results.rds")
optNC <- readRDS("nc_final_results.rds")

# Get list of unique SNPs in intersection between the pc cis-eQTLs and nc trans-eQTLs
uniq_snps <- intersect(optPC$eQTLs$trans$eqtls$snps, optNC$eQTLs$cis$eqtls$snps)

# barchart of uniques
df1 <- data.frame("number" = c(length(setdiff(optNC$eQTLs$cis$eqtls$snps, optPC$eQTLs$trans$eqtls$snps)), 
                               length(uniq_snps), 
                               length(setdiff(optPC$eQTLs$trans$eqtls$snps, optNC$eQTLs$cis$eqtls$snps))),
                 "category" = c("NC cis", "Both", "PC trans"))
pdf("UniqueSNPs.pdf")
ggplot(data=df1, aes(x = category, y = number)) + geom_bar(stat="identity") + scale_y_break(c(4e2, 7.984e5)) + xlab("SNPs in significant eQTLs") + ylab("Number of unique SNPs in eQTLs")
dev.off()

# Table of top p value trans qtls with overlapping cis
trans_w_overlapping_cis <- optPC$eQTLs$trans$eqtls[which(optPC$eQTLs$trans$eqtls$snps %in% uniq_snps),]
rownames(trans_w_overlapping_cis) <- NULL
# library(kableExtra)
topresults <- head(trans_w_overlapping_cis[,-c(3,6)], n=10)
colnames(topresults) <- c("SNP", "ENSEMBL Gene", "p-value", "FDR")
kable(topresults)
# kable(topresults, format = "html", row.names = TRUE) %>%
#  kable_styling(bootstrap_options = c("striped", "hover"),
#                full_width = F,
#                font_size = 12,
#                position = "left")
