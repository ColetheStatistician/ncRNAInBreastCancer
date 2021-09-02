# Solution from https://speciationgenomics.github.io/pca/
library(data.table)
library(ggplot2)

# Read in data
pca  <- data.table::fread("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/plink.eigenvec")
eigenval <- scan("/u/scratch/n/nolanc/Breast_Carcinoma/MeQTL/plink.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

# plot
pdf("pca.pdf")
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percent variance explained") + theme_light() + theme(aspect.ratio = 1, axis.text = element_text(size = 15), axis.title = element_text(size = 20)) 
dev.off()
