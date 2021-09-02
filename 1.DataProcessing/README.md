
# Script: Purpose

### Formatting to [MatrixeQTL format](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html)

1. FormatCovariates.R: Extract covariate information from massive covariates file, subset to European ancestry individuals, read in and format Principle Components from plink, code resulting dataframe to [MatrixeQTL format](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/Covariates.txt), and output as file.

2. FormatExpression_GeneLoc.R: Manipulate expression profile data in the mRNA and miRNA folders into [MatrixEQTL format](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/GE.txt). Identify location of ENSEMBL annotated genes and output in [MatrixeQTL format](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/geneloc.txt). Involves using hg19 genome from biomaRt, converting the miRNA (hsa-mri annotations) and mRNA genes (entrez annotations) to ENSEMBL annotations, differentiate between coding and non-coding genes, and output these files.

3. FormatGenotype.R: Use .bim, .bed and .fam files in conjunction with bigSNPR R package to format the observed genotypes according to [MatrixeQTL format](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/SNP.txt). Currently, it formats genotype for CHR22 and the total genome with the intention of using the CHR22 as a sort of training set.

4. FormatSNPloc.R: Use .bim, .bed and .fam files split by chromosome in conjunction with bigSNPR R package to format the SNP locations according to [MatrixeQTL format](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/Sample_Data/snpsloc.txt). Currently, it formats the SNP location for CHR22 and the total genome with the intention of using the CHR22 as a sort of training set.

5. GetEDPatients.R: Finds which patients in datasets are of European descent (ED) before finding which patients data are present in the expression (both mRNA and miRNA), genotype, and covariate datasets. Outputs a list of these ED patients, *uniquePatients.csv*.

6. IterateMatrixeQTL.R: Gets HCPs (via Rhcpp), formats data correctly, and writes function to iterate over MatrixEQTL's *Matrix_eQTL_engine* function with varying numbers of HCPs and PCs to maximize the number of significant eQTLs.

7.  OldFormatGeneLoc.R: outdated script to format Gene locations

8. ParseEnsemblAnnotations.R: Extract desired informatio from a downloaded ENSEMBL annotated [hg.19 genome](http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/). Ultimately not necessary as biomaRt was eventually able to be installed on Hoffman2.
