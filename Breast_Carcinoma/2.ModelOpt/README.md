# Script: Purpose

1. IterateMatrixeQTLFunctions.R: Functions used to prepare data for and iterate MatrixEQTL
2. IterateMatrixeQTL_AB.R: Debugging help from Arjun
3. IterateNCMatrixeQTL.R: Iterate MatrixEQTL to find which combination of PCs and HCPs optimize the number of cis and trans eQTLs for non-coding RNA.
4. IteratePCMatrixeQTL.R: Iterate MatrixEQTL to find which combination of PCs and HCPs optimize the number of cis and trans eQTLs for coding genes.
5. optimize_nc_lambdas.R: Script that fixes number of PCs to 4 and varies other parameters to find the maximum number of significant eQTLs for non-coding RNA. 
6. optimize_pc_lambdas.R: Script that fixes number of PCs to 4 and varies other parameters to find the maximum number of significant eQTLs for coding genes. 

# Iterate Results

## iterate non-coding RNA

Using the `optimize_nc_lambdas.R` with a fixed 4 PCs and CHR22 data, I maximized the number of cis-eQTLs (3 total) and trans-eQTLs (10397 total). For both cis and trans, the maximizers were PCs=4 and k = lambda1 = lambda2 = lambda3 = 1. Using `ncIterateMatrixeQTL.R` and CHR22 data to vary the principal components, I found the number of principal components did not change the number of significant cis- and trans-eQTLs. Preliminary results:

| eQTL | Principal Components (PCs) | HCPs (k) | lambda1 | lambda2 | lambda3 | Num. Sig. |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| cis | 4 | 1 | 1 | 1 | 1 | 1,263 |
| trans | 4 | 1 | 1 | 1 | 1 | 670,462 |

## iterate coding genes

Using the `optimize_pc_lambdas.R` with a fixed 4 PCs and CHR22 data, I maximized the number of cis-eQTLs (68 total) and trans-eQTLs (63,998 total). For cis-eQTLs, the maximizers were PCs=4, k = 3. and lambda1 = lambda2 = lambda3 = 1. For trans-eQTLs, the maximizers were PCs=4, k=7, and lambda1 = lambda2 = lambda3 = 1. Using `pcIterateMatrixe_cis_QTL.R` with CHR22 data to vary the principal components, I found the number of principal components did not change the number of significant cis- and trans-eQTLs. For `pcIterateMatrixe_trans_QTL.R`, PCs=6 optmized the trans-eQTLs. Preliminary results:

| eQTL | Principal Components (PCs) | HCPs (k) | lambda1 | lambda2 | lambda3 | Num. Sig. |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| cis | 4 | 3 | 1 | 1 | 1 | X |
| trans | 4 | 7 | 1 | 2 | 1 | 142,595 |
