# Script: Purpose

1. finalNC.R: Use optimized parameter values and total genome to get cis- and trans-eQTLs for non-coding RNA.
2. finalPC.R: Use optimized parameter values and total genome to get cis- and trans-eQTLs for coding genes.
3. prelimNC.R: Get non-optimized example results for BIG summer presentation for non-coding RNAs.
4. prelimPC.R: Get non-optimized example results for BIG summer presentation for coding genes.

# Notable Results

## non-coding RNA

| eQTL | p-value Threshold | Principal Components (PCs) | HCPs (k) | lambda1 | lambda2 | lambda3 | Num. Sig. |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| cis | 1e-6 | 4 | 1 | 1 | 1 | 1 | 1,298 |

## coding genes

| eQTL | p-value Threshold | Principal Components (PCs) | HCPs (k) | lambda1 | lambda2 | lambda3 | Num. Sig. |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| trans | 1e-6| 3 | 8 | 2 | 1 | 1 | 3,719,308 |
