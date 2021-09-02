# Script :Purpose
1. MakeTriplets.R: Script that takes the resulting **.rds** files from the final eQTL models to find triplets (SNP, pc-gene, nc-RNA). It also finds which SNP/pc-gene combinations have single/multiple nc-RNA(s) associated with them.
2. MediAnalysisFunctions.R: Functions used to perform mediation analysis on triplets with functions borrowed from MOSTWAS.
3. SinglesMediationAnalysis.R: Performs mediation analysis on triplets with a single nc-RNA.
4. formatCovariatesforMA.R: Formats HCPs, PCs, and continuous covariates in preparation to regress them out in the mediation analysis step.
5. multiNCMediationAnalysis.R: Performs mediation analysis on triplets with a multiple nc-RNA.
6. subsetMediationAnalysis.R: Test script
