# PBMC_covariation
Code repository for reproducing analyses in 2026_Khoury_et _al.

This codebase covers all raw data processing, outputs, analysis, and figure generation.

Raw data accessible via: 
Meta data, processed abundance matrices, and .rds files accessible via: 

When executing any script or markdown file, always run Install_and_load_packages.R and Functions_for_mapping_and_analysis.R first.
To reproduce the exact analysis, run:
1. Install_and_load_packages.R
2. Functions_for_mapping_and_analysis.R
3. n6p2_preprocessing_and_QC.R
4. n6p1_preprocessing_and_QC.R
5. n3p1_preprocessing_and_QC.R
6. n2p2_preprocessing_and_QC.R
7. n2p1_preprocessing_and_QC.R
8. Protein_data_integration_annotation.R
9. SS3xpress_preprocessing_and_QC.R
10. PBMC_Covariation.Rmd

Meta data, processed abundance matrices, and .rds files are loaded in each script where needed (within PBMC_Covariation.Rmd).
Script 4 in PBMC_Covariation.Rmd is the only script within the markdown file that needs each protein preprocessing script run first.
Executing all sections of PBMC_Covariation.Rmd (main analysis) requires first running:
1. Protein_data_integration_annotation.R
2. SS3xpress_preprocessing_and_QC.R
