# PBMC_covariation
Code repository for reproducing analyses in [Khoury et al., 2026](https://doi.org/10.64898/2026.09.11.751068), [PDF](https://slavovlab.net/Slavov-Lab-Publications/2026_cell-type-specific-functional-coordination-in-PBMCs.pdf), [Website](https://scp.slavovlab.net/Khoury_et_al_2026)

This codebase covers all raw data processing, outputs, analysis, and figure generation.

Raw sequencing and LC-MS/MS data will be made available.

All meta data, processed data, and files required for reproducing all analyses in this repository are accessible via Zenodo DOI: [10.5281/zenodo.22649483](https://zenodo.org/records/22649483)

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
10. PBMC_Covariation.Rmd (scripts 1 through 8)

Meta data, raw data, processed abundance matrices, and .rds files are loaded in each script where needed.

Script 4 in PBMC_Covariation.Rmd is the only script within the markdown file that needs each protein preprocessing script run first, but otherwise can be skipped to proceed with scripts 5-8.

Executing scripts sequentially in PBMC_Covariation.Rmd (main analysis) requires first running:
1. Protein_data_integration_annotation.R
2. SS3xpress_preprocessing_and_QC.R

