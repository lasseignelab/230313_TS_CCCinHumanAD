README
================
2023-11-27

## Data Directory Structure

Files in this directory are either generated using code from this
project or downloaded from sources specified in our scripts. The
contents of this directory are also deposited on zenodo. Details
(incl. DOIs) can be found in the main repository’s README.

The data directory should include the following files:

    ## .
    ## +-- CellRangerCounts
    ## |   +-- GSE157827
    ## |   |   +-- post_soupX
    ## |   |   |   +-- SAMN16100276_S10_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100278_S09_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100280_S08_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100282_S07_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100284_S06_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100286_S05_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100289_S04_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100290_S01_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100291_S03_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100292_S02_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100293_S18_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100294_S17_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100295_S16_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100296_S15_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100297_S14_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100299_S13_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100302_S12_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100304_S11_AD
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100306_S21_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   +-- SAMN16100307_S20_CTRL
    ## |   |   |   |   +-- barcodes.tsv
    ## |   |   |   |   +-- genes.tsv
    ## |   |   |   |   \-- matrix.mtx
    ## |   |   |   \-- SAMN16100308_S19_CTRL
    ## |   |   |       +-- barcodes.tsv
    ## |   |   |       +-- genes.tsv
    ## |   |   |       \-- matrix.mtx
    ## |   |   \-- pre_soupX
    ## |   |       +-- SAMN16100276_S10_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100278_S09_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100280_S08_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100282_S07_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100284_S06_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100286_S05_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100289_S04_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100290_S01_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100291_S03_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100292_S02_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100293_S18_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100294_S17_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100295_S16_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100296_S15_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100297_S14_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100299_S13_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100302_S12_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100304_S11_AD
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100306_S21_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       +-- SAMN16100307_S20_CTRL
    ## |   |       |   +-- filtered_feature_bc_matrix.h5
    ## |   |       |   \-- raw_feature_bc_matrix.h5
    ## |   |       \-- SAMN16100308_S19_CTRL
    ## |   |           +-- filtered_feature_bc_matrix.h5
    ## |   |           \-- raw_feature_bc_matrix.h5
    ## |   \-- GSE174367
    ## |       +-- SAMN19128593_S19_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128594_S18_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128595_S17_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128596_S16_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128597_S15_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128598_S14_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128599_S13_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128600_S12_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128601_S11_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128602_S10_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128603_S9_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128604_S8_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128605_S6_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128606_S5_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128607_S4_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128608_S3_AD
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128609_S2_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       +-- SAMN19128610_S1_CTRL
    ## |       |   +-- barcodes.tsv.gz
    ## |       |   +-- features.tsv.gz
    ## |       |   \-- matrix.mtx.gz
    ## |       \-- SAMN19128611_S7_AD
    ## |           +-- barcodes.tsv.gz
    ## |           +-- features.tsv.gz
    ## |           \-- matrix.mtx.gz
    ## +-- README.Rmd
    ## +-- README.md
    ## +-- ccc
    ## |   +-- geo_multinichenet_output.rds
    ## |   +-- geo_signaling_igraph_objects.rds
    ## |   +-- gse_multinichenet_output.rds
    ## |   +-- gse_signaling_igraph_objects.rds
    ## |   +-- nichenet_grn
    ## |   |   +-- gr_network_human_21122021.rds
    ## |   |   +-- ligand_tf_matrix_nsga2r_final.rds
    ## |   |   +-- signaling_network_human_21122021.rds
    ## |   |   \-- weighted_networks_nsga2r_final.rds
    ## |   +-- nichenet_prior
    ## |   |   +-- ligand_target_matrix.rds
    ## |   |   \-- lr_network.rds
    ## |   \-- nichenet_v2_prior
    ## |       +-- ligand_target_matrix_nsga2r_final.rds
    ## |       \-- lr_network_human_21122021.rds
    ## \-- seurat_preprocessing
    ##     +-- geo_clustered_seurat.rds
    ##     +-- geo_filtered_seurat.rds
    ##     +-- geo_integrated_seurat.rds
    ##     +-- geo_processed_seurat.rds
    ##     +-- gse_clustered_seurat.rds
    ##     +-- gse_filtered_seurat.rds
    ##     +-- gse_integrated_seurat.rds
    ##     \-- gse_processed_seurat.rds
