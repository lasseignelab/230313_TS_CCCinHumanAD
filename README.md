README
================
2023-11-28

# Altered Glia-Neuron Communication in Alzheimer’s Disease Affects WNT, p53, and NFkB Signaling Determined by snRNA-seq

## Authors

**Tabea M. Soelter, Timothy C. Howton, Amanda D. Clark, Vishal H. Oza,
Brittany N. Lasseigne**

The University of Alabama at Birmingham, Heersink School of Medicine,
Department of Cell, Developmental and Integrative Biology

## Project Overview

We used publicly available snRNA-seq AD data (Morabito et al., 2021)
generated from postmortem human PFC to study altered glia-neuron
interactions and their downstream effects in AD. We inferred
differential CCC interactions between astrocytes, microglia,
oligodendrocytes, or OPCs (sender cell types) and inhibitory or
excitatory neurons (receiver cell types). We also investigated whether
CCC shows high similarity across cell types by calculating the Jaccard
Similarity Index (JI) of ligands, receptors, and target genes across
cell types. Since CCC inference methodologies are known to produce false
positives, we validated our interactions using an independent human PFC
AD snRNA-seq dataset (Lau et al., 2020). We further investigated the
resulting high-confidence ligand-receptor pairs from both data sets,
their predicted downstream target genes, and signaling modulators
through transcription factor (TF) and canonical signaling pathway
activity.

![alt
text](https://github.com/lasseignelab/230313_TS_CCCinHumanAD/blob/main/results/figures/figure1.png)

## Datasets

The datasets used in this study can be found on GEO:  
1. [Morabito et al.,
2021](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367)  
2. [Lau et al.,
2020](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157827)

## Scripts

Data Alignment (Cell Ranger):

    ## src/cellranger/
    ## +-- geo
    ## |   +-- SAMN19128593_S19_CTRL.sh
    ## |   +-- SAMN19128594_S18_CTRL.sh
    ## |   +-- SAMN19128595_S17_CTRL.sh
    ## |   +-- SAMN19128596_S16_CTRL.sh
    ## |   +-- SAMN19128597_S15_CTRL.sh
    ## |   +-- SAMN19128598_S14_CTRL.sh
    ## |   +-- SAMN19128599_S13_AD.sh
    ## |   +-- SAMN19128600_S12_AD.sh
    ## |   +-- SAMN19128601_S11_AD.sh
    ## |   +-- SAMN19128602_S10_AD.sh
    ## |   +-- SAMN19128603_S9_AD.sh
    ## |   +-- SAMN19128604_S8_AD.sh
    ## |   +-- SAMN19128605_S6_AD.sh
    ## |   +-- SAMN19128606_S5_AD.sh
    ## |   +-- SAMN19128607_S4_AD.sh
    ## |   +-- SAMN19128608_S3_AD.sh
    ## |   +-- SAMN19128609_S2_CTRL.sh
    ## |   +-- SAMN19128610_S1_CTRL.sh
    ## |   \-- SAMN19128611_S7_AD.sh
    ## \-- gse
    ##     +-- GSE157827_all_array.sh
    ##     +-- id_sheet.csv
    ##     \-- sample_sheet.csv

Ambient RNA removal (SoupX):

    ## src/soupX/
    ## \-- 01_ambientRNA_removal.Rmd

Pre-processing (Seurat):

    ## src/seurat_preprocessing/
    ## +-- 01_geo_seurat_preprocessing.Rmd
    ## \-- 02_gse_seurat_preprocessing.Rmd

Cell-cell communication inference (MultiNicheNet) and JI calculations:

    ## src/ccc/
    ## +-- 01_differential_ccc_multinichenet.Rmd
    ## +-- 02_jaccard_index_geo.Rmd
    ## +-- 03_jaccard_index_gse.Rmd
    ## \-- 04_gene_regulatory_networks.Rmd

Biological Activity (decoupleR):

    ## src/biological_activity/
    ## +-- 01_pseudobulk_dea.Rmd
    ## +-- 02_tf_activity.Rmd
    ## \-- 03_pathway_activity.Rmd

Manuscript figures:

    ## src/manuscript_figures/
    ## +-- figure_2.Rmd
    ## +-- figure_3.Rmd
    ## +-- figure_4.Rmd
    ## +-- figure_5.Rmd
    ## +-- figure_S1.Rmd
    ## +-- figure_S2.Rmd
    ## +-- figure_S3.Rmd
    ## \-- figure_S4.Rmd

## Code and Data Availability

Data:
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.10214497-blue)](https://doi.org/10.5281/zenodo.10214497)  
Docker Images:
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.10214660-blue)](https://doi.org/10.5281/zenodo.10214660)  
GitHub Repository:
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.10211622-blue)](https://doi.org/10.5281/zenodo.10211622)

## Docker

We performed all analyses in docker with R version 4.1.3. The repository
with all the Docker images used for this project can be found on Docker
Hub at
[tsoelter/rstudio_ccc_ad](https://hub.docker.com/repository/docker/tsoelter/rstudio_ccc_ad/general).
While individual docker image tags are noted in every script, we provide
an overview of analyses and their associated tags below:

| Tag   | Associated Analyses                                          |
|-------|--------------------------------------------------------------|
| 1.0.1 | Pre-processing of Morabito et al., 2021                      |
| 1.0.3 | Pre-processing of Lau et al., 2020, CCC, and JI              |
| 1.0.5 | GRNs, pseudo-bulking, DEA, biological activity, and plotting |

## Lasseigne Lab

[What is Happening in the Lasseigne Lab?](https://www.lasseigne.org/)

<img src="https://www.lasseigne.org/img/main/lablogo.png" width="75" height="75">

## Funding

This work was supported in part by the UAB Lasseigne Lab funds, the NIA
R00HG009678-04S1, the Alzheimer’s of Central Alabama Lindy Harrell
Predoctoral Scholar Program.

## Acknowledgements

We would also like to thank the members of the Lasseigne Lab,
specifically Jordan H. Whitlock, Emma F. Jones, and Elizabeth J. Wilk
for their valuable input throughout this study.

## License

[![MIT
License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)
