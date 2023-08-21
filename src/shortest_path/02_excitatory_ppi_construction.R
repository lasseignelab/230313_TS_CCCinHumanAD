#### PPI network generation and edge weight calculation for Dijkstra's in Excitatory Neurons.
### Author: Tabea M. Soelter
### Date: 2023-05-16

## Goal: Create igraph objects of AD and CTRL PPI networks for multiple datasets

## time tracking
ptm <- proc.time()

## set seed
set.seed(42)
print("seed set")

## load packages
library(tidyverse)
library(STRINGdb)
library(igraph)
library(Seurat)
print("packages loaded")

## enable usage of args
args <- R.utils::commandArgs(trailingOnly = TRUE)
print("enabled args usage")

wd <- args[4]

## load my functions
source(paste0(wd, "/src/functions_CCCinAD_tms.R"))
print("functions loaded")

## get dataset name
name <- sub("_processed_seurat.rds.*", "", basename(args[1]))

## rename vine cortex for saving purposes
if (name == "vinecortex") {
  name <- "vine_cortex"
  print(name)
} else {
  print(name)
}

## load in data
# Seurat object
object <- readRDS(args[1])
print("loaded seurat object")

# mapped genes generated in 01_ppi_input_generation.Rmd
mapped_genes <- readRDS(args[2])
print("loaded mapped gene input")

# tmp ppi generated in 01_ppi_input_generation.Rmd
ppi_tmp <- readRDS(args[3])
print("loaded tmp ppi")

# AD-risk gene list compiled
ad_gene_list <- read.csv(paste0(
  wd, "/data/ccc/ad_gene_list.csv"),
  header = FALSE
)
print("loaded AD-risk gene list")

## Create igraph objects
# Wrapper function to:
#   - Create igraph object by condition of a Seurat object
# Calculate edge weights
#   - While STRINGdb provides us with scores which indicate the confidence in the interaction based on prior knowledge, we are interested in interactions directly tied to the gene expression (GEx) of our data.
#   - In order to account for GEx, we calculate edge weights based on the sum of the gene expression of 2 connected nodes (by an edge) and scale the value by the STRINGdb score to continue accounting for the prior knowledge.
#     - An annotated Seurat Object is used to obtain GEx values.
#   - Values are inversed to allow for shortest path calculation prioritizing genes with highest GEx
# Outputs a list of the condition specific igraph objects
object_list <- create_igraph_object(condition = c("AD", "CTRL"),
                                    receiver = "Excitatory Neurons",
                                    seurat_object = object,
                                    mapped_genes = mapped_genes,
                                    ppi_tmp = ppi_tmp)
print("created igraph objects")

# pull out individual objects
AD_igraph <- object_list[[1]]
CTRL_igraph <- object_list[[2]]
print("split objects")

# save objects
print("saving objects")
#check if data/shortest_path exists, create if not
if (!dir.exists(paste0(wd, "/data/shortest_path/"))) {
  dir.create(paste0(wd, "/data/shortest_path/")) 
} 
saveRDS(
  AD_igraph,
  file = paste0(wd, "/data/shortest_path/", name, "_AD_ex_igraph.rds")
)
saveRDS(
  CTRL_igraph,
  file = paste0(wd, "/data/shortest_path/", name, "_CTRL_ex_igraph.rds")
)

# session info
sessionInfo()

# time tracking
fptm <- proc.time() - ptm
fptm <- (fptm[3] / 60) / 60
print(paste0("Run time: ", fptm, " hours"))

# reproducibility:
# This script was styled and linted.
# Code excluded here, as script is submitted as an array.
