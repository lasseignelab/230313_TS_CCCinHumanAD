#### Shortest path calculations using Dijkstra's algorithm.
### Author: Tabea M. Soelter
### Date: 2023-05-16

## Goal: Calculating the shortest path between all nodes in the PPI using Dijkstra's. 

## set seed
set.seed(42)
print("seed set")

## load packages
library(tidyverse)
library(igraph)
print("packages loaded")

## enable usage of args
args <- R.utils::commandArgs(trailingOnly = TRUE)
print("enabled args usage")

## get dataset name
name <- sub("_igraph.rds.*", "", basename(args[1]))
print(name)

## load in data
igraph <- readRDS(args[1])
print("loaded igraph object")

## calculate shortest path
# I am using Dijkstra's shortest path algorithm to calculate the shortest paths between all possible nodes. 
# As Dijkstra's calculates the shortest path, we previously inverted the weights we calculated, so the edges with the highest GEx will become the lowest values and therefore be favored during the shortest path calculation.
distance_matrix <- distances(igraph,
                             algorithm = "dijkstra",
                             v = V(igraph),
                             to = V(igraph))
print("calculated distances")

## save distance matrix
print("saving distance matrix")
saveRDS(distance_matrix, file = paste0("/data/user/tsoelter/projects/230313_TS_CCCinHumanAD/data/shortest_path/", name, "_distance_matrix.rds"))