### Functions for CCC in AD project 
# Tabea M. Soelter 

## make_seurat_object
# A function which takes a path to sample folders with the three CellRanger output files and creates a merged seurat object
make_seurat_object <- function(path){
  counts_list <- list.dirs(path, 
                           full.names = TRUE, 
                           recursive = FALSE)
  object_list <- vector('list')
  for (i in counts_list){
    counts <- Read10X(i)
    print(i)
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    seurat_object <- CreateSeuratObject(counts = counts, project = sample_name, min.features = 200)
    if (str_sub(sample_name, - 2, - 1) == "AD") {
      seurat_object$orig.ident <- "AD"
    } else {
      seurat_object$orig.ident <- "CTRL"
    }
    object <- seurat_object
    object_list[[i]] <- object
  } 
  print("Making diseased Seurat Object")
  AD_list <- object_list[grepl("AD", names(object_list))]
  for (i in names(AD_list)) {
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    AD_list[[i]] <- RenameCells(AD_list[[i]],
                                add.cell.id = sample_name)
  }
  diseased <- Merge_Seurat_List(AD_list)
  
  print("Making control Seurat Object")
  CTRL_list <- object_list[grepl("CTRL", names(object_list))]
  for (i in names(CTRL_list)) {
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    CTRL_list[[i]] <- RenameCells(CTRL_list[[i]],
                                  add.cell.id = sample_name)
  }
  control <- Merge_Seurat_List(CTRL_list)
  
  print("Making merged Seurat Object")
  merged_seurat <- merge(x = control,
                         y = diseased,
                         add.cell.id = c("CTRL", "AD"))
  return(merged_seurat)
}

## calculate_qc
# A function which calculates quality control metrics for a merged seurat object
calculate_qc <- function(seurat_object){
  seurat_object$log10GenesPerUMI <- log10(seurat_object$nFeature_RNA) / log10(seurat_object$nCount_RNA)
  seurat_object$mitoRatio <- PercentageFeatureSet(object = seurat_object, 
                                                  pattern = "^MT-")
  seurat_object$mitoRatio <- seurat_object@meta.data$mitoRatio / 100
  return(seurat_object)
}

## format_metadata
# A function which extracts and formats metadata of a seurat object
format_metadata <- function(seurat_object){
  metadata <- seurat_object@meta.data
  metadata$cells <- rownames(metadata)
  # Create sample column -----
  metadata$sample <- NA
  metadata$sample[which(str_detect(metadata$cells, "^CTRL_"))] <- "ctrl"
  metadata$sample[which(str_detect(metadata$cells, "^AD_"))] <- "AD"
  # Rename columns -----
  metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident,
                                         nUMI = nCount_RNA,
                                         nGene = nFeature_RNA)
  return(metadata)
}

## plot_qc
# A function which takes seurat metadata and plots quality control metrics for filtering purposes
plot_qc <- function(metadata) {
  # Visualize the number of cell counts per sample
  number_of_cells <- metadata %>% 
    ggplot(aes(x = sample, fill = sample)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells") 
  # Visualize the number UMIs/transcripts per cell
  number_of_umis <- metadata %>% 
    ggplot(aes(color = sample, x = nUMI, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  # Visualize the distribution of genes detected per cell
  dist_genes_per_cell <- metadata %>% 
    ggplot(aes(color = sample, x = nGene, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  novelty_score <- metadata %>%
    ggplot(aes(x = log10GenesPerUMI, color = sample, fill = sample)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  # Visualize the distribution of mitochondrial gene expression detected per cell
  dist_mito_gex <- metadata %>% 
    ggplot(aes(color = sample, x = mitoRatio, fill = sample)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2)
  # Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells   with low numbers of genes/UMIs
  cor <- metadata %>% 
    ggplot(aes(x = nUMI, y = nGene, color = mitoRatio)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method = lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~sample)
  # Plot QC metrics
  plot(number_of_cells) 
  plot(number_of_umis)
  plot(dist_genes_per_cell)
  plot(novelty_score)
  plot(dist_mito_gex)
  plot(cor)
}

## cell_cycle_effects
# A function which calculates and plots the effect of cell cycle on the data using a filtered seurat object as input. It also performs log normalization, scaling, and dimension reduction using PCA
cell_cylce_effects <- function(filtered_seurat){
  # log normalize -----
  filtered_seurat <- NormalizeData(filtered_seurat)
  # cell cylce markers -----
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  # score cells based in gex of genes -----
  filtered_seurat <- CellCycleScoring(filtered_seurat,
                                      g2m.features = g2m.genes,
                                      s.features = s.genes)
  filtered_seurat <- FindVariableFeatures(filtered_seurat,
                                          selection.method = "vst",
                                          verbose = FALSE)
  # scale data -----
  filtered_seurat <- ScaleData(filtered_seurat)
  # run pca -----
  filtered_seurat <- RunPCA(filtered_seurat, approx = FALSE)
  # plot pca -----
  elbow <- ElbowPlot(filtered_seurat, reduction = "pca", ndims = 50)
  # plot cell cycle scoring -----
  cell_cycle_plot <- DimPlot(filtered_seurat,
                             reduction = "pca",
                             group.by = "Phase",
                             split.by = "Phase")
  plot(cell_cycle_plot)
  plot(elbow)
  return(filtered_seurat)
}

## harmony_integration
# A function which integrates a seurat object using harmony
harmony_integration <- function(seurat_object, dims){
  seurat_object <- RunHarmony(seurat_object,
                              group.by.vars = "sample",
                              reduction = "pca",
                              dims.use = dims, assay.use = "RNA")
  # Here we use all PCs computed from Harmony for UMAP calculation -----
  seurat_object <- RunUMAP(seurat_object, dims = dims, reduction = "harmony", reduction.name = "umap_harmony")
  return(seurat_object)
}

## find_clusters
# A function which finds clusters for a Seurat Object at user determined resolutions
find_clusters <- function(object, dims, reduction, resolutions) {
  # set reduction method to harmony -----
  object <- FindNeighbors(object, 
                          dims = dims, 
                          reduction = reduction)
  # clustering (Leiden aka algorithm 4)
  for (res in resolutions) {
    object <- FindClusters(object,
                           graph.name = "RNA_snn",
                           resolution = res,
                           algorithm = 4,
                           method = "igraph")
  } 
  return(object)
}

## find_markers
# A function allowing for individual FindMarkers analyses for a defined list of clusters
find_markers <- function(object, resolution, identities, value){
  # set resolution ----------
  object <- SetIdent(object, value = resolution)
  # create empty df ----------
  top_markers <- data.frame()
  # iterating through all unidentified clusters ----------
  for (i in identities) {
    # find marker genes ----------
    markers <- FindMarkers(object, ident.1 = i, max.cells.per.ident = 100, logfc.threhold = 0.25, only.pos = TRUE)
    # print out which cluster was completed ----------
    print(paste0("Markers for cluster ", i))
    # filter markers by specified value and adjusted p-value ----------
    markers_filt <- markers %>% top_n(-value, p_val_adj) %>% add_column(cluster = i)
    # bind empty df and filtered markers together ----------
    top_markers <- rbind(top_markers, markers_filt)
  }
  return(top_markers)
}

## prep_liana
# Ensures the "RNA" assay is used before plotting each object's UMAP split by condition (AD, CTRL) and stacked bar plots of cell type proportions across conditions. Plots are saved to their respective "results/final_outputs/" directories.
prep_liana <- function(object_list) {
  for(name in names(object_list)) {
    object <- get(name)
    DefaultAssay(object) <- "RNA"
    # UMAP split by condition ----------
    umap <- DimPlot(object,
                    label = FALSE,
                    group.by = "orig.ident",
                    cols = "Paired",
                    shuffle = TRUE)
    ggsave(filename = "UMAP_condition.png",
           path = paste0(here("results", "final_outputs", name, "/")),
           plot = umap)
    # stacked barplot of cell type proportions ----------
    df <- table(Idents(object), object$orig.ident)
    df <- as.data.frame(df)
    df$Var1 <- as.character(df$Var1)
    barplot <- ggplot(df, aes(x = fct_rev(Var1), y = Freq, fill = Var2)) +
      theme_bw(base_size = 15) +
      geom_col(position = "fill", width = 0.5) +
      xlab("Cell Type") +
      ylab("Proportion") +
      scale_fill_manual(values = brewer.pal(12, "Paired")) +
      theme(legend.title = element_blank()) +
      coord_flip()
    ggsave(filename = "celltype_proportions_stackedbarplot.png",
           path = paste0(here("results", "final_outputs", name, "/")),
           plot = barplot)
  }
}

## split_objects
# Split multiple Seurat objects by their condition (orig.ident in my data) and return a list with split objects.
split_objects <- function(object, active_assay = "RNA") {
  if(object@active.assay == active_assay) {
    split_object <- SplitObject(object, split.by = "orig.ident")
  } else {
    print(paste0("Active assay is not RNA in ", i))
  }
  return(split_object)
}

## liana_prioritization
# A function which identified ligands between sender and receiver cells for a seurat object. 
liana_prioritization <- function(object, sender, receiver) {
  results <- liana_wrap(object) %>%
    liana_aggregate() %>% 
    filter(source %in% sender & target == receiver) %>%
    dplyr::rename(ligand = ligand.complex, receptor = receptor.complex) %>%
    arrange(aggregate_rank) %>%
    mutate(id = fct_inorder(paste0(ligand, " -> ", receptor)))
  return(results)
}

## filter_ligands
# A function which filters duplicate ligands and those present in NicheNet prior
filter_ligands <- function(df){
  ligands <- unique(df$ligand)
  ligands <- ligands[ligands %in% colnames(ligand_target_matrix)]
  print(length(ligands))
  return(ligands)
}

## prep_NicheNet
# A function which prepares the prioritized NicheNet outputs for mapping to STRINGdb PPI and returns a dataframe with a target and a sender column 
prep_NicheNet <- function(prioritizedNicheNet, cond_niche){
  top_ligand_niche_df <- prioritizedNicheNet$prioritization_tbl_ligand_receptor %>% 
    select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
    group_by(ligand) %>% 
    top_n(1, prioritization_score) %>% 
    ungroup() %>% 
    select(ligand, receptor, niche) %>% 
    rename(top_niche = niche)
  ligand_prioritized_tbl_oi <- prioritizedNicheNet$prioritization_tbl_ligand_receptor %>% 
    select(niche, sender, receiver, ligand, prioritization_score) %>% 
    group_by(ligand, niche) %>% 
    top_n(1, prioritization_score) %>% 
    ungroup() %>% 
    distinct() %>% 
    inner_join(top_ligand_niche_df) %>% 
    filter(niche == top_niche) %>% 
    group_by(niche) %>% 
    top_n(50, prioritization_score) %>% 
    ungroup() 
  targets <- ligand_prioritized_tbl_oi %>% 
    inner_join(prioritizedNicheNet$prioritization_tbl_ligand_target,
               by = c("niche",
                      "receiver",
                      "sender",
                      "ligand",
                      "receptor",
                      "prioritization_score"),
               multiple = "all") %>% 
    filter(niche == cond_niche) %>%
    select(sender, target, receiver) %>% 
    group_by (sender, target, receiver) %>% 
    distinct(target) %>% 
    ungroup()
  return(targets)
}

## map_genes
# A function which joins a disease gene list with a prioritized gene df in order to map the genes to the STRINGdb PPI
map_genes <- function(prioritized_targets, gene_list){
  colnames(gene_list) <- "target"
  genes <- prioritized_targets %>% 
    select(target) %>% 
    full_join(gene_list) %>% 
    unique() %>% 
    rename(gene = "target") %>% 
    as.data.frame()
  mapped_genes <- ppi$map(genes, 
                          "gene", 
                          removeUnmappedRows = TRUE)
  return(mapped_genes)
}

## normalize
# A function which performs min-max scaling on a df 
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

## compile_gex
# A function to read in a Seurat object, filter it by receiver cell type and condition, as well as pull GEx values. GEx values are then summed and normalized.
compile_gex <- function(object, receiver, condition, genes_ppi){
  seurat_obj <- readRDS(here("data", 
                             object))
  target_all <- subset(seurat_obj, 
                       idents = receiver)
  target_cond <- subset(target_all, 
                        subset = orig.ident == condition)
  counts_ppi <- FetchData(target_cond, 
                          slot = "counts", 
                          vars = genes_ppi)
  counts_m <- as.data.frame(colMeans(counts_ppi))
  counts_norm <- normalize(counts_m)
  return(counts_norm)
}

## calculate_weights
# A function to calculate edge weights for PPI based on GEx and STRINGdb score.
calculate_weights <- function(x, counts){
  z <- x[1]
  z_count <- counts %>% filter(row.names(counts) %in% z) %>% as.numeric()
  y <- x[2]
  y_count <- counts %>% filter(row.names(counts) %in% y) %>% as.numeric()
  edge_weight <- (z_count + y_count)*(as.numeric(x[4]))
  return(edge_weight)
}

## inverse_weights
# A function to inverse weights from Djikstra's shortest path algorithm and calculation
inverse_weights <- function(df){
  edge_weight_t <- max(df[5]) - df[5] + 1 
  return(edge_weight_t)
} 

## create_igraph_object
# A function which serves as a wrapper function to filter target genes, map them to STRING identifiers, add gex values, and calculate edge weights before returning a list of the different igraph objects returned.
# This function will do everything necessary to create an igraph object for target genes from 1 seurat object across conditions.
create_igraph_object <- function(ppi, prioritized_targets, disease_gene_list, condition, receiver, seurat_object) {
  for (cond in condition) {
    if(cond == "AD") {
      prioritized_targets_AD <- subset(prioritized_targets, 
                                       sub(".*_", "", prioritized_targets$receiver) == "AD", 
                                       drop = TRUE)
      mapped_genes_AD <- map_genes(prioritized_targets_AD,
                                   gene_list = disease_gene_list)
      # get interactions for genes of interest -----
      ppi_tmp <- ppi$get_interactions(mapped_genes_AD$STRING_id)
      # create df of ppi -----
      ppi_df <- data.frame(from = mapped_genes_AD[match(ppi_tmp$from, mapped_genes_AD$STRING_id), 1],
                           to = mapped_genes_AD[match(ppi_tmp$to, mapped_genes_AD$STRING_id), 1],
                           score = ppi_tmp$combined_score) %>% unique()
      # compile genes from ppi -----
      genes_ppi <- ppi_df %>% select(from, to) %>% unlist() %>% unique()
      # pull GEx, filter and normalize it -----
      counts_norm <- compile_gex(object = seurat_object, 
                                 receiver = receiver, 
                                 condition = cond, 
                                 genes_ppi = genes_ppi)
      # Scale StringDB scores -----
      ppi_df$score_scaled <- ppi_df$score/1000
      # Calculate edge weights for interactions -----
      edge_weight <- apply(ppi_df, 1, calculate_weights, counts = counts_norm)
      # Add new column with edge weights after removing NAs, as well as a new column with the transformed weights through inversion (the largest value will become the smallest and vice versa, as Dijkstra's will calculate the shortest path) -----
      ppi_weights <- cbind(ppi_df, edge_weight) %>% 
        filter(!is.na(edge_weight))
      
      ppi_weights <- ppi_weights %>% 
        mutate(edge_weight_t = as.numeric(unlist(inverse_weights(ppi_weights))))
      # Create igraph object -----
      AD_igraph <- graph_from_data_frame(ppi_weights)
      # Set edge attributes as the calculated edge weight -----
      E(AD_igraph)$weight <- ppi_weights$edge_weight_t
    } else {
      prioritized_targets_control <- subset(prioritized_targets, 
                                            sub(".*_", "", prioritized_targets$receiver) == "Control", 
                                            drop = TRUE)
      mapped_genes_control <- map_genes(prioritized_targets_control,
                                        gene_list = disease_gene_list)
      ppi_tmp <- ppi$get_interactions(mapped_genes_control$STRING_id)
      ppi_df <- data.frame(from = mapped_genes_control[match(ppi_tmp$from, mapped_genes_control$STRING_id), 1],
                           to = mapped_genes_control[match(ppi_tmp$to, mapped_genes_control$STRING_id), 1],
                           score = ppi_tmp$combined_score) %>% unique()
      genes_ppi <- ppi_df %>% select(from, to) %>% unlist() %>% unique()
      counts_norm <- compile_gex(object = seurat_object, 
                                 receiver = receiver, 
                                 condition = cond, 
                                 genes_ppi = genes_ppi)
      ppi_df$score_scaled <- ppi_df$score/1000
      edge_weight <- apply(ppi_df, 1, calculate_weights, counts = counts_norm)
      ppi_weights <- cbind(ppi_df, edge_weight) %>%
        filter(!is.na(edge_weight))
      ppi_weights <- ppi_weights %>%
        mutate(edge_weight_t = as.numeric(unlist(inverse_weights(ppi_weights))))
      CTRL_igraph <- graph_from_data_frame(ppi_weights)
      E(CTRL_igraph)$weight <- ppi_weights$edge_weight_t
    }
  } 
  return(list(AD_igraph, CTRL_igraph))
}

## CellTypeHeatmaps
# A function which subsets prioritized targets by sender cell and plots a heatmap for it
CellTypeHeatmaps <- function(distMatrix, prioritized_targets, cell_type){
  cell <- prioritized_targets %>%
    filter(sender == cell_type) %>% 
    select(target)
  distMatrix_filt <- distMatrix %>% 
    as.data.frame() %>% 
    filter(rownames(distMatrix) %in% cell$target) %>% 
    as.matrix()
  plot(Heatmap(distMatrix_filt, 
               col = viridis(100),
               cluster_columns = TRUE, 
               cluster_rows = TRUE,
               show_column_dend = FALSE,
               show_row_dend = FALSE,
               row_names_gp = gpar(fontsize = 10, 
                                   fontcolor = "red"),
               column_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(title = "distance")))
}