### Functions for CCC in AD project 
# Tabea M. Soelter 

## remove_ambientRNA
# A function which removes ambient RNA from h5 files outputted from Cell Ranger for single cell data.
remove_ambientRNA <- function(inputs, outputs, plots) {
  print("Making list of objects")
  counts_list <- list.dirs(inputs, 
                           full.names = TRUE, 
                           recursive = FALSE)
  pdf(paste0(plots, "rho_density_plots.pdf"))
  for (i in counts_list) {
    set.seed(42)
    # load in cell ranger h5 outputs
    print("Loading cell ranger h5 objects")
    filt_matrix <- Read10X_h5(paste0(i, "/filtered_feature_bc_matrix.h5"))
    raw_matrix <- Read10X_h5(paste0(i, "/raw_feature_bc_matrix.h5"))
    # create seurat object
    print("Making seurat object")
    object <- CreateSeuratObject(counts = filt_matrix)
    # make soup channel object
    print("Making soup channel object")
    sco <- SoupChannel(raw_matrix, filt_matrix)
    # get cluster info
    print("Get cluster info")
    object <- SCTransform(object, verbose = FALSE)
    object <- RunPCA(object, approx = FALSE, verbose = FALSE)
    object <- RunUMAP(object, dims = 1:30, verbose = FALSE)
    object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
    object <- FindClusters(object, verbose = FALSE)
    # ading metadata to soup channel object
    meta <- object@meta.data
    umap <- object@reductions$umap@cell.embeddings
    sco <- setClusters(sco, setNames(meta$seurat_clusters, rownames(meta)))
    sco <- setDR(sco, umap)
    # Analyzing the soup
    print("Profiling the soup")
    sco <- autoEstCont(sco)
    # Create integer matrix
    adjusted_matrix <- adjustCounts(sco, roundToInt = TRUE)
    # save
    print("Saving filtered objects")
    sample_name <- basename(i)
    DropletUtils::write10xCounts(paste0(outputs, sample_name), adjusted_matrix) 
  }
  dev.off()
}

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
  AD_list <- object_list[grepl("_AD", names(object_list))]
  for (i in names(AD_list)) {
    sample_name <- basename(i)
    sample_name <- gsub("[[:punct:]]", "", sample_name)
    AD_list[[i]] <- RenameCells(AD_list[[i]],
                                add.cell.id = sample_name)
  }
  diseased <- Merge_Seurat_List(AD_list)
  
  print("Making control Seurat Object")
  CTRL_list <- object_list[grepl("_CTRL", names(object_list))]
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
  # Rename columns -----
  metadata <- metadata %>% dplyr::rename(seq_folder = orig.ident,
                                         nUMI = nCount_RNA,
                                         nGene = nFeature_RNA)
  return(metadata)
}

## plot_qc
# A function which takes seurat metadata and plots quality control metrics for filtering purposes
plot_qc <- function(metadata) {
  # Visualize the number of cell counts per condition
  number_of_cells <- metadata %>% 
    ggplot(aes(x = seq_folder, fill = seq_folder)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    ggtitle("NCells") 
  # Visualize the number UMIs/transcripts per cell
  number_of_umis <- metadata %>% 
    ggplot(aes(color = seq_folder, x = nUMI, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = 500)
  # Visualize the distribution of genes detected per cell
  dist_genes_per_cell <- metadata %>% 
    ggplot(aes(color = seq_folder, x = nGene, fill = seq_folder)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    geom_vline(xintercept = 300)
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
  novelty_score <- metadata %>%
    ggplot(aes(x = log10GenesPerUMI, color = seq_folder, fill = seq_folder)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8)
  # Visualize the distribution of mitochondrial gene expression detected per cell
  dist_mito_gex <- metadata %>% 
    ggplot(aes(color = seq_folder, x = mitoRatio, fill = seq_folder)) + 
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
    facet_wrap(~seq_folder)
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
    markers <- FindMarkers(object, ident.1 = i, max.cells.per.ident = 100, logfc.threshold = 0.25, only.pos = TRUE)
    # print out which cluster was completed ----------
    print(paste0("Markers for cluster ", i))
    # filter markers by specified value and adjusted p-value ----------
    markers_filt <- markers %>% top_n(-value, p_val_adj) %>% add_column(cluster = i)
    # bind empty df and filtered markers together ----------
    top_markers <- rbind(top_markers, markers_filt)
  }
  return(top_markers)
}

## make_names_valid
# Making names of interested columns valid (no spacing for example)
# Code adapted from MultiNicheNet
make_names_valid <- function(object) {
  SummarizedExperiment::colData(object)$ident <-
    SummarizedExperiment::colData(object)$ident %>%
    make.names()
  SummarizedExperiment::colData(object)$orig.ident <-
    SummarizedExperiment::colData(object)$orig.ident %>%
    make.names()
  SummarizedExperiment::colData(object)$sample <-
    SummarizedExperiment::colData(object)$sample %>%
    make.names()
  return(object)
}

## multinichenet_wrapper
# A wrapper for original MultiNicheNet code
# Requires a list of single cell experiment objects as input and returns a list of multinichenet outputs, which can be used for visualization. 
multinichenet_wrapper <- function(object_list, results_path, celltype_id, sample_id, group_id, lr_network, batches, contrasts_oi, contrast_tbl, covariates, empirical_pval, p_val_adj, cores_system, ligand_target_matrix, prioritizing_weights) {
  multinichenet_objects <- tibble::lst()
  for (name in names(object_list)) {
    object <- get(name)
    print(name)
    # get sender and receiver cell types from object ---------------
    senders_oi <- colData(object)[,celltype_id] %>% unique()
    receivers_oi <- colData(object)[,celltype_id] %>% unique()
    object <- object[, colData(object)[,celltype_id] %in% c(senders_oi, receivers_oi)]
    print("Grabbed senders and receivers")
    # determine whether all cell types have enough cells -----------
    min_cells <- 10
    abundance_expression_info <- get_abundance_expression_info(sce = object,
                                                               sample_id = sample_id,
                                                               group_id = group_id,
                                                               celltype_id = celltype_id,
                                                               min_cells = min_cells,
                                                               senders_oi = senders_oi,
                                                               receivers_oi = receivers_oi,
                                                               lr_network = lr_network,
                                                               batches = batches)
    plot1 <- abundance_expression_info$abund_plot_sample
    plot2 <- abundance_expression_info$abund_plot_group
    pdf(here(paste0(results_path, name, "/ccc/abundance_expression.pdf")))
    print(plot1)
    print(plot2)
    dev.off()
    print("Plotted abundance expression")
    # Get differential expression information
    DE_info <- get_DE_info(sce = object,
                           sample_id = sample_id,
                           group_id = group_id,
                           celltype_id = celltype_id,
                           batches = batches,
                           covariates = covariates,
                           contrasts_oi = contrasts_oi,
                           min_cells = min_cells)
    print("Calculated DE")
    # Calculate p-values
    if(empirical_pval == FALSE) {
      celltype_de <- DE_info$celltype_de$de_output_tidy
    } else {
      celltype_de <- DE_info_emp$de_output_tidy_emp %>% dplyr::select(-p_val, -p_adj) %>%
        dplyr::rename(p_val = p_emp, p_adj = p_adj_emp)
    }
    sender_receiver_de = combine_sender_receiver_de(
      sender_de = celltype_de,
      receiver_de = celltype_de,
      senders_oi = senders_oi,
      receivers_oi = receivers_oi,
      lr_network = lr_network)
    print("Calculated p-values")
    # Identify ligand activities
    fraction_cutoff <- 0.05
    n.cores <- min(cores_system, sender_receiver_de$receiver %>% unique() %>% length())
    if(length(receivers_oi) > n.cores) {
      print(paste0("Core req not met. Minimum number of cores needed: ",
                   length(receivers_oi)))
      stop()
    } else {
      print("Calculating ligand activities")
      ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
        get_ligand_activities_targets_DEgenes(
          receiver_de = celltype_de,
          receivers_oi = receivers_oi,
          ligand_target_matrix = ligand_target_matrix,
          logFC_threshold = 0.50,
          p_val_threshold = 0.05,
          p_val_adj = p_val_adj,
          top_n_target = 250,
          verbose = FALSE, 
          n.cores = n.cores
        )))
    }
    # Prepare for prioritization
    print("Prioritizing interactions")
    sender_receiver_tbl <- sender_receiver_de %>%
      dplyr::distinct(sender, receiver)
    metadata_combined <- colData(object) %>% tibble::as_tibble()
    if(!is.na(batches)){
      grouping_tbl <- metadata_combined[,c(sample_id, group_id, batches)] %>%
        tibble::as_tibble() %>%
        dplyr::distinct()
      colnames(grouping_tbl) <- c("sample", "group", batches)
    } else {
      grouping_tbl <- metadata_combined[,c(sample_id, group_id)] %>%
        tibble::as_tibble() %>%
        dplyr::distinct()
      colnames(grouping_tbl) <- c("sample", "group")
    }
    # Prioritize interactions
    prioritization_tables <- suppressMessages(generate_prioritization_tables(
      sender_receiver_info = abundance_expression_info$sender_receiver_info,
      sender_receiver_de = sender_receiver_de,
      ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
      contrast_tbl = contrast_tbl,
      sender_receiver_tbl = sender_receiver_tbl,
      grouping_tbl = grouping_tbl,
      prioritizing_weights = prioritizing_weights,
      fraction_cutoff = fraction_cutoff, 
      abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
      abundance_data_sender = abundance_expression_info$abundance_data_sender
    ))
    # Prioritize by correlation coefficients
    print("Prior knowledge correlations")
    lr_target_prior_cor <- lr_target_prior_cor_inference(
      prioritization_tables$group_prioritization_tbl$receiver %>%
        unique(),
      abundance_expression_info,
      celltype_de,
      grouping_tbl,
      prioritization_tables,
      ligand_target_matrix,
      logFC_threshold = 0.50,
      p_val_threshold = 0.05,
      p_val_adj = p_val_adj)
    # Create combined output
    print("Creating multinichenet output")
    multinichenet_output = list(
      celltype_info = abundance_expression_info$celltype_info,
      celltype_de = celltype_de,
      sender_receiver_info = abundance_expression_info$sender_receiver_info,
      sender_receiver_de =  sender_receiver_de,
      ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
      prioritization_tables = prioritization_tables,
      grouping_tbl = grouping_tbl,
      lr_target_prior_cor = lr_target_prior_cor
    )
    new_name <- paste0(name, "_multinichenet_output")
    multinichenet_objects[new_name] <- list(multinichenet_output)
  }
  return(multinichenet_objects)
}

## fea
# A function to perform pathway nnalysis using gprofiler2 and filters for the top 50 pathways for plotting purposes.
# Adapted from Lizzy Wilk
fea <- function(genes){
  set.seed(42)
  # create gprofiler2 query ----------
  fea_result <- gost(query = genes,
                     organism = "hsapiens",
                     ordered_query = FALSE,
                     multi_query = FALSE,
                     significant = TRUE,
                     exclude_iea = FALSE,
                     measure_underrepresentation = FALSE,
                     evcodes = TRUE,
                     user_threshold = 0.05,
                     correction_method = "bonferroni",
                     domain_scope = "annotated",
                     numeric_ns = "",
                     sources = NULL,
                     as_short_link = FALSE) 
  # remove arbitrary pathways ----------
  fea_result <- fea_result$result %>% filter(term_size < 1000 | term_size > 10)
  # keep only top 50 pathways for plotting purposes ---------
  fea_result_filt <- fea_result %>% top_n(n = 15)
  return(fea_result_filt)
}

## bubbleplot
# A function to create a DotPlot for gprofiler2 results
# Adapted from Lizzy Wilk
bubbleplot <- function(fea_result, region, file_path, file_name){
  plot <- ggplot(fea_result,
                 aes(x = direction,
                     y = reorder(term_name, -p_value),
                     size = intersection_size,
                     fill = p_value)) +
    geom_point(alpha = 0.7, shape = 21) +
    scale_size(range = c(2, 10), name = "Intersection Size") + 
    scale_fill_distiller(palette = "Purples") + 
    labs(x = "Direction", y = "Functional Enrichment Terms") +
    theme_minimal() + 
    ggtitle(paste0("Top Enriched Terms for Predicted Target\nGenes in ", region)) +
    theme(axis.text = element_text(face = "bold"))
  ggsave(filename = file_name,
         path = paste0(here(file_path)),
         plot = plot,
         width = 10,
         height = 9,
         bg = "white")
  return(plot)
}

## combined_fea
# A function to get pathways for up- and down-regulated information of genes submitted 
combined_fea <- function(genes, receiver){
  genes_AD <- genes %>% filter(receiver == receiver)
  upgenes <- as.list(genes_AD %>% filter(direction == "up"))
  up_fea_filt <- fea(genes = upgenes) %>% mutate(direction = "upregulated")
  downgenes <- as.list(genes_AD %>% filter(direction == "down"))
  down_fea_filt <- fea(genes = downgenes) %>% mutate(direction = "downregulated")
  comb_fea <- rbind(up_fea_filt, down_fea_filt)
  fea_result <- list(downregulated = down_fea_filt, upregulated = up_fea_filt, combined = comb_fea)
  return(fea_result)
}

## pathway_analysis
# A wrapper function to prepare inputs for pathway analysis, perform pathway analysis using gprofiler2 on up/down regulated target genes from nichenet
pathway_analysis <- function(nichenet_output, receiver_type) {
  # subset output for df with target expression info -----
  target_expression <- nichenet_output$exprs_tbl_target
  # filter expression df for necessary columns and cell type as well as add up/down gex info -----
  neuro_targets_direction <- target_expression %>%
    select(receiver, target, target_expression_scaled) %>%
    filter(receiver == receiver_type) %>%
    mutate(direction = ifelse(target_expression_scaled > 0, "up", "down")) %>%
    unique()
  # pathway analysis of up/down genes
  fea_res <- combined_fea(neuro_targets_direction, receiver = receiver)
  return(fea_res)
}

## map_genes
# A function which joins a disease gene list with a prioritized gene df in order to map the genes to the STRINGdb PPI
map_genes <- function(targets, gene_list){
  colnames(gene_list) <- "target"
  genes <- targets %>% 
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
#compile_gex <- function(object, receiver, condition, genes_ppi){
#  target_all <- subset(object, 
#                       idents = receiver)
#  target_cond <- subset(target_all, 
#                        subset = orig.ident == condition)
#  counts_ppi <- FetchData(target_cond, 
#                          slot = "data", 
#                          vars = genes_ppi)
#  counts_m <- as.data.frame(colMeans(counts_ppi))
#  counts_norm <- normalize(counts_m)
#  return(counts_norm)
#}

compile_gex <- function(object, receiver, condition, genes_ppi){
  target_all <- subset(object, 
                       idents = receiver)
  target_cond <- subset(target_all, 
                        subset = orig.ident == condition)
  counts_ppi <- FetchData(target_cond, 
                          slot = "data", 
                          vars = genes_ppi)
  counts_norm <- as.data.frame(colMeans(counts_ppi))
  return(counts_norm)
}

## calculate_weights
# A function to calculate edge weights for PPI based on GEx and STRINGdb score.
#calculate_weights <- function(x, counts){
#  z <- x[1]
#  z_count <- counts %>% filter(row.names(counts) %in% z) %>% as.numeric()
#  y <- x[2]
#  y_count <- counts %>% filter(row.names(counts) %in% y) %>% as.numeric()
#  edge_weight <- (z_count + y_count)*(as.numeric(x[4]))
#  return(edge_weight)
#}

calculate_weights <- function(x, counts){
  z <- x[1]
  z_count <- counts %>% filter(row.names(counts) %in% z) %>% as.numeric()
  y <- x[2]
  y_count <- counts %>% filter(row.names(counts) %in% y) %>% as.numeric()
  max_weight <- max(z_count * y_count)
  edge_weight <- ((z_count * y_count)/max_weight)*(as.numeric(x[4]))
  return(edge_weight)
}

## inverse_weights
# A function to inverse weights from Djikstra's shortest path algorithm and calculation
#inverse_weights <- function(df){
#  edge_weight_t <- max(df[5]) - df[5] + 1 
#  return(edge_weight_t)
#} 

inverse_weights <- function(df){
  edge_weight_t <- 1 - df[5] 
  return(edge_weight_t)
} 

## create_igraph_object
# A function which serves as a wrapper function to filter target genes, map them to STRING identifiers, add gex values, and calculate edge weights before returning a list of the different igraph objects returned.
# This function will do everything necessary to create an igraph object for target genes from 1 seurat object across conditions.
create_igraph_object <- function(condition, receiver, seurat_object, mapped_genes, ppi_tmp) {
  # create df of ppi -----
  ppi_df <- data.frame(from = mapped_genes[match(ppi_tmp$from, mapped_genes$STRING_id), 1],
                       to = mapped_genes[match(ppi_tmp$to, mapped_genes$STRING_id), 1],
                       score = ppi_tmp$combined_score) %>% unique()
  # compile genes from ppi -----
  genes_ppi <- ppi_df %>% select(from, to) %>% unlist() %>% unique()
  for (cond in condition) {
    if(cond == "AD") {
      # pull GEx, filter and normalize it -----
      counts_norm <- compile_gex(object = seurat_object, 
                                 receiver = receiver, 
                                 condition = cond, 
                                 genes_ppi = genes_ppi)
      # Scale StringDB scores -----
      ppi_df$score_scaled <- ppi_df$score/1000
      # Calculate edge weights for interactions -----
      print("calculating AD edge weights")
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
      counts_norm <- compile_gex(object = seurat_object, 
                                 receiver = receiver, 
                                 condition = cond, 
                                 genes_ppi = genes_ppi)
      ppi_df$score_scaled <- ppi_df$score/1000
      print("calculating CTRL edge weights")
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

## filter_nichenet
# A function that will filter a multinichenet output by pearson and spearman correlation values above/below 0.33/-0.33
filter_nichenet <- function(object) {
  # Filter correlated object by pearson and spearman correlations
  lr_target_prior_cor_filtered <- object$lr_target_prior_cor %>%
    inner_join(object$ligand_activities_targets_DEgenes$ligand_activities %>%
                 distinct(ligand, target, direction_regulation, contrast)) %>%
    inner_join(contrast_tbl) %>%
    filter(group == group_oi, receiver %in% receiver_oi, sender %in% sender_oi)
  
  lr_target_prior_cor_filtered_up <- lr_target_prior_cor_filtered %>%
    filter(direction_regulation == "up") %>%
    filter((rank_of_target < top_n_target) & (pearson > 0.33 | spearman > 0.33))
  
  lr_target_prior_cor_filtered_down <- lr_target_prior_cor_filtered %>%
    filter(direction_regulation == "down") %>%
    filter((rank_of_target < top_n_target) & (pearson < -0.33 | spearman < -0.33))
  
  lr_target_prior_cor_filtered <- bind_rows(lr_target_prior_cor_filtered_up,
                                            lr_target_prior_cor_filtered_down)
  # Create new column with ligand-receptor-target only info
  lr_target_prior_cor_filtered$ligand_receptor_target <- gsub(
    "(_[^_]+_)(.*?)(_[^_]+_)",
    "\\1\\4",
    lr_target_prior_cor_filtered$id_target,
    perl = TRUE)
  # Create new column with ligand-receptor only info
  lr_target_prior_cor_filtered$ligand_receptor <- gsub(
    "^([^_]*_[^_]*).*",
    "\\1",
    lr_target_prior_cor_filtered$id,
    perl = TRUE)
  return(lr_target_prior_cor_filtered)
}

## jaccard
# A function calculating jaccard similarity index values
# Not written by me, this formula can be found in multiple places and forums
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

## calculate_jaccard
# a wrapper for code written by Vishal H. Oza. I modified his code for my own usage.
calculate_jaccard <- function(df, senders, receivers, type) {
  if(type == "ligands") {
    # empty list
    jaccard_results <- list()
    # calculate JI for ligands across all senders by receivers
    for (i in receivers) {
      filtered_receiver <- df %>% filter(receiver == i) %>% select(c("sender", "ligand"))
      
      for (cell_type in senders) {
        genes <- unique(filtered_receiver$ligand[filtered_receiver$sender == cell_type]) 
        
        for (other_cell_type in senders[senders != cell_type]) {
          other_genes <- unique(filtered_receiver$ligand[filtered_receiver$sender == other_cell_type])
          
          jaccard_index <- jaccard(genes, other_genes) 
          
          result_key <- paste(i, cell_type, other_cell_type, sep = "_")
          jaccard_results[[result_key]] <- jaccard_index
        }
      }
    }
  } else if(type == "receptors") {
    # empty list
    jaccard_results <- list()
    # calculate JI for receptors across all senders by receivers
    for (i in receivers) {
      filtered_receiver <- df %>% filter(receiver == i) %>% select(c("sender", "receptor"))
      
      for (cell_type in senders) {
        genes <- unique(filtered_receiver$receptor[filtered_receiver$sender == cell_type]) 
        
        for (other_cell_type in senders[senders != cell_type]) {
          other_genes <- unique(filtered_receiver$receptor[filtered_receiver$sender == other_cell_type])
          
          jaccard_index <- jaccard(genes, other_genes) 
          
          result_key <- paste(i, cell_type, other_cell_type, sep = "_")
          jaccard_results[[result_key]] <- jaccard_index
        }
      }
    }
  } else if(type == "targets") {
    # empty list
    jaccard_results <- list()
    # # calculate JI for targets across all senders by receivers
    for (i in receivers) {
      filtered_receiver <- df %>% filter(receiver == i) %>% select(c("sender", "target"))
      
      for (cell_type in senders) {
        genes <- unique(filtered_receiver$target[filtered_receiver$sender == cell_type]) 
        
        for (other_cell_type in senders[senders != cell_type]) {
          other_genes <- unique(filtered_receiver$target[filtered_receiver$sender == other_cell_type])
          
          jaccard_index <- jaccard(genes, other_genes) 
          
          result_key <- paste(i, cell_type, other_cell_type, sep = "_")
          jaccard_results[[result_key]] <- jaccard_index
        }
      }
    }
  } else {
    print("Please set either ligands, receptors, or targets equal to TRUE")
  } 
  return(jaccard_results)
}

## calculate_jaccard2
# Modified wrapper of calculate_jaccard wrapper. Original code adapted from Vishal H. Oza
calculate_jaccard2 <- function(df, senders, receivers, type) {
  jaccard_results <- list()
  
  for (i in senders) {
    filtered_receiver <- df %>% filter(sender == i) %>% select(c("receiver", type))
    
    if(type == "receptor") {
      for (cell_type in receivers) {
        genes <- unique(filtered_receiver$receptor[filtered_receiver$receiver == cell_type]) 
        
        for (other_cell_type in receivers[receivers != cell_type]) {
          other_genes <- unique(filtered_receiver$receptor[filtered_receiver$receiver == other_cell_type])
        }
      }
    } else if(type == "target") {
      for (cell_type in receivers) {
        genes <- unique(filtered_receiver$target[filtered_receiver$receiver == cell_type]) 
        
        for (other_cell_type in receivers[receivers != cell_type]) {
          other_genes <- unique(filtered_receiver$target[filtered_receiver$receiver == other_cell_type])
        }
      }
    }
    jaccard_index <- jaccard(genes, other_genes) 
    result_key <- paste(i, cell_type, other_cell_type, sep = "_")
    jaccard_results[[result_key]] <- jaccard_index
  }
  
  
  return(jaccard_results)
}

## make_sce
# A function which requires a list of Seurat objects in order to create single cell experiment objects. It also appends the necessary metadata for pseudo-bulking by cell type and sample.
make_sce <- function(object_list) {
  sce_objects <- tibble::lst()
  for(name in names(object_list)) {
    object <- object_list[[name]]
    # raw data
    counts <- object@assays$RNA@counts
    # metadata
    metadata <- object@meta.data
    # add sample_id column as class 'factor'
    metadata$sample_id <- metadata$sample %>% as.factor()
    # add condition information
    metadata$group_id <- metadata$orig.ident
    metadata$group_id <- relevel(metadata$group_id, "CTRL")
    # add cell type information
    metadata$cluster_id <- factor(object@active.ident)
    # make sce object
    sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
    
    new_name <- gsub("processed_seurat", "sce", name) 
    sce_objects[new_name] <- list(sce)}
  
  return(sce_objects)
}

## pseudobulk
# A wrapper function which generates pseudo-bulked data from single cell experiment objects by cell type. The code was adapted form the HBC pseudobulk tutorial.
pseudobulk <- function(object_list) {
  counts_ls_list <- tibble::lst()
  for(name in names(object_list)) {
    print(name)
    sce <- object_list[[name]]
    # prepare for pseudo-bulking
    cluster_names <- levels(colData(sce)$cluster_id)
    print(length(cluster_names))
    sample_names <- levels(colData(sce)$sample_id)
    print(length(sample_names))
    groups <- colData(sce)[, c("cluster_id", "sample_id")]
    aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                    groupings = groups, fun = "sum")
    aggr_counts <- t(aggr_counts)
    # Loop over cell types and extract counts (pseudo-bulk)
    counts_ls <- list()
    for (i in 1:length(cluster_names)) {
      column_idx <- which(tstrsplit(colnames(aggr_counts), "_")[[1]] == cluster_names[i])
      counts_ls[[i]] <- aggr_counts[, column_idx]
      names(counts_ls)[i] <- cluster_names[i]
    }
    # save counts_ls with dataset name
    new_name <- gsub("sce", "counts_ls", name) 
    counts_ls_list[new_name] <- list(counts_ls)}
  
  return(counts_ls_list)
}

## cts_metadata
# A wrapper function to create cell-type-specific metadata for the previously pseudo-bulked count data. This code was adapted from the HBC pseudobulk tutorials.
cts_metadata <- function(object_list, counts_list) {
  metadata_ls_list <- tibble::lst()
  for(name in names(object_list)) {
    print(name)
    sce <- object_list[[name]]
    metadata <- colData(sce) %>% 
      as.data.frame() %>% 
      dplyr::select(group_id, sample_id)
    metadata <- metadata[!duplicated(metadata), ]
    print(dim(metadata))
    rownames(metadata) <- metadata$sample_id
    t <- table(colData(sce)$sample_id,
               colData(sce)$cluster_id)
    match_name <- gsub("sce", "counts_ls", name)
    counts_ls <- counts_list[[match_name]]
    
    metadata_ls <- list()
    
    for (i in 1:length(counts_ls)) {
      df <- data.frame(cluster_sample_id = colnames(counts_ls[[i]]))
      df$cluster_id <- tstrsplit(df$cluster_sample_id, "_")[[1]]
      df$sample_id  <- tstrsplit(df$cluster_sample_id, "_")[[2]]
      idx <- which(colnames(t) == unique(df$cluster_id))
      cell_counts <- t[, idx]
      cell_counts <- cell_counts[cell_counts > 0]
      sample_order <- match(df$sample_id, names(cell_counts))
      cell_counts <- cell_counts[sample_order]
      df$cell_count <- cell_counts
      df <- plyr::join(df, metadata, 
                       by = intersect(names(df), names(metadata)))
      rownames(df) <- df$cluster_sample_id
      metadata_ls[[i]] <- df
      names(metadata_ls)[i] <- unique(df$cluster_id)
    }
    # save counts_ls with dataset name
    new_name <- gsub("sce", "metadata_ls", name) 
    metadata_ls_list[new_name] <- list(metadata_ls)
  }
  return(metadata_ls_list)
}

## deseq2_dea
# A wrapper function which does DEA using DESeq2 for pseudo-bulked single cell data. It automatically saves both significant and all DEGs in a 'pseudobulk' directory at the specified path. This code is originally form the HBC training guide, but was heavily adapted.
deseq2_dea <- function(cell_types, counts_ls, metadata_ls, group_oi, B, padj_cutoff = 0.05, path) {
  cell_type <- cell_types[1]
  print(cell_type)
  ifelse(!dir.exists(here(paste0(path, "pseudobulk/"))),
         dir.create(here(paste0(path, "pseudobulk/"))),
         print("Info: pseudobulk directory already exists"))
  idx <- which(names(counts_ls) == cell_type)
  cluster_counts <- counts_ls[[idx]]
  cluster_metadata <- metadata_ls[[idx]]
  dds <- DESeqDataSetFromMatrix(cluster_counts, 
                                colData = cluster_metadata, 
                                design = ~ group_id)
  dds <- DESeq(dds)
  contrast <- paste(c("group_id", group_oi, "vs", B), collapse = "_")
  res <- results(dds, name = contrast, alpha = 0.05)
  res <- lfcShrink(dds, coef = contrast, res = res, type = "normal")
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()
  # save all results
  write.csv(res_tbl,
            here(paste0(path, "pseudobulk/", cell_type, "_", contrast, "_all_genes.csv")),
            quote = FALSE, 
            row.names = FALSE)
  
  # save sig results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
    dplyr::arrange(padj)
  
  write.csv(sig_res,
            here(paste0(path, "pseudobulk/", cell_type, "_", contrast, "_signif_genes.csv")),
            quote = FALSE, 
            row.names = FALSE)
}

## filter_multinichenet
# A function which filters multinichenet outputs
filter_multinichenet <- function(object_list, sender_oi, receiver_oi, contrast_tbl, top_n_target = 250){
  filtered_objects <- tibble::lst()
  for(name in names(object_list)) {
    multinichenet_output <- object_list[[name]]
    
    
    lr_target_prior_cor_filtered <- multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>%
      unique() %>%
      lapply(function(group_oi) {
        print(group_oi)
        lr_target_prior_cor_filtered <- multinichenet_output$lr_target_prior_cor %>%
          inner_join(multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
                       distinct(ligand, target, direction_regulation, contrast)) %>%
          inner_join(contrast_tbl) %>%
          filter(group == group_oi)
        lr_target_prior_cor_filtered_up <- lr_target_prior_cor_filtered %>%
          filter(direction_regulation == "up") %>%
          filter((rank_of_target < top_n_target) & (pearson > 0.33 | spearman > 0.33))
        lr_target_prior_cor_filtered_down <- lr_target_prior_cor_filtered %>%
          filter(direction_regulation == "down") %>%
          filter((rank_of_target < top_n_target) & (pearson < -0.33 | spearman < -0.33))
        lr_target_prior_cor_filtered <- bind_rows(lr_target_prior_cor_filtered_up,
                                                  lr_target_prior_cor_filtered_down)}) %>%
      bind_rows()
    
    lr_target_prior_cor_filtered <- lr_target_prior_cor_filtered %>%
      filter(sender %in% sender_oi,
             receiver %in% receiver_oi)
    
    new_name <- gsub("multinichenet_output", "lr_target_prior_cor_filtered", name)
    filtered_objects[new_name] <- list(lr_target_prior_cor_filtered)
  }
  return(filtered_objects)
}

## signaling_igraph
# A function that plots and saves signaling GRNs for each LRT pair in multiple datasets, while also returning igraph objects for downstream analysis.
signaling_igraph <- function(lrt_filtered_list, overlap, ccc_combined, ligand_tf_matrix, weighted_networks,
                             lr_network, sig_network, gr_network, plots) {
  igraph_objects_list <- tibble::lst()
  for(name in names(lrt_filtered_list)) {
    lrt_filtered <- lrt_filtered_list[[name]]
    dataset <- gsub("_lr_target_prior_cor_filtered", "", name)
    
    igraph_objects <- tibble::lst()
    for (i in overlap) {
      ccc_combined_sub <- ccc_combined[ccc_combined$id_target == i,]
      ligand_oi <- ccc_combined_sub$ligand
      receptor_oi <- ccc_combined_sub$receptor
      sender_oi <- ccc_combined_sub$sender
      receiver_oi <- ccc_combined_sub$receiver
      id <- ccc_combined_sub$id_target
      print(id)
      
      targets_all <- lrt_filtered %>% 
        filter(ligand == ligand_oi &
                 receiver == receiver_oi &
                 sender == sender_oi &
                 receptor == receptor_oi) %>%
        pull(target) %>%
        unique()
      
      active_signaling_network <-
        nichenetr::get_ligand_signaling_path_with_receptor(
          ligand_tf_matrix = ligand_tf_matrix,
          ligands_all = ligand_oi,
          receptors_all = receptor_oi,
          targets_all = targets_all,
          weighted_networks = weighted_networks,
          top_n_regulators = 2)
      
      
      data_source_network <- nichenetr::infer_supporting_datasources(
        signaling_graph_list = active_signaling_network,
        lr_network = lr_network,
        sig_network = sig_network,
        gr_network = gr_network)
      
      
      active_signaling_network_min_max <- active_signaling_network
      active_signaling_network_min_max$sig <- active_signaling_network_min_max$sig %>%
        mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
      active_signaling_network_min_max$gr <- active_signaling_network_min_max$gr %>%
        mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
      
      
      colors <- c("ligand" = "purple",
                  "receptor" = "orange",
                  "target" = "royalblue",
                  "mediator" = "grey60")
      
      ggraph_signaling_path <- suppressWarnings(
        make_ggraph_signaling_path(active_signaling_network_min_max,
                                   colors,
                                   ligand_oi,
                                   receptor_oi,
                                   targets_all))
      
      png(here(paste0(plots, dataset, "/multinichenet_grn/", id, ".png")))
      plot <- ggraph_signaling_path$plot
      print(plot)
      dev.off()
      
      active_signaling_network_min_max$sig$type <- "ppi"
      active_signaling_network_min_max$gr$type <- "grn"
      
      df <- rbind(active_signaling_network_min_max$sig,
                  active_signaling_network_min_max$gr)
      
      igraph <- graph_from_data_frame(df, directed = TRUE)
      E(igraph)$weight <- df$weight
      
      igraph_objects[[id]] <- igraph
    }
    new_name <- gsub("lr_target_prior_cor_filtered", "igraph_objects", name)
    igraph_objects_list[new_name] <- list(igraph_objects)
  }
  return(igraph_objects_list)
}

## analyze_network
# A function to analyze network properties, especially between receptors and targets of interest. Adapted from Jordan H. Whitlock.
analyze_network <- function(igraph_object, receptor, target, name){
  # calculate network properties
  diameter <- diameter(igraph_object)
  edge_number <- length(E(igraph_object))
  
  # Calculate degree centrality for receptor
  degree_centrality <- degree(igraph_object, v = receptor, mode = "out")
  
  # Calculate betweeness centrality for receptor
  betweenness_centrality <- betweenness(igraph_object, v = receptor)
  
  # Calculate closeness centrality for receptor
  closeness_centrality <- closeness(igraph_object, v = receptor)
  
  # all nodes attached to receptor
  neighbors_receptor <- neighbors(igraph_object, v = receptor, mode = "out")
  
  # receptor to target shortest path
  shortest_path <- shortest_paths(igraph_object,
                                  from = receptor,
                                  to = target,
                                  output = "both")
  
  
  
  all_shortest_path <- all_shortest_paths(igraph_object,
                                          from = receptor,
                                          to = target,
                                          mode = "out")
  
  dijkstra <- distances(igraph_object,
                        algorithm = "dijkstra")
  
  
  results <- list(id = name,
                  diameter = diameter,
                  edge_number = edge_number,
                  degree_centrality = degree_centrality,
                  betweenness_centrality = betweenness_centrality,
                  closeness_centrality = closeness_centrality,
                  neighbors_receptor = neighbors_receptor,
                  shortest_path = shortest_path,
                  all_shortest_path = all_shortest_path,
                  dijkstra = dijkstra)
  
  return(results)
}

## network_topology
# A wrapper function to allow for automatic application of the analyze_network function across multiple igraph objects.
network_topology <- function(object_list) {
  network_proporties_list <- list()
  for (name in names(object_list)) {
    igraph <- object_list[[name]]
    string <- unlist(strsplit(name, "_"))
    receptor <- string[2]
    target <- string[5]
    
    network_proporties_list[[name]] <- analyze_network(igraph = igraph, receptor = receptor, target = target, name = name)
  }
  return(network_proporties_list)
}






