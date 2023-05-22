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
prep_liana <- function(object_list, file_path) {
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
           path = paste0(here(file_path, name, "/")),
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
           path = paste0(here(file_path, name, "/")),
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

## prep_nichenet
# A function which prepares seurat objects for NicheNet analysis by creating cell type aggregate columns and assigning them as active.ident. It also plots and saves UMAPs showing the split by condition and cell type.
prep_nichenet <- function(object_list, file_path) {
  objects <- tibble::lst()
  for (name in names(object_list)) {
    object <- get(name)
    print(name)
    # Check number of cells per cell type and condition ----------
    print("Checking number of cells per cell type and condition")
    print(table(object@active.ident, object@meta.data$orig.ident))
    # Change cell type names based on condition origin
    print("Adding celltype_aggregate column to metadata")
    object@meta.data$celltype_aggregate <- paste(object@active.ident,
                                                 object@meta.data$orig.ident,
                                                 sep = "_")
    print("Plot cell_type aggregate UMAP and save to intermediate_outputs")
    umap <- DimPlot(object, group.by = "celltype_aggregate")
    plot(umap)
    ggsave(filename = "UMAP_celltype_aggregate.png",
           path = paste0(here(file_path, name, "/ccc/")),
           plot = umap)
    # change metadata column and set as identity ----------
    print("Setting celltype_aggregate column as identity of the object")
    celltype_id <- "celltype_aggregate"
    object <- SetIdent(object, value = object[[celltype_id]])
    objects[name] <- object
  }
  return(objects)
}

## diff_nichenet
# A function which performs differential gene expression analysis and cross-references the NicheNet lr_network to identify DE ligands and receptors
# The original code is from the vignette, but has been made more legible and functional
diff_nichenet <- function(object, niches, expression_pct, lr_network, assay_oi = "RNA") {
  # Differential analysis in sender cells
  DE_sender <- calculate_niche_de(seurat_obj = object %>% 
                                    subset(features = lr_network$ligand %>%
                                             unique()),
                                  niches = niches,
                                  type = "sender",
                                  assay_oi = assay_oi)
  # Differential analysis in receiver cells
  DE_receiver <- calculate_niche_de(seurat_obj = object %>% 
                                      subset(features = lr_network$receptor %>% 
                                               unique()),
                                    niches = niches,
                                    type = "receiver",
                                    assay_oi = assay_oi)
  # Filter log2FC for sender
  DE_sender <- DE_sender %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf,
                               max(avg_log2FC[is.finite(avg_log2FC)]),
                               ifelse(avg_log2FC == -Inf,
                                      min(avg_log2FC[is.finite(avg_log2FC)]),
                                      avg_log2FC)))
  # Filter log2FC for receiver
  DE_receiver <- DE_receiver %>%
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf,
                               max(avg_log2FC[is.finite(avg_log2FC)]),
                               ifelse(avg_log2FC == -Inf,
                                      min(avg_log2FC[is.finite(avg_log2FC)]),
                                      avg_log2FC)))
  # Process sender and receiver niches
  DE_sender_processed <- process_niche_de(DE_table = DE_sender,
                                          niches = niches,
                                          expression_pct = expression_pct,
                                          type = "sender")
  DE_receiver_processed <- process_niche_de(DE_table = DE_receiver,
                                            niches = niches,
                                            expression_pct = expression_pct,
                                            type = "receiver")
  # Combine sender and receiver
  DE_sender_receiver <- combine_sender_receiver_de(DE_sender_processed,
                                                   DE_receiver_processed,
                                                   lr_network,
                                                   specificity_score = "min_lfc")
  return(DE_sender_receiver)
}

## make_mock_spatial
# A function which makes objects needed for downstream NicheNet analyses as a list called mock_spatial_data_list, which can be unlisted for future use. 
make_mock_spatial <- function(include_spatial_info_sender,
                              include_spatial_info_receiver,
                              niches,
                              expression_pct,
                              specificity_score_spatial = "lfc") {
  mock_spatial_data_list <- list()
  if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE) {
    spatial_info <- tibble(celltype_region_oi = NA,
                           celltype_other_region = NA) %>%
      mutate(niche = niches %>%
               names() %>%
               head(1),
             celltype_type = "sender")
  }
  # sender spatial info
  if(include_spatial_info_sender == TRUE) {
    sender_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>% 
                                                subset(features = lr_network$ligand %>%
                                                         unique()),
                                              spatial_info = spatial_info %>%
                                                filter(celltype_type == "sender"))
    sender_spatial_DE_processed <- process_spatial_de(DE_table = sender_spatial_DE,
                                                      type = "sender",
                                                      lr_network = lr_network,
                                                      expression_pct = 0.1,
                                                      specificity_score = specificity_score_spatial)
    # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
    sender_spatial_DE_others <- get_non_spatial_de(niches = niches,
                                                   spatial_info = spatial_info,
                                                   type = "sender",
                                                   lr_network = lr_network)
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
      bind_rows(sender_spatial_DE_others)
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  } else {
    # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
    sender_spatial_DE_processed <- get_non_spatial_de(niches = niches,
                                                      spatial_info = spatial_info,
                                                      type = "sender",
                                                      lr_network = lr_network)
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  }
  # receiver spatial info 
  if(include_spatial_info_receiver == TRUE) {
    receiver_spatial_DE <- calculate_spatial_DE(seurat_obj = seurat_obj %>%
                                                  subset(features = lr_network$receptor %>%
                                                           unique()),
                                                spatial_info = spatial_info %>%
                                                  filter(celltype_type == "receiver"))
    receiver_spatial_DE_processed <- process_spatial_de(DE_table = receiver_spatial_DE,
                                                        type = "receiver",
                                                        lr_network = lr_network,
                                                        expression_pct = 0.1,
                                                        specificity_score = specificity_score_spatial)
    # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
    receiver_spatial_DE_others <- get_non_spatial_de(niches = niches,
                                                     spatial_info = spatial_info,
                                                     type = "receiver",
                                                     lr_network = lr_network)
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
      bind_rows(receiver_spatial_DE_others)
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  } else {
    # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
    receiver_spatial_DE_processed <- get_non_spatial_de(niches = niches,
                                                        spatial_info = spatial_info,
                                                        type = "receiver",
                                                        lr_network = lr_network)
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }
  mock_spatial_data_list <- tibble::lst(sender_spatial_DE_processed, receiver_spatial_DE_processed)
  return(mock_spatial_data_list)
}

## calculate_ligand_activity
# A function that creates an object called ligand_activities_targets needed for downstream NicheNet analyses.
calculate_ligand_activity <- function(object, niches, top_n_targets, lfc_cutoff, assay_oi = "RNA") {
  DE_receiver_targets <- calculate_niche_de_targets(seurat_obj = object, 
                                                    niches = niches,
                                                    lfc_cutoff = lfc_cutoff,
                                                    expression_pct = 0.1,
                                                    assay_oi = assay_oi) 
  DE_receiver_processed_targets <- process_receiver_target_de(DE_receiver_targets = DE_receiver_targets,
                                                              niches = niches,
                                                              expression_pct = 0.1,
                                                              specificity_score = "min_lfc")
  background <- DE_receiver_processed_targets %>% pull(target) %>% unique()
  # Determine gene set of niche 1 ----------------------------------------------
  geneset_niche1 <- DE_receiver_processed_targets %>% 
    filter(receiver == niches[[1]]$receiver &
             target_score >= lfc_cutoff &
             target_significant == 1 &
             target_present == 1) %>%
    pull(target) %>%
    unique()
  # Determine gene set of niche 2 ----------------------------------------------
  geneset_niche2 <- DE_receiver_processed_targets %>%
    filter(receiver == niches[[2]]$receiver &
             target_score >= lfc_cutoff &
             target_significant == 1 &
             target_present == 1) %>%
    pull(target) %>%
    unique()
  # check which genes are excluded from the gene sets --------------------------
  geneset_niche1 %>% setdiff(rownames(ligand_target_matrix)) 
  geneset_niche2 %>% setdiff(rownames(ligand_target_matrix)) 
  length1 <- length(geneset_niche1)
  print(paste0("geneset_niche1 has ", length1, " genes"))
  length2 <- length(geneset_niche2)
  print(paste0("geneset_niche2 has ", length2, " genes"))
  # Make niche gene set list ---------------------------------------------------
  niche_geneset_list <- list(
    "AD_niche" = list(
      "receiver" = niches[[1]]$receiver,
      "geneset" = geneset_niche1,
      "background" = background),
    "CTRL_niche" = list(
      "receiver" = niches[[2]]$receiver,
      "geneset" = geneset_niche2,
      "background" = background)
  )
  # Get ligand-target interactions and activity --------------------------------
  ligand_activities_targets <- get_ligand_activities_targets(niche_geneset_list = niche_geneset_list,
                                                             ligand_target_matrix = ligand_target_matrix,
                                                             top_n_target = top_n_targets)
  outs <- tibble::lst(DE_receiver_processed_targets, ligand_activities_targets)
  return(outs)
}

## calculate_scaled_gex
# A function which uses previously generated outputs and a seurat object as input and outputs a list of expression tables of ligands, receptors, and targets, which needs to be unlisted after running this function. 
calculate_scaled_gex <- function(object,
                                 lr_network,
                                 ligand_activities_targets,
                                 niches = user_niches,
                                 assay_oi = "RNA",
                                 expression_pct = 0.1) {
  # get ligands, receptors, and targets from ligand_activities_targets and lr_networks
  features_oi <- union(lr_network$ligand, lr_network$receptor) %>%
    union(ligand_activities_targets$target) %>%
    setdiff(NA)
  dotplot <- suppressWarnings(Seurat::DotPlot(object %>%
                                                subset(idents = niches %>%
                                                         unlist() %>%
                                                         unique()),
                                              features = features_oi,
                                              assay = assay_oi))
  # create expression table
  exprs_tbl <- dotplot$data %>% as_tibble()
  exprs_tbl <- exprs_tbl %>%
    rename(celltype = id,
           gene = features.plot,
           expression = avg.exp,
           expression_scaled = avg.exp.scaled,
           fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>%
    as_tibble() %>%
    select(celltype, gene, expression, expression_scaled, fraction) %>%
    distinct() %>%
    arrange(gene) %>% 
    mutate(gene = as.character(gene))
  # create ligand expression df
  exprs_tbl_ligand <- exprs_tbl %>%
    filter(gene %in% lr_network$ligand) %>%
    rename(sender = celltype,
           ligand = gene,
           ligand_expression = expression,
           ligand_expression_scaled = expression_scaled,
           ligand_fraction = fraction)
  # create receptor expression df
  exprs_tbl_receptor <- exprs_tbl %>%
    filter(gene %in% lr_network$receptor) %>%
    rename(receiver = celltype,
           receptor = gene,
           receptor_expression = expression,
           receptor_expression_scaled = expression_scaled,
           receptor_fraction = fraction)
  # create target expression df
  exprs_tbl_target <- exprs_tbl %>%
    filter(gene %in% ligand_activities_targets$target) %>%
    rename(receiver = celltype,
           target = gene,
           target_expression = expression,
           target_expression_scaled = expression_scaled,
           target_fraction = fraction)
  # create final ligand expression df
  exprs_tbl_ligand <- exprs_tbl_ligand %>%
    mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>%
    mutate(ligand_fraction_adapted = ligand_fraction) %>%
    mutate_cond(ligand_fraction >= expression_pct,
                ligand_fraction_adapted = expression_pct) %>%
    mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  # create final receptor expression df
  exprs_tbl_receptor <- exprs_tbl_receptor %>%
    mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>%
    mutate(receptor_fraction_adapted = receptor_fraction) %>%
    mutate_cond(receptor_fraction >= expression_pct,
                receptor_fraction_adapted = expression_pct) %>%
    mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  # return list with all relevant objects
  exprs_tbl_list <- tibble::lst(exprs_tbl_receptor, exprs_tbl_ligand, exprs_tbl_target)
  return(exprs_tbl_list)
}

## score_interactions
# A function which uses expression tables, DE information between niches, and the lr network to score the l-r interactions based on their expression strength of the receptor
score_interactions <- function(lr_network, exprs_tbl_ligand, exprs_tbl_receptor, DE_sender_receiver) {
  # Combine expression tables for ligands, receptors, and DE between sender and receiver niches 
  exprs_sender_receiver <- lr_network %>% 
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>%
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>%
    inner_join(DE_sender_receiver %>%
                 distinct(niche,
                          sender,
                          receiver)
    )
  # score interactions
  ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>%
    group_by(ligand, receiver) %>%
    mutate(rank_receptor_expression = dense_rank(receptor_expression),
           rank_receptor_fraction = dense_rank(receptor_fraction)) %>%
    mutate(ligand_scaled_receptor_expression_fraction = 0.5*((rank_receptor_fraction)) + 
             ((rank_receptor_expression / max(rank_receptor_expression))) ) %>%
    distinct(ligand,
             receptor,
             receiver,
             ligand_scaled_receptor_expression_fraction,
             bonafide) %>%
    distinct() %>%
    ungroup()
  # return object
  return(ligand_scaled_receptor_expression_fraction_df)
}

## nichenet_wrapper
# A wrapper function that includes all functions to run differential nichenet and generate the output matrix which includes all ligands, receptors, and target genes between sender and receiver cell types.
# The function has the following capabilities:
#   * cross-reference LIANA ligands to NicheNet ligands
#   * adding mock spatial data
#   * calculate ligand activities and infer active ligand-target links
#   * calculate the scaled expression of ligands, receptors, and targets across cell types of interest
#   * scoring the l-r interctions based on their expression strength of the receptor
nichenet_wrapper <- function(seurat_obj, user_niches, lr_network, ligands, file_path, file_name) {
  # DE analysis between niches --------------------
  DE_sender_receiver <- diff_nichenet(object = seurat_obj,
                                      niches = user_niches,
                                      expression_pct = 0.10,
                                      lr_network = lr_network)
  print("Differential gene expression analysis finished")
  # Filter ligands by LIANA ligands 
  DE_sender_receiver_LIANA <- subset(DE_sender_receiver, ligand %in% ligands)
  # Create mock spatial data --------------------
  mock_spatial_data_list <- make_mock_spatial(include_spatial_info_sender = FALSE,
                                              include_spatial_info_receiver = FALSE,
                                              niches = user_niches)
  print("Mock spatial data created")
  # unlist mock spatial data list to get objects needed --------------------
  list2env(mock_spatial_data_list, environment())
  # calculate ligand activity --------------------
  ligand_activity_list <- calculate_ligand_activity(seurat_obj,
                                                    niches = user_niches,
                                                    top_n_targets = 250,
                                                    lfc_cutoff = 0.15)
  print("Calculated ligand activity")
  # unlist ligand activity list list to get objects needed -------------------
  list2env(ligand_activity_list, environment())
  # calculate scaled expression --------------------
  exprs_tbl_list <- calculate_scaled_gex(seurat_obj,
                                         lr_network = lr_network,
                                         ligand_activities_targets = ligand_activities_targets)
  print("Calculated scaled expression")
  # unlist expression table list to get objects needed --------------------
  list2env(exprs_tbl_list, environment())
  # score interactions --------------------
  ligand_scaled_receptor_expression_fraction_df <-
    score_interactions(lr_network = lr_network,
                       exprs_tbl_ligand = exprs_tbl_ligand,
                       exprs_tbl_receptor = exprs_tbl_receptor,
                       DE_sender_receiver = DE_sender_receiver_LIANA)
  print("Scored interactions")
  # combine all outputs into one large list --------------------
  output <- list(DE_sender_receiver = DE_sender_receiver_LIANA,
                 ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
                 sender_spatial_DE_processed = sender_spatial_DE_processed,
                 receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                 ligand_activities_targets = ligand_activities_targets,
                 DE_receiver_processed_targets = DE_receiver_processed_targets,
                 exprs_tbl_ligand = exprs_tbl_ligand,
                 exprs_tbl_receptor = exprs_tbl_receptor,
                 exprs_tbl_target = exprs_tbl_target)
  print("Created output")
  # save output
  saveRDS(output, file = here(paste0(file_path, file_name)
  )
  )
  print("Saved output")
}

## prioritize_interactions
# A function which inputs a list of NicheNet outputs and prioritizes interactions. The resulting prioritization tables output is saved to respective directory and returned as a list for downstream analyses.
prioritize_interactions <- function(output_list, file_path, weights = prioritizing_weights) {
  prioritization_tables_ls <- tibble::lst()
  for(i in names(output_list)) {
    output <- get(i)
    name <- sub("output", "prioritization_tbl", i)
    path <- sub("_prioritization_tbl", "/ccc/", name)
    prioritization_tables <- get_prioritization_tables(output, weights)
    saveRDS(prioritization_tables, file = here(file_path, path, "prioritization_tables.rds"))
    prioritization_tables_ls[name] <- list(prioritization_tables)
  }
  return(prioritization_tables_ls)
}


## make_nichenet_plot
# A function which prioritizes ligands based on prioritization score for plotting and generate the exprs_activity_target_plot basic for NicheNet.
make_nichenet_plot <- function(prioritization_tables, output, receiver_oi, lfc_cutoff) {
  top_ligand_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    select(niche,
           sender,
           receiver,
           ligand,
           receptor,
           prioritization_score
    ) %>%
    group_by(ligand) %>%
    top_n(1, prioritization_score) %>%
    ungroup() %>%
    select(ligand,
           receptor,
           niche
    ) %>%
    rename(top_niche = niche)
  top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    select(niche,
           sender,
           receiver,
           ligand,
           receptor,
           prioritization_score
    ) %>%
    group_by(ligand,
             receptor) %>%
    top_n(1, prioritization_score) %>%
    ungroup() %>%
    select(ligand,
           receptor,
           niche
    ) %>%
    rename(top_niche = niche)
  ligand_prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    select(niche,
           sender,
           receiver,
           ligand,
           prioritization_score
    ) %>%
    group_by(ligand,
             niche
    ) %>%
    top_n(1, prioritization_score) %>%
    ungroup() %>%
    distinct() %>%
    inner_join(top_ligand_niche_df) %>%
    filter(niche == top_niche) %>%
    group_by(niche) %>%
    top_n(50, prioritization_score) %>% # get the top50 ligands per niche
    ungroup() 
  
  filtered_ligands <- ligand_prioritized_tbl_oi %>%
    filter(receiver == receiver_oi) %>%
    top_n(20, prioritization_score) %>%
    pull(ligand) %>%
    unique()
  
  prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    filter(ligand %in% filtered_ligands) %>%
    select(niche,
           sender,
           receiver,
           ligand,
           receptor,
           ligand_receptor,
           prioritization_score
    ) %>% distinct() %>%
    inner_join(top_ligand_receptor_niche_df) %>%
    group_by(ligand) %>%
    filter(receiver == receiver_oi) %>%
    top_n(2, prioritization_score) %>%
    ungroup() 
  
  exprs_activity_target_plot <-
    make_ligand_activity_target_exprs_plot(receiver_oi,
                                           prioritized_tbl_oi,
                                           prioritization_tables$prioritization_tbl_ligand_receptor,
                                           prioritization_tables$prioritization_tbl_ligand_target,
                                           output$exprs_tbl_ligand,
                                           output$exprs_tbl_target,
                                           lfc_cutoff,
                                           ligand_target_matrix,
                                           plot_legend = FALSE,
                                           heights = NULL,
                                           widths = NULL)
  lfc_plot <- make_ligand_receptor_lfc_plot(receiver_oi,
                                            prioritized_tbl_oi,
                                            prioritization_tables$prioritization_tbl_ligand_receptor,
                                            plot_legend = FALSE,
                                            heights = NULL,
                                            widths = NULL)
  nichenet_plot <- tibble::lst(exprs_activity_target_plot, lfc_plot)
  return(nichenet_plot)
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
bubbleplot <- function(fea_result, region, file_path){
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
  ggsave(filename = "bubbleplot_pathways.png",
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
  combined_fea <- rbind(up_fea_filt, down_fea_filt)
  fea_result <- list(downregulated = down_fea_filt, upregulated = up_fea_filt, combined = combined_fea)
  return(fea_result)
}

## pathway_analysis
# A wrapper function to prepare inputs for pathway analysis, perform pathway analysis using gprofiler2 on up/down regulated target genes from nichenet
pathway_analysis <- function(nichenet_output, receiver) {
  # subset output for df with target expression info -----
  target_expression <- nichenet_output$exprs_tbl_target
  # filter expression df for necessary columns and cell type as well as add up/down gex info -----
  neuro_targets_direction <- target_expression %>%
    select(receiver, target, target_expression_scaled) %>%
    filter(receiver == "Excitatory Neurons_AD") %>%
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
compile_gex <- function(object, receiver, condition, genes_ppi){
  target_all <- subset(object, 
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