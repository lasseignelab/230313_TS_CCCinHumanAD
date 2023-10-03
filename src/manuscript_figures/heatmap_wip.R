library(ComplexHeatmap)

de_genes <- output[["ligand_activities_targets_DEgenes"]][["de_genes_df"]]

de_genes_filt <- de_genes %>%
  filter(receiver %in% receiver_oi, contrast == "AD-CTRL")

geo_targets <- lr_target_prior_cor_filtered_geo$target %>% unique()

gse_targets <- lr_target_prior_cor_filtered_gse$target %>% unique()

targets <- c(geo_targets,gse_targets) %>% unique()

de_genes_filt <- de_genes_filt %>%
  filter(gene %in% targets) %>%
  select(gene, receiver, logFC) %>%
  unique()

de_genes_filt_ex <- de_genes_filt %>%
  filter(receiver == "Inhibitory.Neurons")

# prep df
heatmap_matrix <- de_genes_filt %>%
  pivot_wider(names_from = receiver, values_from = logFC, values_fill = 0) %>%
  column_to_rownames("gene")

# prep anno
meta <- as.data.frame(colnames(heatmap_matrix))
rownames(meta) <- meta$`colnames(heatmap_matrix)`
colnames(meta) <- "Receiver"

anno_cols <- list("Receiver" = c("Excitatory.Neurons" = "slateblue3", "Inhibitory.Neurons" = "darkslategray3"))

anno <- HeatmapAnnotation(df = meta, show_annotation_name = TRUE, col = anno_cols)

heatmap_matrix <- heatmap_matrix[, rownames(meta), drop = FALSE]

mat <- as.matrix(heatmap_matrix)

Heatmap(mat,
        top_annotation = anno,
        cluster_rows = TRUE,
        show_column_names = FALSE,
        show_row_dend = FALSE)



### testing
h1 <- Heatmap(mat,
              cluster_rows = TRUE,
              show_column_dend = FALSE,
              show_column_names = FALSE) 

h2 <- Heatmap(mat,
              cluster_rows = TRUE,
              show_column_dend = FALSE,
              show_column_names = FALSE) 

h2 <- Heatmap(t(mat),
              cluster_columns = TRUE,
              show_column_dend = TRUE,
              show_column_names = TRUE,
              show_row_names = FALSE) 

h3 <- Heatmap(mat,
              top_annotation = anno,
              cluster_rows = TRUE,
              show_column_names = FALSE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              row_names_gp = gpar(fontsize = 10))
