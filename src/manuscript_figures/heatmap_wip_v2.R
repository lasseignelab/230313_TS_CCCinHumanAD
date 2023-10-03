library(circlize)
library(reshape2)
library(ComplexHeatmap)

geo_output <- readRDS("/data/user/tsoelter/projects/230313_TS_CCCinHumanAD/data/ccc/geo_multinichenet_output.rds")
gse_output <- readRDS("/data/user/tsoelter/projects/230313_TS_CCCinHumanAD/data/ccc/gse_multinichenet_output.rds")

de_genes_geo <- geo_output[["ligand_activities_targets_DEgenes"]][["de_genes_df"]]
de_genes_gse <- gse_output[["ligand_activities_targets_DEgenes"]][["de_genes_df"]]

de_genes_filt_geo <- de_genes_geo %>%
  filter(receiver %in% receiver_oi, contrast == "AD-CTRL")

de_genes_filt_gse <- de_genes_gse %>%
  filter(receiver %in% receiver_oi, contrast == "AD-CTRL")

geo_targets <- lr_target_prior_cor_filtered_geo$target %>% unique()

gse_targets <- lr_target_prior_cor_filtered_gse$target %>% unique()

de_genes_filt_geo <- de_genes_filt_geo %>%
  filter(gene %in% geo_targets) %>%
  select(gene, receiver, logFC) %>%
  unique() %>%
  mutate(dataset = "geo")

de_genes_filt_gse <- de_genes_filt_gse %>%
  filter(gene %in% gse_targets) %>%
  select(gene, receiver, logFC) %>%
  unique() %>%
  mutate(dataset = "gse")

de_genes_filt <- rbind(de_genes_filt_geo, de_genes_filt_gse)

# prep df
heatmap_matrix <- dcast(de_genes_filt, gene ~ receiver + dataset, value.var = "logFC") %>%
  column_to_rownames("gene")

heatmap_matrix[is.na(heatmap_matrix)] <- 0

# prep anno
meta <- as.data.frame(colnames(heatmap_matrix))
rownames(meta) <- meta$`colnames(heatmap_matrix)`

Receiver <- c("Excitatory.Neurons", "Excitatory.Neurons", "Inhibitory.Neurons", "Inhibitory.Neurons")
Dataset <- c("geo", "gse", "geo", "gse")
met <- cbind(Dataset, Receiver)
rownames(met) <- rownames(meta)
met <- as.data.frame(met)

#anno <- HeatmapAnnotation(
  #Receiver = c("Excitatory.Neurons", "Inhibitory.Neurons"),
  #Dataset = c("geo", "gse"),
  #col = list(Receiver = c("Excitatory.Neurons" = "slateblue3", "Inhibitory.Neurons" = "darkslategray3"))
#)



anno_cols <- list("Receiver" = c("Excitatory.Neurons" = "slateblue3", "Inhibitory.Neurons" = "darkslategray3"),
                  "Dataset" = c("geo" = "firebrick3", "gse" = "goldenrod3"))

anno <- HeatmapAnnotation(df = met, show_annotation_name = TRUE, col = anno_cols)

heatmap_matrix <- heatmap_matrix[, rownames(met), drop = FALSE]

mat <- as.matrix(heatmap_matrix)

h1 <- Heatmap(mat,
        top_annotation = anno,
        show_column_names = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        cluster_columns = FALSE)

anno2 <- rowAnnotation(df = met, show_annotation_name = TRUE, col = anno_cols)
h2 <- Heatmap(t(mat),
              right_annotation = anno2,
              show_column_names = TRUE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              show_row_names = FALSE)

####### only plot overlapping targets gex
targets <- intersect(geo_targets, gse_targets)

de_genes_filt_geo <- de_genes_filt_geo %>%
  filter(gene %in% targets) %>%
  select(gene, receiver, logFC) %>%
  unique() %>%
  mutate(dataset = "geo")

de_genes_filt_gse <- de_genes_filt_gse %>%
  filter(gene %in% targets) %>%
  select(gene, receiver, logFC) %>%
  unique() %>%
  mutate(dataset = "gse")

de_genes_filt <- rbind(de_genes_filt_geo, de_genes_filt_gse)

# prep df
heatmap_matrix <- dcast(de_genes_filt, gene ~ receiver + dataset, value.var = "logFC") %>%
  column_to_rownames("gene")

heatmap_matrix[is.na(heatmap_matrix)] <- 0

# prep anno
meta <- as.data.frame(colnames(heatmap_matrix))
rownames(meta) <- meta$`colnames(heatmap_matrix)`

Receiver <- c("Excitatory.Neurons", "Excitatory.Neurons", "Inhibitory.Neurons", "Inhibitory.Neurons")
Dataset <- c("geo", "gse", "geo", "gse")
met <- cbind(Dataset, Receiver)
rownames(met) <- rownames(meta)
met <- as.data.frame(met)

#anno <- HeatmapAnnotation(
#Receiver = c("Excitatory.Neurons", "Inhibitory.Neurons"),
#Dataset = c("geo", "gse"),
#col = list(Receiver = c("Excitatory.Neurons" = "slateblue3", "Inhibitory.Neurons" = "darkslategray3"))
#)



anno_cols <- list("Receiver" = c("Excitatory.Neurons" = "slateblue3", "Inhibitory.Neurons" = "darkslategray3"),
                  "Dataset" = c("geo" = "firebrick3", "gse" = "goldenrod3"))

anno <- HeatmapAnnotation(df = met, show_annotation_name = FALSE, col = anno_cols,
                          annotation_legend_param = list(direction = "horizontal", labels_gp = gpar("bold")))

anno <- HeatmapAnnotation(df = met, show_annotation_name = FALSE, col = anno_cols,
                          annotation_legend_param = list(Dataset = list(direction = "horizontal", labels_gp = gpar(fontface = "bold")),
                                                         Receiver = list(labels = c("Excitatory Neurons", "Inhibitory Neurons"), direction = "horizontal", labels_gp = gpar(fontface = "bold"))))
                            
                          

heatmap_matrix <- heatmap_matrix[, rownames(met), drop = FALSE]

mat <- as.matrix(heatmap_matrix)

h2 <- Heatmap(mat,
              top_annotation = anno,
              show_column_names = FALSE,
              show_row_dend = TRUE,
              show_column_dend = FALSE,
              cluster_columns = FALSE,
              row_title = "Overlapping Target Genes",
              row_title_gp = gpar(fontface = "bold"),
              row_names_gp = gpar(fontface = "bold"),
              heatmap_legend_param = list(direction = "horizontal", title = "Pseudobulk Expression",
                                          labels_gp = gpar(fontface = "bold")))


# plot only 3 targets that had overlapping l-r pairs
####### only plot overlapping targets gex
targets <- c("SMAD7", "LMO1", "CCND1")

de_genes_filt_geo <- de_genes_filt_geo %>%
  filter(gene %in% targets) %>%
  select(gene, receiver, logFC) %>%
  unique() %>%
  mutate(dataset = "geo")

de_genes_filt_gse <- de_genes_filt_gse %>%
  filter(gene %in% targets) %>%
  select(gene, receiver, logFC) %>%
  unique() %>%
  mutate(dataset = "gse")

de_genes_filt <- rbind(de_genes_filt_geo, de_genes_filt_gse)

# prep df
heatmap_matrix <- dcast(de_genes_filt, gene ~ receiver + dataset, value.var = "logFC") %>%
  column_to_rownames("gene")

heatmap_matrix[is.na(heatmap_matrix)] <- 0

# prep anno
meta <- as.data.frame(colnames(heatmap_matrix))
rownames(meta) <- meta$`colnames(heatmap_matrix)`

Receiver <- c("Excitatory.Neurons", "Excitatory.Neurons", "Inhibitory.Neurons", "Inhibitory.Neurons")
Dataset <- c("geo", "gse", "geo", "gse")
met <- cbind(Dataset, Receiver)
rownames(met) <- rownames(meta)
met <- as.data.frame(met)

#anno <- HeatmapAnnotation(
#Receiver = c("Excitatory.Neurons", "Inhibitory.Neurons"),
#Dataset = c("geo", "gse"),
#col = list(Receiver = c("Excitatory.Neurons" = "slateblue3", "Inhibitory.Neurons" = "darkslategray3"))
#)



anno_cols <- list("Receiver" = c("Excitatory.Neurons" = "slateblue3", "Inhibitory.Neurons" = "darkslategray3"),
                  "Dataset" = c("geo" = "firebrick3", "gse" = "goldenrod3"))

anno <- HeatmapAnnotation(df = met, show_annotation_name = FALSE, col = anno_cols,
                          annotation_legend_param = list(direction = "horizontal", labels_gp = gpar("bold")))

anno <- HeatmapAnnotation(df = met, show_annotation_name = FALSE, col = anno_cols,
                          annotation_legend_param = list(Dataset = list(direction = "horizontal", labels_gp = gpar(fontface = "bold", fontsize = 12)),
                                                         Receiver = list(labels = c("Excitatory Neurons", "Inhibitory Neurons"), direction = "horizontal", labels_gp = gpar(fontface = "bold", fontsize = 12))))



heatmap_matrix <- heatmap_matrix[, rownames(met), drop = FALSE]

mat <- as.matrix(heatmap_matrix)

h2 <- Heatmap(mat,
              top_annotation = anno,
              show_column_names = FALSE,
              show_row_dend = TRUE,
              show_column_dend = FALSE,
              cluster_columns = FALSE,
              row_title = "Overlapping Target Genes",
              row_title_gp = gpar(fontface = "bold"),
              row_names_gp = gpar(fontface = "bold"),
              heatmap_legend_param = list(direction = "horizontal", title = "Pseudobulk Expression",
                                          labels_gp = gpar(fontface = "bold")))


## heatmap of ligand activity of overlaps
geo_ligand_activity <- geo_output[["ligand_activities_targets_DEgenes"]][["ligand_activities"]]

geo_ligand_activity <- geo_ligand_activity %>%
  filter(contrast == "AD-CTRL", receiver %in% receiver_oi)

geo_overlaps <- lr_target_prior_cor_filtered_geo %>%
  filter(id_target %in% i5) %>%
  select(ligand, target, receiver, ligand_receptor_target)

geo_ligand_act <- inner_join(geo_overlaps, geo_ligand_activity, by = c("ligand", "target")) %>%
  mutate(dataset = "geo")

gse_ligand_activity <- gse_output[["ligand_activities_targets_DEgenes"]][["ligand_activities"]]

gse_ligand_activity <- gse_ligand_activity %>%
  filter(contrast == "AD-CTRL", receiver %in% receiver_oi)

gse_overlaps <- lr_target_prior_cor_filtered_gse %>%
  filter(id_target %in% i5) %>%
  select(ligand, target, receiver, ligand_receptor_target)

gse_ligand_act <- inner_join(gse_overlaps, gse_ligand_activity, by = c("ligand", "target")) %>%
  mutate(dataset = "gse")

ligand_activity <- rbind(geo_ligand_act, gse_ligand_act)

ligand_activity_filt <- ligand_activity %>%
  select(ligand_receptor_target, activity_scaled, dataset)

heatmap_matrix <- ligand_activity_filt %>%
  pivot_wider(names_from = dataset, values_from = activity_scaled) %>%
  column_to_rownames("ligand_receptor_target")

meta <- as.data.frame(colnames(heatmap_matrix))
rownames(meta) <- meta$`colnames(heatmap_matrix)`
colnames(meta) <- "Dataset"

anno_cols <- list("Dataset" = c("geo" = "firebrick3", "gse" = "goldenrod3"))

anno <- HeatmapAnnotation(df = meta, show_annotation_name = TRUE, col = anno_cols)

heatmap_matrix <- heatmap_matrix[, rownames(meta), drop = FALSE]

mat <- as.matrix(heatmap_matrix)

cols <- colorRamp2(c(-4,-2, 0, 2, 4), c("blue", "purple", "white", "yellow", "red"))

h4 <- Heatmap(mat,
              col = cols,
              top_annotation = anno,
              show_column_names = FALSE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              cluster_columns = FALSE,
              cluster_rows = TRUE)

