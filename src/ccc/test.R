make_mock_spatial <- function(include_spatial_info_sender,
                              include_spatial_info_receiver,
                              specificity_score_spatial = "lfc") {
  if (include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE) {
    spatial_info <- tibble(
      celltype_region_oi = NA,
      celltype_other_region = NA
    ) %>%
      mutate(
        niche = niches %>%
          names() %>%
          head(1),
        celltype_type = "sender"
      )
  }
  # sender spatial info
  if (include_spatial_info_sender == TRUE) {
    sender_spatial_DE <- calculate_spatial_DE(
      seurat_obj = seurat_obj %>%
        subset(features = lr_network$ligand %>%
          unique()),
      spatial_info = spatial_info %>%
        filter(celltype_type == "sender")
    )
    sender_spatial_DE_processed <- process_spatial_de(
      DE_table = sender_spatial_DE,
      type = "sender",
      lr_network = lr_network,
      expression_pct = expression_pct,
      specificity_score = specificity_score_spatial
    )
    # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
    sender_spatial_DE_others <- get_non_spatial_de(
      niches = niches,
      spatial_info = spatial_info,
      type = "sender",
      lr_network = lr_network
    )
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
      bind_rows(sender_spatial_DE_others)
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  } else {
    # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
    sender_spatial_DE_processed <- get_non_spatial_de(
      niches = niches,
      spatial_info = spatial_info,
      type = "sender",
      lr_network = lr_network
    )
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  }
  # receiver spatial info
  if (include_spatial_info_receiver == TRUE) {
    receiver_spatial_DE <- calculate_spatial_DE(
      seurat_obj = seurat_obj %>%
        subset(features = lr_network$receptor %>%
          unique()),
      spatial_info = spatial_info %>%
        filter(celltype_type == "receiver")
    )
    receiver_spatial_DE_processed <- process_spatial_de(
      DE_table = receiver_spatial_DE,
      type = "receiver",
      lr_network = lr_network,
      expression_pct = expression_pct,
      specificity_score = specificity_score_spatial
    )
    # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
    receiver_spatial_DE_others <- get_non_spatial_de(
      niches = niches,
      spatial_info = spatial_info,
      type = "receiver",
      lr_network = lr_network
    )
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
      bind_rows(receiver_spatial_DE_others)
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  } else {
    # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
    receiver_spatial_DE_processed <- get_non_spatial_de(
      niches = niches,
      spatial_info = spatial_info,
      type = "receiver",
      lr_network = lr_network
    )
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }
  return(c(sender_spatial_DE_processed, receiver_spatial_DE_processed))
}






style_file(here(
  "src",
  "ccc",
  "test.R"
))
