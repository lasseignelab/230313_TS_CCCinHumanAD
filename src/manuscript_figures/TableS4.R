
geo_t <- lr_target_prior_cor_filtered_geo %>%
  filter(id %in% i4) %>%
  select(id, sender, receiver, target, ligand_receptor) %>%
  unique()

geo_t$dataset <- "geo"

gse_t <- lr_target_prior_cor_filtered_gse %>%
  filter(id %in% i4) %>%
  select(id, sender, receiver, target, ligand_receptor) %>%
  unique()

gse_t$dataset <- "gse"


overlap <- bind_rows(geo_t, gse_t)

shared_targets <- overlap %>%
  pivot_wider(
    names_from = dataset,
    values_from = target,
    names_prefix = "targets_"
  ) %>%
  unnest_wider(c(starts_with("targets_"))) %>%
  select(ligand_receptor, sender, receiver, starts_with("targets_"))

install.packages("openxlsx")
library(openxlsx)

write.xlsx(shared_targets, here("results", "final_outputs", "comparison", "overlappingLR_targets.xlsx"))
