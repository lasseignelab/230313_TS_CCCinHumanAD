ex_neurons_nichenet_output <- readRDS("/data/user/tsoelter/projects/230313_TS_CCCinHumanAD/results/intermediate_outputs/vine_cortex/ccc/ex_neurons_nichenet_output.rds")

ex_neurons_8v8_ctx_output <- readRDS("/data/user/tsoelter/projects/230313_TS_CCCinHumanAD/data/test/ex_neurons_8v8_ctx_output.rds")

ligands_vine_c <- ex_neurons_nichenet_output$exprs_tbl_ligand$ligand %>% unique()

ligands_8v8 <- ex_neurons_8v8_ctx_output$exprs_tbl_ligand$ligand %>% unique()

ligs <- intersect(ligands_vine_c, ligands_8v8)

receptors_vine_c <- ex_neurons_nichenet_output$exprs_tbl_receptor$receptor %>% unique()

receptors_8v8 <- ex_neurons_8v8_ctx_output$exprs_tbl_receptor$receptor %>% unique()

recs <- intersect(receptors_vine_c, receptors_8v8)


targets_vine_c <- ex_neurons_nichenet_output$exprs_tbl_target$target %>% unique()

targets_8v8 <- ex_neurons_8v8_ctx_output$exprs_tbl_target$target %>% unique()

targs <- intersect(targets_vine_c, targets_8v8)
