library(scDesign3)
library(SingleCellExperiment)

sce <- readRDS('/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/E10_E14_2E12.rds')
sce <- readRDS('/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/E12_CS7_2E14.rds')
sce <- readRDS('/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/CS7_CS9_2CS8.rds')
sce <- readRDS("/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/CS9_CS11_2CS10.rds")
sce <- readRDS("/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/CS10_CS12_2CS11.rds")

set.seed(123)
re_sim_s <- Sys.time()
simu <- scdesign3(
  sce = sce,
  assay_use = "counts",
  celltype = "merge_type",
  pseudotime = NULL,
  other_covariates = NULL,
  mu_formula = "1",
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 18,  #8 for E12\E14\CS8, 18 for CS11, 23 for CS10
  corr_formula = "1"
)
re_sim_e <- Sys.time()

simu_sce <- SingleCellExperiment(list(counts = simu$new_count), colData = simu$new_covariate,
                                 rowData = row.names(list(counts = simu$new_count)))

saveRDS(simu_sce, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/E12f_scDesign3.rds")
saveRDS(simu_sce, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/E14f_scDesign3.rds")
saveRDS(simu_sce, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/CS8f_scDesign3.rds")
saveRDS(simu_sce, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/CS10f_scDesign3.rds")
saveRDS(simu_sce, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/CS11f_scDesign3.rds")

set.seed(123)
compare_figure <- plot_reduceddim(ref_sce = sce, 
                                  sce_list = list(simu_sce), 
                                  name_vec = c("Reference", "scDesign3"),
                                  assay_use = "counts", 
                                  if_plot = TRUE, 
                                  color_by = "merge_type", 
                                  n_pc = 20)
plot(compare_figure$p_umap)