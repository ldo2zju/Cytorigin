library(splatter)
library(SingleCellExperiment)

sce <- readRDS('/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/E10_E14_2E12.rds')
sce <- readRDS('/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/E12_CS7_2E14.rds')
sce <- readRDS('/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/CS7_CS9_2CS8.rds')
sce <- readRDS("/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/CS9_CS11_2CS10.rds")
sce <- readRDS("/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/source_matrix/CS10_CS12_2CS11.rds")

re_splEs_s <- Sys.time()
params <- splatEstimate(sce)
re_splEs_e <- Sys.time()

porp_count <- data.frame(table(colData(sce)$merge_type)/length(colData(sce)$merge_type))

set.seed(123)
re_sim_s <- Sys.time()
sim.groups <- splatSimulate(
  params = params,
  group.prob = porp_count$Freq,
  method = "groups",
  verbose = FALSE
)
re_sim_e <- Sys.time()
rowData(sim.groups) = rownames(sce)

saveRDS(sim.groups, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/E12f_splatter.rds")
saveRDS(sim.groups, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/E14f_splatter.rds")
saveRDS(sim.groups, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/CS8f_splatter.rds")
saveRDS(sim.groups, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/CS10f_splatter.rds")
saveRDS(sim.groups, file = "/slurm/home/yrd/liaolab/caohaoxue/embryo_work/benchmark/simulated_matrix/CS11f_splatter.rds")

sim.groups <- logNormCounts(sim.groups)
sim.groups <- runPCA(sim.groups)
plotPCA(sim.groups, colour_by = "Group")

