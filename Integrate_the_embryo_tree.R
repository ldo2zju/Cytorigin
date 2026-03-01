library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(RColorBrewer)
#library(plot1cell)
library(anndata)
library(Matrix)
library(reticulate)
use_python("/slurm/home/yrd/liaolab/wangtianhao/anaconda3/envs/sc3.8.5/bin/python", required = T)
py_config()

g <- read_h5ad("~/data/embryo/combined_new_inner_2.h5ad")
obj <- CreateSeuratObject(t(g$X), meta.data = g$obs)

obj$celltype <- as.character(obj$celltype)
obj$developmental.system <- as.character(obj$developmental.system)
obj$subCT <- as.character(obj$subCT)
obj@meta.data[is.na(obj$celltype),]$celltype <- "NA"
obj@meta.data[is.na(obj$developmental.system),]$developmental.system <- "NA"
obj@meta.data[is.na(obj$subCT),]$subCT <- "NA"
#obj@meta.data[obj$stage=="CS8_CytOrgin",]$celltype <- obj@meta.data[obj$stage=="CS8_CytOrgin",]$cell_type
meta <- read.csv("/slurm/home/yrd/liaolab/wangtianhao/data/embryo/CS12_16_meta.csv")
rownames(meta) <- meta$X
meta <- meta[colnames(obj[,obj$stage %in% c("CS12","CS13-14","CS15-16")]),]
obj@meta.data[obj$stage %in% c("CS12","CS13-14","CS15-16"),]$celltype <- meta[,"annotation"]
obj@meta.data[obj$stage=="CS8_CytOrgin",]$celltype <- as.character(obj@meta.data[obj$stage=="CS8_CytOrgin",]$cell_type)

# Filter miscellaneous cells, high-MT cells, and undefined cells
obj <- obj[,!(obj$developmental.system=="miscellaneous")]
obj <- obj[,!(obj$celltype %in% c("Other_MT_high", "Other_undefined", "Trophoblast_MT_high", "EPI.PrE.INT", "EPI.PrE.TE"))]
obj <- obj[,!(obj$celltype %in% c("Unknown","Ambiguous","undefined immature neuron-2","undefined-1","undefined-2","undefined-3","undefined-4","undefined-5","undefined (EN1)","undefined (GADD45G)","undefined immature neuron-1"))]

# Rename synonyms
obj@meta.data[obj$celltype %in% c("early_TE","medium_TE","late_TE","EarlyTE"),]$celltype <- "TE"
obj@meta.data[obj$celltype=="Epiblast",]$celltype <- "EPI"
obj@meta.data[obj$celltype %in% c("EightCells","L8C"),]$celltype <- "L08C"
obj@meta.data[obj$celltype %in% c("Somites","Somite"),]$celltype <- "somite"
obj@meta.data[obj$celltype %in% c("Cardiac Myocyte"),]$celltype <- "Cardiomyocyte"
obj@meta.data[obj$celltype %in% c("Neuron"),]$celltype <- "neuron"
obj@meta.data[obj$celltype %in% c("Endoderm"),]$celltype <- "endoderm"
obj@meta.data[obj$celltype %in% c("AxMes"),]$celltype <- "Axial Mesoderm"
obj@meta.data[obj$celltype %in% c("Fibroblast"),]$celltype <- "fibroblast"
obj@meta.data[obj$celltype %in% c("Epithelium"),]$celltype <- "epithelium"
obj@meta.data[obj$celltype %in% c("PriS"),]$celltype <- "Primitive Streak"
obj@meta.data[obj$celltype %in% c("Spinal cord"),]$celltype <- "spinal cord"
obj@meta.data[obj$celltype %in% c("Endothelium"),]$celltype <- "endothelium"
obj@meta.data[obj$celltype %in% c("YS.Endoderm"),]$celltype <- "YS endoderm"
obj@meta.data[obj$celltype %in% c("AdvMes"),]$celltype <- "Advanced Mesoderm"
obj@meta.data[obj$celltype %in% c("Head mesoderm"),]$celltype <- "head mesoderm"
obj@meta.data[obj$celltype %in% c("Splanchnic LPM"),]$celltype <- "splanchnic LPM"
obj@meta.data[obj$celltype %in% c("Epidermis","Epidermis-female"),]$celltype <- "epidermis"
obj@meta.data[obj$celltype %in% c("Neural Crest","Neural crest"),]$celltype <- "neural crest"
obj@meta.data[obj$celltype %in% c("Erythroblasts","Erythroblast","Ery"),]$celltype <- "erythroblast"
obj@meta.data[obj$celltype %in% c("Lateral plate Mesoderm 1","Lateral plate Mesoderm 2"),]$celltype <- "Lateral plate mesoderm"

# Rename cell types after inspection
obj@meta.data[obj$stage=="E5" & obj$celltype %in% c("B1.EPI"),]$celltype <- "EPI"
obj@meta.data[obj$stage=="E5" & obj$celltype %in% c("PE"),]$celltype <- "Hypoblast"
obj@meta.data[obj$stage=="E5" & obj$celltype %in% c("ICM-TE_trans","B1_B2"),]$celltype <- "Prelineage"
obj@meta.data[obj$stage=="E6" & obj$celltype %in% c("B1.EPI"),]$celltype <- "EPI"
obj@meta.data[obj$stage=="E6" & obj$celltype %in% c("EPI.PrE", "PE", "PrE"),]$celltype <- "Hypoblast"
obj@meta.data[obj$stage=="E7" & obj$celltype %in% c("PE","PrE"),]$celltype <- "Hypoblast"
obj@meta.data[obj$stage=="E8" & obj$celltype %in% c("PrE"),]$celltype <- "Hypoblast"
obj@meta.data[obj$stage=="E9" & obj$celltype %in% c("PrE"),]$celltype <- "YS endoderm"
obj@meta.data[obj$stage=="E10" & obj$celltype %in% c("PrE"),]$celltype <- "YS endoderm"
obj@meta.data[obj$stage=="E12" & obj$celltype %in% c("PrE"),]$celltype <- "YS endoderm"
obj@meta.data[obj$stage=="E12" & obj$celltype %in% c("PrE_derived"),]$celltype <- "ExE mesoderm"
obj@meta.data[obj$stage=="E14" & obj$celltype %in% c("PSA-EPI"),]$celltype <- "EPI"
obj@meta.data[obj$stage=="E14" & obj$celltype %in% c("PrE"),]$celltype <- "YS endoderm"
obj@meta.data[obj$stage=="E14" & obj$celltype %in% c("PrE_derived"),]$celltype <- "ExE mesoderm"
obj@meta.data[obj$stage=="E14" & obj$celltype %in% c("ExE_Mes"),]$celltype <- "ExE mesoderm"
obj@meta.data[obj$stage=="CS7" & obj$celltype %in% c("Amnion.Ecto"),]$celltype <- "Amnion"
obj@meta.data[obj$stage=="CS7" & obj$celltype %in% c("ExE_Mes"),]$celltype <- "ExE mesoderm"
obj@meta.data[obj$stage=="CS10" & obj$celltype %in% c("Neural tbue 1"),]$celltype <- "FB/MB"
obj@meta.data[obj$stage=="CS10" & obj$celltype %in% c("Neural tbue 2"),]$celltype <- "HB"
obj@meta.data[obj$stage=="CS10" & obj$celltype %in% c("Neural tbue 3"),]$celltype <- "SC"

# Use mixed level as the final layer
obj@meta.data[obj$stage %in% c("CS12","CS13-14","CS15-16") & obj$developmental.system %in% c("neuron","schwann","epidermis","sensory neuron","craniofacial","head mesoderm","somite","IM","limb","somatic LPM","splanchnic LPM","endothelium","blood","endoderm","PGC","epithelium","fibroblast"),]$celltype <- obj@meta.data[obj$stage %in% c("CS12","CS13-14","CS15-16") & obj$developmental.system %in% c("neuron","schwann","epidermis","sensory neuron","craniofacial","head mesoderm","somite","IM","limb","somatic LPM","splanchnic LPM","endothelium","blood","endoderm","PGC","epithelium","fibroblast"),]$developmental.system
obj@meta.data[obj$developmental.system=="neural progenitor" & grepl("rhombomere", obj$celltype),]$celltype <- "hindbrain"
obj@meta.data[obj$developmental.system=="neural progenitor" & (grepl("vesicle", obj$celltype) | grepl("olfactory", obj$celltype) | grepl("retinal", obj$celltype)),]$celltype <- "placode"
obj@meta.data[obj$developmental.system=="neural progenitor" & (grepl("plate", obj$celltype) | obj$celltype %in% c("pMNv-1","pA2","pMNv-2","pA1","pB2","pB1","pB3","pB4","pA3","dp6","dp5","pMN","p1","dp2","p0","dp3","dp4","p2","pMNs","p3","dp1")),]$celltype <- "spinal cord"
obj@meta.data[obj$developmental.system=="neural progenitor" & (grepl("midhindbrain", obj$celltype) | grepl("mesencephalon", obj$celltype)),]$celltype <- "midbrain"
obj <- obj[,!(obj$developmental.system=="neural progenitor" & grepl("GABAergic neuron precursor", obj$celltype))]
obj@meta.data[obj$developmental.system=="neural progenitor" & !(obj$celltype %in% c("hindbrain","midbrain","placode","spinal cord")),]$celltype <- "forebrain"

#saveRDS(obj, "~/data/embryo/combined_clean_v1.rds")
#obj <- readRDS("~/data/embryo/combined_clean_v1.rds")

# downsample to 6000 cells per stage
Idents(obj) <- obj$stage
obj <- subset(obj, downsample = 5000)
#saveRDS(obj, "~/data/embryo/combined_clean_v2_sub.rds")
#obj <- readRDS("~/data/embryo/combined_clean_v2_sub.rds")

obj$stage <- factor(obj$stage, levels = c("Zygote","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E12","E14","CS7","CS8_CytOrgin","CS9","CS10","CS11","CS12","CS13-14","CS15-16"))
obj$study <- as.character(obj$batch)
obj@meta.data[grepl("Yan et al", obj$study),]$study <- "Yan et al"
obj@meta.data[grepl("Meistermann et al", obj$study),]$study <- "Meistermann et al"
obj@meta.data[grepl("Petropoulos et al", obj$study),]$study <- "Petropoulos et al"
obj@meta.data[grepl("Yanagida et al", obj$study),]$study <- "Yanagida et al"
obj@meta.data[grepl("Xiang et al", obj$study),]$study <- "Xiang et al"
obj@meta.data[grepl("Tyser et al", obj$study),]$study <- "CS7"
obj@meta.data[grepl("CS8_CytOrgin_CS8_CytOrgin", obj$study),]$study <- "CS8_CytOrgin"
obj@meta.data[grepl("CS9", obj$study) | grepl("CS11", obj$study),]$study <- "ours"
obj@meta.data[grepl("CS10", obj$study),]$study <- "CS10"
obj@meta.data[grepl("CS12", obj$study) | grepl("CS12", obj$study) | grepl("CS13", obj$study) | grepl("CS15", obj$study),]$study <- "Xu et al"

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study)
obj <- NormalizeData(obj)
obj@assays$RNA$data.CS8_CytOrgin <- obj@assays$RNA$counts.CS8_CytOrgin
obj.list <- SplitObject(obj, split.by = "batch")
HVGs <- SelectIntegrationFeatures(obj.list, nfeatures = 4000)
HVGs <- c(HVGs, "TBXT")

obj <- ScaleData(obj, features = HVGs)
obj <- RunPCA(obj, features = HVGs)

obj <- IntegrateLayers(
  object = obj, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE,
  features = HVGs,
  k = 15
)
obj <- RunUMAP(obj, reduction = "integrated.mnn", dims = 1:50, reduction.name = "umap.mnn")
DimPlot(
  obj,
  reduction = "umap.mnn",
  group.by = c("stage"),
  combine = FALSE, label.size = 5, label = F,
  cols = colorRampPalette(brewer.pal(12,"Paired")[c(6,1:5,7:12)])(21), raster = F
)
saveRDS(obj, "~/data/embryo/combined_clean_sub_v2_extra_mnn.rds")

# Generate Circlized plot with tracks
obj <- JoinLayers(obj)
write.csv(obj@assays$RNA$data, "~/data/embryo/data_clean_sub_v2_mnn.csv")
write.csv(obj@meta.data, "~/data/embryo/meta_clean_sub_v2_mnn.csv")
write.csv(obj@reductions$umap@cell.embeddings, "~/data/embryo/umap_clean_sub_v2_mnn.csv")
write.csv(obj@reductions$integrated.mnn@cell.embeddings, "~/data/embryo/reduction_clean_sub_v2_mnn.csv")

stage_celltype_levels <- c('Zygote_Zygote','E1_L02C','E2_L04C','E3_L08C',
                           'E4_L08C','E4_Morula','E4_Prelineage',
                           'E5_Morula','E5_Prelineage','E5_TE','E5_ICM','E5_Hypoblast','E5_EPI',
                           'E6_TE','E6_CTB','E6_Hypoblast','E6_EPI',
                           'E7_TE','E7_CTB','E7_Hypoblast','E7_EPI',
                           'E8_TE','E8_CTB','E8_STB','E8_Hypoblast','E8_EPI',
                           'E9_CTB','E9_STB','E9_EPI','E9_YS endoderm',
                           'E10_CTB','E10_STB','E10_EPI','E10_YS endoderm',
                           'E12_CTB','E12_STB','E12_EVT','E12_EPI','E12_ExE mesoderm','E12_YS endoderm',
                           'E14_CTB','E14_STB','E14_EVT','E14_EPI','E14_ExE mesoderm','E14_YS endoderm',
                           'CS7_EPI','CS7_Primitive Streak','CS7_NasMes','CS7_EmMes','CS7_Advanced Mesoderm','CS7_HEP','CS7_erythroblast','CS7_Axial Mesoderm','CS7_endoderm','CS7_ExE mesoderm','CS7_Amnion',
                           
                           'CS8_CytOrgin_Epiblast','CS8_CytOrgin_PS','CS8_CytOrgin_Node','CS8_CytOrgin_Axial mesoderm','CS8_CytOrgin_Paraxial mesoderm','CS8_CytOrgin_LPM','CS8_CytOrgin_Endothelium','CS8_CytOrgin_HEP','CS8_CytOrgin_Erythroblast','CS8_CytOrgin_Endoderm','CS8_CytOrgin_YS endoderm','CS8_CytOrgin_YS mesoderm','CS8_CytOrgin_ExE mesoderm',
                           
                           'CS9_Surface ectoderm','CS9_Neural ectoderm','CS9_Notochord','CS9_PSM','CS9_Paraxial mesoderm','CS9_Heart progenitor','CS9_Mesenchyme','CS9_endothelium','CS9_HEP','CS9_erythroblast','CS9_endoderm','CS9_YS endoderm','CS9_YS mesoderm','CS9_ExE mesoderm','CS9_Amnion','CS9_Trophoblast',
                           
                           'CS10_epithelium','CS10_neural crest','CS10_FB/MB','CS10_HB','CS10_SC','CS10_NMP','CS10_Notochord','CS10_Intermediate Mesoderm','CS10_Paraxial mesoderm','CS10_Presomitic Mesoderm','CS10_Lateral plate mesoderm','CS10_Myocyte Progenitor','CS10_somite','CS10_Mesenchyme','CS10_Cardiomyocyte','CS10_endothelium','CS10_Erythroid','CS10_endoderm',
                           
                           'CS11_epidermis','CS11_neural crest','CS11_Neural progenitor-PAX6+','CS11_Roof plate','CS11_Floor plate','CS11_Forebrain/Midbrain','CS11_Hindbrain','CS11_spinal cord','CS11_Neural progenitor-proliferating','CS11_neuron','CS11_NMP','CS11_Notochord','CS11_PSM','CS11_somite','CS11_Dermomyotome','CS11_Myotome','CS11_Sclerotome','CS11_splanchnic LPM','CS11_Renal epithelium','CS11_Cardiomyocyte','CS11_Mesenchyme','CS11_fibroblast','CS11_Limb progenitor','CS11_Endothelium/HEP','CS11_erythroblast','CS11_Macrophage','CS11_head mesoderm','CS11_Cranifacial mesoderm','CS11_Intestine','CS11_Hepatocyte',
                           
                           'CS12_epidermis','CS12_epithelium','CS12_forebrain','CS12_midbrain','CS12_hindbrain','CS12_spinal cord','CS12_placode','CS12_neuron','CS12_sensory neuron','CS12_schwann','CS12_craniofacial','CS12_head mesoderm','CS12_IM','CS12_splanchnic LPM','CS12_somite','CS12_endothelium','CS12_blood','CS12_somatic LPM','CS12_limb','CS12_fibroblast','CS12_endoderm',
                               
                           'CS13-14_epidermis','CS13-14_epithelium','CS13-14_forebrain','CS13-14_hindbrain','CS13-14_midbrain','CS13-14_placode','CS13-14_spinal cord','CS13-14_neuron','CS13-14_sensory neuron','CS13-14_schwann','CS13-14_craniofacial','CS13-14_head mesoderm','CS13-14_IM','CS13-14_splanchnic LPM','CS13-14_somite','CS13-14_endothelium','CS13-14_blood','CS13-14_somatic LPM','CS13-14_limb','CS13-14_fibroblast','CS13-14_endoderm','CS13-14_PGC',
                           
                           'CS15-16_epidermis','CS15-16_epithelium','CS15-16_forebrain','CS15-16_hindbrain','CS15-16_midbrain','CS15-16_placode','CS15-16_spinal cord','CS15-16_neuron','CS15-16_sensory neuron','CS15-16_schwann','CS15-16_craniofacial','CS15-16_head mesoderm','CS15-16_IM','CS15-16_splanchnic LPM','CS15-16_somite','CS15-16_endothelium','CS15-16_blood','CS15-16_somatic LPM','CS15-16_limb','CS15-16_fibroblast','CS15-16_endoderm')
saveRDS(obj1, "combined_clean_sub_v2_extra_mnn_revised_cs8_revised.rds")
obj$stage_celltype <- factor(paste(obj$stage, obj$celltype, sep = "_"))
obj$stage_celltype <- factor(paste(obj$stage, obj$celltype, sep = "_"), level = stage_celltype_levels)


TF_data <- as.data.frame(t(as.matrix(AverageExpression(obj, assays = "RNA", features = c('BMP2', 'BMP4', 'BMP7', 'NODAL', 'WNT3', 'WNT6', 'FGF2', 'FGF8', 'SHH', 'ALDH1A2'), group.by = "stage_celltype", layer = "scale.data")$RNA)))
TF_data <- as.data.frame(t(as.matrix(AverageExpression(obj, assays = "RNA", features = c('NANOG', 'SOX2', 'OTX2', 'PAX6', 'TBXT', 'EOMES', 'SOX17', 'FOXA2', 'GATA3', 'GATA6'), group.by = "stage_celltype", layer = "scale.data")$RNA)))
for (gene in colnames(TF_data)) {
  TF_data[[gene]] <- (TF_data[[gene]] - min(TF_data[[gene]]))/(max(TF_data[[gene]]) - min(TF_data[[gene]]))
}
TF_data$group <- unlist(lapply(strsplit(stage_celltype_levels, "_"), FUN = function(x) {x[1]}))
TF_data$celltype <- unlist(lapply(strsplit(stage_celltype_levels, "_"), FUN = function(x) {x[2]}))
write.csv(TF_data, "~/data/embryo/TF_data_avg.csv")

ratio_data <- prop.table(table(obj$batch, obj$stage_celltype), margin = 2)
ratio_data <- dcast(as.data.frame(ratio_data), Var1~Var2)
rownames(ratio_data) <- ratio_data$Var1
ratio_data <- ratio_data[,-1]
ratio_data <- as.data.frame(t(ratio_data))
ratio_data$group <- unlist(lapply(strsplit(stage_celltype_levels, "_"), FUN = function(x) {x[1]}))
ratio_data$celltype <- unlist(lapply(strsplit(stage_celltype_levels, "_"), FUN = function(x) {x[2]}))
write.csv(ratio_data, "~/data/embryo/ratio_data.csv", row.names = T, col.names = T)


stage_ct_ratio_data <- prop.table(table(obj$stage, obj$stage_celltype), margin = 1)
stage_ct_ratio_data <- dcast(as.data.frame(stage_ct_ratio_data), Var1~Var2)
rownames(stage_ct_ratio_data) <- stage_ct_ratio_data$Var1
stage_ct_ratio_data <- stage_ct_ratio_data[,-1]
stage_ct_ratio_data <- as.data.frame(t(stage_ct_ratio_data))
for (i in colnames(stage_ct_ratio_data)) {
  stage_ct_ratio_data[[i]] <- stage_ct_ratio_data[[i]] / max(stage_ct_ratio_data[[i]])
}
stage_ct_ratio_data$group <- unlist(lapply(strsplit(stage_celltype_levels, "_"), FUN = function(x) {x[1]}))
stage_ct_ratio_data$celltype <- unlist(lapply(strsplit(stage_celltype_levels, "_"), FUN = function(x) {x[2]}))
write.csv(stage_ct_ratio_data, "~/data/embryo/stage_ct_ratio_data.csv", row.names = T, col.names = T)

