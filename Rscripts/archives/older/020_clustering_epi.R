# 
analysis_step <- "020_clustering_epi"

# load packages ----
library(tidyverse)
library(ggrepel)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path_root <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Helper function ----
cluster_harm = function(seu_obj, npcs, res){
  seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
  seu <- RunPCA(seu, npcs = npcs)
  seu <- IntegrateLayers(
    object = seu, method = HarmonyIntegration,
    orig.reduction = "pca", 
    new.reduction = "harmony")
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
  return(seu)
}

# Load data ----
seu_all <- readRDS("RDSfiles/seu_012_cellgroup.RDS")
load("RDSfiles/cellgroup_names.Rdata")

# ITER1: cluster w/ harmony ----
plot_path <- file.path(plot_path_root, "iter1")
fs::dir_create(c(plot_path))
seu <- seu_all[, epi_names]
seu <- JoinLayers(seu)
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$pt_id) 
seu <- cluster_harm(seu, npcs = 30, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "tissue_type") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("type.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)
seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE)

# ITER2: remove probable heterotypic multiplets, then re-cluster ----
plot_path <- file.path(plot_path_root, "iter2")
fs::dir_create(c(plot_path))
seu <- subset(seu, idents = c(17,18), invert = TRUE)
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$pt_id) 
seu <- cluster_harm(seu, npcs = 30, res = 1)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 2))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "orig.ident") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "Lauren.Classification") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))
ggsave("type.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)

#### NOT well clustered in terms of markers and tissue type (T or NT)
#### probably related to original data (raw count estimated from transformed data)
#### the following are tentative annotation

# Add cellgroup annotation ----
prog_names <- colnames(seu)[seu$seurat_clusters %in% c(17,9)]
pit_names <- colnames(seu)[seu$seurat_clusters %in% c(12,1,15,11)]
neck_names <- colnames(seu)[seu$seurat_clusters %in% c(5,8,18)]
chief_names <- colnames(seu)[seu$seurat_clusters %in% c(10)]
pit2_names <- colnames(seu)[seu$seurat_clusters %in% c(0,2,16)]
neck2_names <- colnames(seu)[seu$seurat_clusters %in% c(4,3,14)]
pariet_names <- colnames(seu)[seu$seurat_clusters %in% c(13)]
si_names <- colnames(seu)[seu$seurat_clusters %in% c(6,7)]
save(prog_names, pit_names, neck_names, chief_names, pit2_names, neck2_names, pariet_names, si_names, 
     file = file.path("RDSfiles", "epi_type_names.Rdata"))

seu$celltype <- ""
seu$celltype[prog_names] <- "Prog."
seu$celltype[pit_names] <- "Pit"
seu$celltype[pit2_names] <- "Pit2"
seu$celltype[neck_names] <- "Neck"
seu$celltype[neck2_names] <- "Neck2"
seu$celltype[chief_names] <- "Chief"
seu$celltype[pariet_names] <- "Pariet."
seu$celltype[si_names] <- "SI"
seu$celltype <- factor(seu$celltype, levels = c("Prog.", "Pit", "Pit2", "Neck", "Neck2", "Chief", "Pariet.", "SI"))
DimPlot(seu, group.by = "celltype", cols = "polychrome") & NoAxes()
ggsave("celltype.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# save RDS
saveRDS(seu, file = file.path("RDSfiles", "seu_020_epi.RDS"))

# additional plots ----
add_feat <- "rna_MUC6"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

add_feat <- "LRG1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            max.cutoff = "q80", min.cutoff = "q0"
) + NoAxes() 
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3.6, height = 3, units = "in", dpi = 150)

DotPlot(seu, features = "LRG1", group.by = "tissue_type")
ggsave("LRG1_dotplot.png", path = plot_path, width = 4, height = 4, units = "in", dpi = 150)

seu_cp <- seu
seu <- subset(seu, subset = tissue_type %in% c("T"))
DotPlot(seu, features = "LRG1", group.by = "Lauren.Classification")
ggsave("LRG1_dotplot_T_Laurens.png", path = plot_path, width = 5, height = 4, units = "in", dpi = 150)
DotPlot(seu, features = "LRG1", group.by = "TCGA.Molecular.Grouping")
ggsave("LRG1_dotplot_T_TCGA.png", path = plot_path, width = 4, height = 4, units = "in", dpi = 150)
