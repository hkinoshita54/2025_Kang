analysis_step <- "011_all_cells"

# Make directories ----
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

# load packages ----
library(tidyverse)
library(Seurat)
library(decontX)

# load data ----
seu <- readRDS("RDSfiles/seu_012_cellgroup.RDS")
DimPlot(seu, group.by = "cellgroup") & NoAxes()
View(seu[[]])

# refine annotation ----
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
DimPlot(seu, cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes() &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster_res1.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "tissue_type") & NoAxes()

## find markers for annotation
seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

## remove doublet-like or highly contaminated cells
seu <- subset(seu, subset = seurat_clusters %in% c(16,17,24,27,28), invert = TRUE)

# celltype1 annotation ----
seu$celltype1 <- ""
seu$celltype1[seu$seurat_clusters %in% c(20)] <- "Epi-Int"
seu$celltype1[seu$seurat_clusters %in% c(13,14)] <- "Epi-St"
seu$celltype1[seu$seurat_clusters %in% c(7,23,25)] <- "B"
seu$celltype1[seu$seurat_clusters %in% c(2,5,19)] <- "Plasma"
seu$celltype1[seu$seurat_clusters %in% c(0,1,3,6,10,21)] <- "T/NK"
seu$celltype1[seu$seurat_clusters %in% c(8,11,26,29,31)] <- "Myeloid"
seu$celltype1[seu$seurat_clusters %in% c(18)] <- "Mast"
seu$celltype1[seu$seurat_clusters %in% c(4,9,22)] <- "Fibro."
seu$celltype1[seu$seurat_clusters %in% c(15,30,32)] <- "EC"
seu$celltype1[seu$seurat_clusters %in% c(12)] <- "Peri."
seu$celltype1 <- factor(
  seu$celltype1,
  levels = c("Epi-St", "Epi-Int", "B", "Plasma", "T/NK", "Myeloid", "Mast", "Fibro.", "EC", "Peri.")
)
DimPlot(seu, group.by = "celltype1", cols = "polychrome", label = TRUE, repel = TRUE) & NoAxes()
ggsave("celltype1.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

## find markers for each annotated cluster
Idents(seu) <- "celltype1"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_anno.csv"))

# decontx for visualization ----
sce <- as.SingleCellExperiment(seu, assay = "RNA")
sce <- decontX(
  sce,
  z = seu$celltype1,      # or seu$seurat_clusters
  estimateDelta = TRUE,
  verbose = TRUE
)
decont_counts <- decontXcounts(sce)
seu[["decontX"]] <- CreateAssayObject(counts = decont_counts)
DefaultAssay(seu) <- "decontX"
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)

## dotplot
features <- c(
  "CLDN18","MUC5AC",
  "ALDOB","APOA4",
  "MS4A1","CD79A",
  "MZB1","JCHAIN",
  "TRAC","NKG7",
  "TYROBP","FCER1G",
  "TPSAB1","KIT",
  "DCN","COL1A1",
  "VWF","PLVAP",
  "RGS5","ACTA2"
)
DotPlot(seu, group.by = "celltype1", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 7, height = 4, units = "in", dpi = 150)

# additional plots
add_feat <- "LRG1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"),
            max.cutoff = "q25", min.cutoff = "q0", order = TRUE, pt.size = 2,
) + NoAxes() 
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3.6, height = 3, units = "in", dpi = 150)

DotPlot(seu, features = add_feat, group.by = "celltype1")
ggsave(paste0(add_feat, "_dotplot.png"), path = plot_path, width = 3.3, height = 4, units = "in", dpi = 150)

add_feat <- "CD38"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"),
            max.cutoff = "q25", min.cutoff = "q0", order = TRUE, pt.size = 2,
) + NoAxes() 
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3.6, height = 3, units = "in", dpi = 150)

DotPlot(seu, features = add_feat, group.by = "celltype1")
ggsave(paste0(add_feat, "_dotplot.png"), path = plot_path, width = 3.3, height = 4, units = "in", dpi = 150)
