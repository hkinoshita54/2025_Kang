# 2026-03-10
# Make directories ----
analysis_step <- "020_all"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## Helper function
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

## load packages
library(tidyverse)
library(Seurat)
source("Rscripts/theme_panel.R")

# Clustering w/ Harmony integration ----
seu <- readRDS("RDSfiles/seu_011_singlet.RDS")
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$pt_id)
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
hvg <- VariableFeatures(seu)

cc.genes <- Seurat::cc.genes
seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"), block.size = 4000)

npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.pdf", path = plot_path, width = 45, height = 30, units = "mm")

DimPlot(seu, group.by = "pt_id", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue.pdf", path = plot_path, width = 30, height = 30, units = "mm")

## check cell cycle regression
sapply(c("S.Score", "G2M.Score"), save_fp, seu, fp_path)
DimPlot(seu, group.by = "Phase", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Phase")
ggsave("phase.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# # Check markers by feature plots
# features <- readLines("aux_data/gene_set/global_markers.txt")
# sapply(features, save_fp, seu, fp_path)

# check markers ----
seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

add_feat = "percent.mt"
VlnPlot(seu, features = add_feat, pt.size = 0) +
  theme_panel() & NoLegend() &
  labs(title = add_feat, x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 50, height = 35, units = "mm")

add_feat = "nFeature_RNA"
VlnPlot(seu, features = add_feat, pt.size = 0) +
  theme_panel() & NoLegend() &
  labs(title = add_feat, x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 50, height = 35, units = "mm")

# Add celltype1 annotation ----
seu_all <- seu
seu <- subset(seu, subset = seurat_clusters %in% c(12,14,21,22), invert = TRUE)
seu$celltype1 <- "Epi"
seu$celltype1[seu$seurat_clusters %in% c(0,2,7,9,17)] <- "T_NK"
seu$celltype1[seu$seurat_clusters %in% c(5,20)] <- "Bcell"
seu$celltype1[seu$seurat_clusters %in% c(1,16)] <- "Plasma"
seu$celltype1[seu$seurat_clusters %in% c(6,10,19)] <- "Mye"
seu$celltype1[seu$seurat_clusters %in% c(15)] <- "Mast"
seu$celltype1[seu$seurat_clusters %in% c(3,8)] <- "Fibro"
seu$celltype1[seu$seurat_clusters %in% c(11)] <- "Myo"
seu$celltype1[seu$seurat_clusters %in% c(13)] <- "EC"
seu$celltype1 <- factor(seu$celltype1, levels = c("Epi", "T_NK", "Bcell", "Plasma", "Mye", "Mast", "Fibro", "EC", "Myo"))
DimPlot(seu, group.by = "celltype1", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes()
ggsave("celltype1.pdf", path = plot_path, width = 35, height = 31, units = "mm")

DimPlot(seu, group.by = "pt_id", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt2.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue2.pdf", path = plot_path, width = 30, height = 30, units = "mm")

Idents(seu) <- "celltype1"
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers_annotated.csv"))

# Dot plot
Idents(seu) <- "celltype1"
features <- c(
  "EPCAM", 
  "CD3E", 
  "MS4A1", 
  "JCHAIN",
  "TYROBP", 
  "TPSB2",
  "DCN", 
  "PECAM1", 
  "RGS5"
)
DotPlot(seu, group.by = "celltype1", features = features, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot.pdf", path = plot_path, width = 60, height = 40, units = "mm")

# save RDS ----
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
saveRDS(seu, file = "RDSfiles/seu_020_all.RDS")

# EDA ----
seu <- readRDS("RDSfiles/seu_020_all.RDS")
seu <- NormalizeData(seu)
add_feat <- c("BMPR1A", "GREM1", "WNT5A", "RSPO3", "LRG1", "CD38", "ENG")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("GREM1")
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), raster = TRUE, 
            max.cutoff = "q95", min.cutoff = "q0", order = TRUE, pt.size = 2)  &
  theme_panel() & NoAxes()
ggsave(paste0(add_feat, ".pdf"), path = fp_path, width = 35, height = 30, units = "mm")

features <- c("LRG1", "GREM1", "BMPR1A")
DotPlot(seu, group.by = "celltype1", features = features, dot.scale = 2.3) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      title = "% Expr",
      title.position = "top"
    ),
    colour = guide_colorbar(
      title = "Avg Expr",
      title.position = "top",
      barheight = grid::unit(10, "mm"),
      barwidth  = grid::unit(3, "mm")
    )
  )+
  theme(
    axis.text.x = element_text(margin = margin(t = -3, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot.pdf", path = plot_path, width = 30, height = 36, units = "mm")

features <- c("ENG")
DotPlot(seu, group.by = "celltype1", features = features, dot.scale = 2.3) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  guides(
    size = guide_legend(
      title = "% Expr",
      title.position = "top"
    ),
    colour = guide_colorbar(
      title = "Avg Expr",
      title.position = "top",
      barheight = grid::unit(10, "mm"),
      barwidth  = grid::unit(3, "mm")
    )
  )+
  theme(
    axis.text.x = element_text(margin = margin(t = -1, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot_eng.pdf", path = plot_path, width = 25, height = 36, units = "mm")
