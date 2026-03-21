# 2026-03-21
# Setting up ----

## Make directories
analysis_step <- "020_all"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(Seurat)
source("Rscripts/helpers.R")

## Helper functions
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Clustering w/o integration ----
seu <- readRDS("RDSfiles/seu_011_singlet.RDS")

seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
hvg <- VariableFeatures(seu)
cc.genes <- Seurat::cc.genes
seu <- CellCycleScoring(seu, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = FALSE)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:npcs)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 4))
ggsave("cluster.pdf", path = plot_path, width = 50, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# Check markers ----
features <- readLines("aux_data/gene_set/010_cellgroup.txt")
sapply(features, save_fp, seu, fp_path)

markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

VlnPlot(seu, features = c("nFeature_RNA"), pt.size = 0) & 
  theme_classic(base_size = 6) & NoLegend() &
  labs(x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("QC_vln_nFeature.pdf", path = plot_path, width = 45, height = 35, units = "mm")

VlnPlot(seu, features = c("nCount_RNA"), pt.size = 0) & 
  theme_classic(base_size = 6) & NoLegend() &
  labs(x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("QC_vln_nCount.pdf", path = plot_path, width = 45, height = 35, units = "mm")

VlnPlot(seu, features = c("percent.mt"), pt.size = 0) & 
  theme_classic(base_size = 6) & NoLegend() &
  labs(x = NULL, y = NULL) &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("QC_vln_mt.pdf", path = plot_path, width = 45, height = 35, units = "mm")

# Add celltype1 annotation ----
seu <- subset(seu, subset = seurat_clusters %in% c(30), invert = TRUE)
seu$celltype1 <- "Epi"
seu$celltype1[seu$seurat_clusters %in% c(0,1,2,7,8,18,29)] <- "T_NK"
seu$celltype1[seu$seurat_clusters %in% c(6,23,26)] <- "Bcell"
seu$celltype1[seu$seurat_clusters %in% c(3,9,13,14)] <- "Plasma"
seu$celltype1[seu$seurat_clusters %in% c(4,20,25,31)] <- "Mye"
seu$celltype1[seu$seurat_clusters %in% c(19)] <- "Mast"
seu$celltype1[seu$seurat_clusters %in% c(5,15,17,21,27)] <- "Fibro"
seu$celltype1[seu$seurat_clusters %in% c(16,22)] <- "Myo"
seu$celltype1[seu$seurat_clusters %in% c(10)] <- "EC"
seu$celltype1 <- factor(seu$celltype1, levels = c("Epi", "T_NK", "Bcell", "Plasma", "Mye", "Mast", "Fibro", "Myo", "EC"))
DimPlot(seu, group.by = "celltype1", cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 2) &
  theme_panel() & NoAxes()
ggsave("celltype1.pdf", path = plot_path, width = 35, height = 30, units = "mm")

# save RDS ----
# seu <- DietSeurat(seu)
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
seu$RNA_snn_res.0.5 <- NULL
seu$RNA_snn_res.1 <- NULL
saveRDS(seu, file = "RDSfiles/seu_020_all.RDS")

# EDA ----
seu <- readRDS("RDSfiles/seu_020_all.RDS")
seu <- NormalizeData(seu)
add_feat <- c("BMPR1A", "GREM1", "WNT5A", "RSPO3", "LRG1", "CD38", "ENG")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- "LRG1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), raster = TRUE, 
            max.cutoff = "q25", min.cutoff = "q0", order = TRUE, pt.size = 2)  &
  theme_panel() & NoAxes()
ggsave(paste0(add_feat, "q25.pdf"), path = fp_path, width = 35, height = 30, units = "mm")

features <- c("LRG1", "GREM1", "BMPR1A")
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
