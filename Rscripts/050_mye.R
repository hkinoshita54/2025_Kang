# 2026-03-21
# Setting up ----

## Make directories
analysis_step <- "050_mye"
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(openxlsx2)
library(Seurat)
source("Rscripts/helpers.R")

## Helper functions
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 4) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Clustering w/o integration ----
seu <- readRDS("RDSfiles/seu_020_all.RDS")
seu <- subset(seu, subset = celltype1 == "Mye")

# w/o integration
seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
hvg <- VariableFeatures(seu)
seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
npcs <- 50
seu <- RunPCA(seu, npcs = npcs)
seu <- FindNeighbors(seu, dims = 1:npcs)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:npcs)

# w/ harmony
# seu[["RNA"]] <- split(seu[["RNA"]], f = seu$Platform)
# seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 4000)
# hvg <- VariableFeatures(seu)
# seu <- ScaleData(seu, features = hvg, vars.to.regress = c("S.Score", "G2M.Score"))
# 
# npcs <- 50
# seu <- RunPCA(seu, npcs = npcs)
# seu <- IntegrateLayers(
#   object = seu, method = HarmonyIntegration,
#   orig.reduction = "pca", 
#   new.reduction = "harmony")
# seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
# seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
# seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.pdf", path = plot_path, width = 40, height = 30, units = "mm")

# Check markers ----
add_feat <- c("S100A9", "APOE", "MRC1", "IL1B", "CLEC4A")
sapply(add_feat, save_fp, seu, fp_path)

add_feat <- c("EPCAM", "PGA4", "JCHAIN", "CD79A", "CD3D", "PDGFRA", "DCN")
sapply(add_feat, save_fp, seu, fp_path)

# seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

# remove contamination
seu <- subset(seu, subset = seurat_clusters %in% c(2,9,10,15:18), invert = TRUE) # remove obvious contamination

DimPlot(seu, cols = "polychrome", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "seurat_clusters") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster2.pdf", path = plot_path, width = 40, height = 30, units = "mm")

DimPlot(seu, group.by = "patient", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") 
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# save RDS ----
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
seu$RNA_snn_res.0.5 <- NULL
seu$RNA_snn_res.1 <- NULL
saveRDS(seu, file = "RDSfiles/seu_050_mye.RDS")
