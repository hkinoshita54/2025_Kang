# Continued from 020_
# epithelial cells
# 2026-03-10

# Make directories ----
analysis_step <- "030_epi"
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(plot_path, res_path, fp_path))

## load packages
library(tidyverse)
library(openxlsx2)
library(Seurat)
library(msigdbr)
library(UCell)
library(ggpubr)
library(pals)
source("Rscripts/theme_panel.R")

## Helper functions
save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred"), raster = TRUE, pt.size = 2) &
      theme_panel() & NoAxes() & NoLegend()
    ggsave(paste0(feature, ".pdf"), path = path, width = 25, height = 30, units = "mm")
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}

# Clustering ----
seu <- readRDS("RDSfiles/seu_020_all.RDS")

seu <- DietSeurat(seu)
seu <- subset(seu, subset = celltype1 == "Epi")

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
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)

DimPlot(seu, cols = "polychrome", raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Cluster") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.pdf", path = plot_path, width = 45, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "pt_id", raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & NoLegend() & labs(title = "Patient") &
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("pt.pdf", path = plot_path, width = 25, height = 30, units = "mm")

seu <- JoinLayers(seu)
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 20) %>%
  ungroup() -> top20
write_csv(top20, file = file.path(res_path, "top20_markers.csv"))

add_feat <- "EPCAM"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), raster = TRUE, 
            max.cutoff = "q95", min.cutoff = "q0", order = TRUE, pt.size = 2)  &
  theme_panel() & NoAxes()
ggsave(paste0(add_feat, ".pdf"), path = fp_path, width = 35, height = 30, units = "mm")


## annotation
seu <- subset(seu, subset = seurat_clusters %in% c(17,18), invert = TRUE)
seu$celltype2 <- "Tumor"
seu$celltype2[seu$seurat_clusters %in% c(15)] <- "EEC"
seu$celltype2[seu$seurat_clusters %in% c(13)] <- "Pariet"
seu$celltype2[seu$seurat_clusters %in% c(7,8)] <- "IM"
seu$celltype2[seu$celltype2 == "Tumor" & seu$tissue_type == "N"] <- "Non-T"
seu$celltype2 <- factor(seu$celltype2, levels = c("Non-T", "EEC", "Pariet", "IM", "Tumor"))
DimPlot(seu, group.by = "celltype2", cols = "alphabet", raster = TRUE, raster.dpi = c(600, 600), pt.size = 4) &
  theme_panel() & NoAxes()
ggsave("celltype2.pdf", path = plot_path, width = 30, height = 30, units = "mm")

DimPlot(seu, group.by = "tissue_type", raster = TRUE, pt.size = 4) &
  theme_panel() & NoAxes() & labs(title = "Tissue")
ggsave("tissue2.pdf", path = plot_path, width = 30, height = 30, units = "mm")

# save RDS ----
seu[["RNA"]]$scale.data <- NULL
saveRDS(seu, file = "RDSfiles/seu_030_epi.RDS")
seu <- readRDS("RDSfiles/seu_030_epi.RDS")

# UCell scoring ----
seu <- subset(seu, subset = celltype2 == "Tumor")
H <- msigdbr(species    = "Homo sapiens", collection = "H")
H$gs_name <- gsub("HALLMARK_", "", H$gs_name)
H <- split(H$gene_symbol, H$gs_name)

yap_list <- readRDS("RDSfiles/yap_list_human.RDS")
yap_list <- yap_list[4:5]

seu <- AddModuleScore_UCell(
  seu,
  features = H,
  chunk.size = 1000,
  ncores = 8,
)

seu <- AddModuleScore_UCell(
  seu,
  features = yap_list,
  chunk.size = 1000,
  ncores = 8,
)

## make ucell scores into z-scores
### change names for brevity
colnames(seu@meta.data)[colnames(seu@meta.data) == "EPITHELIAL_MESENCHYMAL_TRANSITION_UCell"] <- "EMT_UCell"
colnames(seu@meta.data)[colnames(seu@meta.data) == "INFLAMMATORY_RESPONSE_UCell"] <- "INFL_RESP_UCell"

ucell_cols <- grep("_UCell$", colnames(seu[[]]), value = TRUE)

for (cc in ucell_cols) {
  seu[[cc]] <- as.numeric(scale(seu@meta.data[[cc]]))
}

## create a composite score
seu[["TME_score"]] <- rowMeans(
  seu@meta.data[, c(
    "EMT_UCell",
    "HYPOXIA_UCell",
    "INFL_RESP_UCell"
  )],
  na.rm = TRUE
)

# FeatureScatter with YAP ----
## TME_score
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_TME.pdf", path = plot_path, width = 35, height = 30, units = "mm")

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_TME_no_p.pdf", path = plot_path, width = 35, height = 30, units = "mm", device = cairo_pdf)

## EMT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_EMT.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "EMT_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_EMT_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

## HYPOXIA
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_HYPOXIA.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "HYPOXIA_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_HYPOXIA_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

## INFL_RESP
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_YAP_INFL_RESP.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "YAP_TARGETS_UCell", "INFL_RESP_UCell", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_YAP_INFL_RESP_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

# FeatureScatter with G2M ----
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  stat_cor(method = "spearman") &
  theme_panel() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90)) & NoLegend()
ggsave("fscatter_G2M_TME.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

### remove stats >> add them later in PPT
FeatureScatter(seu, "G2M_CHECKPOINT_UCell", "TME_score", group.by = "celltype2", cols = "#565656", pt.size = 0.1) & 
  theme_panel() & NoLegend() &
  theme(plot.title = element_blank(), axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave("fscatter_G2M_TME_no_p.pdf", path = plot_path, width = 40, height = 30, units = "mm", device = cairo_pdf)

# EDA ----
seu <- readRDS("RDSfiles/seu_030_epi.RDS")
seu <- subset(seu, subset = celltype2 == "Tumor")
add_feat <- c("SPP1", "LRG1", "CD38")
sapply(add_feat, save_fp, seu, fp_path)

DotPlot(seu, group.by = "pt_id", features = add_feat, dot.scale = 2.5) + 
  theme_panel() + RotatedAxis() + labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(margin = margin(t = -2, unit = "mm")),
    plot.margin = margin(t = 2, r = 2, b = 3.5, l = 2, unit = "mm"),
    legend.title = element_text(size = 6),
    legend.text  = element_text(size = 5.5),
    legend.key.size = grid::unit(2.5, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm")
  ) 
ggsave("dotplot_SPP1_LRG1_CD38.pdf", path = plot_path, width = 45, height = 60, units = "mm")

# Intestinal vs Diffuse ----
seu2 <- subset(seu, subset = Lauren.Classification %in% c("Diffuse", "Intestinal"))
seu2$Laurens <- factor(seu2$Lauren.Classification, levels = c("Intestinal", "Diffuse"))
add_feat = "YAP_TARGETS_UCell"
VlnPlot(seu2, group.by = "Laurens", features = add_feat, pt.size = 0) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black") &
  theme_panel() & NoLegend() &
  labs(title = "YAP_TARGETS", x = NULL, y = "z-score") &
  theme(axis.title.y = element_text(angle = 90), axis.text.x = element_text(angle = 90))
ggsave(paste0("vln_", add_feat, ".pdf"), path = plot_path, width = 25, height = 35, units = "mm")
