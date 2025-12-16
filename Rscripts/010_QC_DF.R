analysis_step <- "010_QC"

# load packages ----
library(tidyverse)
library(data.table)
library(Matrix)
library(readxl)
library(Seurat)
library(DoubletFinder)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path))

# Create a Seurat object with some patient data ----
pts_data <- read_excel("data/13059_2022_2828_MOESM2_ESM.xlsx", skip = 1) %>% as.data.frame
pts_data$`Donor ID` <- sub(pattern = "D", replacement = "", pts_data$`Donor ID`)
pts_data <- pts_data[,-c(1,3,4,21:24,26)] # remove unecessary columns

mtx <- readRDS("data/GSE206785_raw_matrix.rds") # this matrix is created in 000_checking_data_format.R
seu <- CreateSeuratObject(counts = t(mtx), min.cells = 3, min.features = 200)
seu$orig.ident <- colnames(seu) %>% strsplit("-") %>% lapply("[",3) %>% unlist()
seu$pt_id <- gsub(pattern = "[A-Z]", replacement = "", seu$orig.ident)
seu$tissue_type <- gsub("[0-9]", "", seu$orig.ident)
seu[[]] <- left_join(seu[[]], pts_data, by = c("pt_id" = "Donor ID"))
seu[[]] <- mutate_if(seu[[]], is.character, as.factor)

saveRDS(seu, "RDSfiles/seu_010_raw.RDS")

# initial QC ----
Idents(seu) <- "orig.ident"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_unfil.png", path = plot_path, width = 25, height = 3, units = "in", dpi = 150)
### looks like it is already filtered based on percent.mt
# mito_cutoff <- quantile(seu$percent.mt, 0.90) # set threshold at 90 percentile
# seu <- subset(seu, subset = percent.mt < mito_cutoff)
seu <- subset(seu, subset = orig.ident == "180305N", invert = TRUE) # remove sample with only 2 cells
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_wo_small.png", path = plot_path, width = 25, height = 3, units = "in", dpi = 150)

# Create Seurat object list for DoubleFinder ----
seu_list <- SplitObject(seu, split.by = "pt_id")

# DoubletFinder ----
sample_name <- seu$pt_id %>% unique()
source(file.path("Rscripts", "011_DoubletFinder.R"))

## merge the list ----
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)])
# rm(seu_list)
names(seu@meta.data)
seu@meta.data <- seu@meta.data[,c(1:23,25,27)]
seu@meta.data <- mutate_if(seu@meta.data, is.character, as.factor)


# Filter doublet ----

seu <- subset(seu, subset = doublet_finder == "Singlet")
Idents(seu) <- "orig.ident"
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_singlet.png", path = plot_path, width = 25, height = 3, units = "in", dpi = 150)

# Save filtered object
seu <- JoinLayers(seu)
seu[["RNA"]]$scale.data <- NULL
seu[["RNA"]]$data <- NULL
saveRDS(seu, file.path("RDSfiles", "seu_011_singlet.RDS"))
seu_raw<-seu    # copy seurat object w/ raw count alone

# Clustering w/ Harmony integration ----
plot_path <- file.path("plot", analysis_step, "harmony")    # set a new plot directory
fs::dir_create(plot_path)
seu[["RNA"]] <- split(seu[["RNA"]], f = seu$pt_id)
seu <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData()
npcs <- 30    # This is arbitrary. Just pick one from e.g. 15, 20, 30 or 50.
seu <- RunPCA(seu, npcs = npcs)
seu <- IntegrateLayers(
  object = seu, method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "harmony")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:npcs, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)
DimPlot(seu, group.by = "pt_id") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("sample.png", path = plot_path, width = 5.5, height = 3, units = "in", dpi = 150)

# Check markers by feature plots
fp_path <- file.path(plot_path, "feature_plot")
fs::dir_create(c(fp_path))
features <- readLines("aux_data/gene_set/010_cellgroup.txt")

save_fp <- function(feature, seu, path){
  tryCatch({
    p <- FeaturePlot(seu, features = feature, cols = c("lightgrey","darkred")) +
      NoAxes() + NoLegend()
    ggsave(paste0(feature, ".png"), plot = p, path = path, 
           width = 3, height = 3, units = "in", dpi = 150)
  }, error = function(e){cat("ERROR :", conditionMessage(e), "\n")})
}
sapply(features, save_fp, seu, fp_path)

# adjust resolution
seu <- FindClusters(seu, resolution = 2, verbose = FALSE)
DimPlot(seu, label = TRUE, repel = TRUE, cols = "polychrome") + NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 3))
ggsave("cluster_res2.png", path = plot_path, width = 5, height = 3, units = "in", dpi = 150)

# Add cellgroup annotation ----
epi_names <- colnames(seu)[seu$seurat_clusters %in% c(14,17,29,32)]
str_names <- colnames(seu)[seu$seurat_clusters %in% c(4,9,13,18,27,30,33,36,37,41,47,48,49)]
b_names <- colnames(seu)[seu$seurat_clusters %in% c(3,7,12,15,19,22,26,35,38,39,42,43,44)]
t_names <- colnames(seu)[seu$seurat_clusters %in% c(0,1,2,5,6,8,10,16,23,25,31)]
mye_names <- colnames(seu)[seu$seurat_clusters %in% c(11,20,21,24,28,34,40,45,46)]
save(epi_names, str_names, b_names, t_names, mye_names, file = file.path("RDSfiles", "cellgroup_names.Rdata"))

seu$cellgroup <- ""
seu$cellgroup[epi_names] <- "Epi."
seu$cellgroup[str_names] <- "Str."
seu$cellgroup[b_names] <- "Bcell"
seu$cellgroup[t_names] <- "Tcell"
seu$cellgroup[mye_names] <- "Mye."
seu$cellgroup <- factor(seu$cellgroup, levels = c("Epi.", "Str.", "Bcell", "Tcell", "Mye."))
DimPlot(seu, group.by = "cellgroup") & NoAxes()
ggsave("cellgroup.png", path = plot_path, width = 4, height = 3, units = "in", dpi = 150)

# add_feat <- "SPARC"
# FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred")) + NoAxes() + NoLegend()
# ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3, height = 3, units = "in", dpi = 150)

# Dot plot
# seu <- JoinLayers(seu)
Idents(seu) <- "cellgroup"
# markers <- FindAllMarkers(seu, only.pos = TRUE)
features <- c("EPCAM", "CLDN18", "DCN", "SPARC", "CD79A", "JCHAIN", "CD3D", "CD3E", "CD14", "CD68")
DotPlot(seu, group.by = "cellgroup", features = features) + RotatedAxis()
ggsave("dotplot.png", path = plot_path, width = 5, height = 3.7, units = "in", dpi = 150)

# save RDS
saveRDS(seu, file = file.path("RDSfiles", "seu_012_cellgroup.RDS"))


# additional plots
add_feat <- "LRG1"
FeaturePlot(seu, features = add_feat, cols = c("lightgrey","darkred"), 
            # pt.size = 1, order = TRUE,
            max.cutoff = "q25", min.cutoff = "q0"
) + NoAxes() 
ggsave(paste0(add_feat, ".png"), path = fp_path, width = 3.6, height = 3, units = "in", dpi = 150)

DotPlot(seu, features = "LRG1", group.by = "cellgroup")
# ggsave("LRG1_dotplot.png", path = plot_path, width = 4, height = 4, units = "in", dpi = 150)