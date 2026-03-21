# continued from 000_
# 2026-03-10
analysis_step <- "010_QC"

# load packages ----
library(tidyverse)
# library(data.table)
# library(Matrix)
library(readxl)
library(Seurat)

# Make directories ----
# fs::dir_create(c("plot", "result", "RDSfiles", "Rscripts"))
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

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
