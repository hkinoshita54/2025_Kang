# continued from 000_
# 2026-03-21

# Settings ----
analysis_step <- "010_QC"
dataset <- "GSE206785"

## load packages
library(tidyverse)
library(Seurat)
library(openxlsx2)

## Make directories
plot_path <- file.path("plot", analysis_step)
fp_path <- file.path(plot_path, "feature_plot")
res_path <- file.path("result", analysis_step)
fs::dir_create(c(plot_path, res_path, fp_path))

# Create a Seurat object with patient data ----
pts_data <- read_xlsx("data/13059_2022_2828_MOESM2_ESM.xlsx", start_row = 3)
pts_data$`Donor ID` <- sub(pattern = "D", replacement = "", pts_data$`Donor ID`)
pts_data <- pts_data[,-c(1,3,4,21:24,26)] # remove unecessary columns
pts_data <- pts_data %>% select(
  patient = `Donor ID`,
  Lauren = `Lauren Classification`,
  TCGA = `TCGA Molecular Grouping`,
  ACRG = ACRG
)

gsm_sample <- read_xlsx("data/gsm_sample.xlsx") # this table is created manually from geo web page
gsm_sample$Sample <- sub(".*D(\\d+)_(\\w)$", "\\1\\2", gsm_sample$Sample)

meta <- read_delim("data/GSE206785_metadata.txt.gz")
### > cell barcodes are not provided in this meta data, so they are not usable
### > keep the information on Platform (SC3P or SC5P)
meta <- meta %>% select(Sample, Platform) %>% distinct() 

mtx <- readRDS("data/GSE206785_raw_matrix.rds")
seu <- CreateSeuratObject(counts = t(mtx), min.cells = 3, min.features = 200)
seu$orig.ident <- colnames(seu) %>% strsplit("-") %>% lapply("[",3) %>% unlist()
seu$dataset <- dataset
seu[[]] <- left_join(seu[[]], gsm_sample, by = c("orig.ident" = "Sample"))
seu[[]] <- left_join(seu[[]], meta, by = c("orig.ident" = "Sample"))
seu$patient <- gsub(pattern = "[A-Z]", replacement = "", seu$orig.ident)
seu$tissue_type <- gsub("[0-9]", "", seu$orig.ident)
seu[[]] <- left_join(seu[[]], pts_data, by = c("patient" = "patient"))

seu[[]] <- mutate_if(seu[[]], is.character, as.factor)
str(seu[[]])

saveRDS(seu, "RDSfiles/seu_010_raw.RDS")

# initial QC ----
Idents(seu) <- "orig.ident"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
ggsave("QC_vln_unfil.png", path = plot_path, width = 25, height = 3, units = "in", dpi = 150)
### > looks like it is already filtered based on percent.mt
table(seu$orig.ident)
# seu <- subset(seu, subset = orig.ident == "180305N", invert = TRUE) # remove sample with only 2 cells

# Create Seurat object list for DoubleFinder ----
seu_cp <- seu
seu_list <- SplitObject(seu, split.by = "orig.ident")

# DoubletFinder ----
source(file.path("Rscripts", "DoubletFinder.R"))

## merge the list ----
seu <- merge(x = seu_list[[1]], y = seu_list[2:length(seu_list)])
# rm(seu_list)
names(seu[[]])
seu@meta.data <- seu@meta.data[,c(1:11,15)]
seu[[]] <- mutate_if(seu[[]], is.character, as.factor)


# Filter doublet ----
seu <- subset(seu, subset = doublet_finder == "Doublet", invert = TRUE)
Idents(seu) <- "orig.ident"
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
ggsave("QC_vln_singlet.png", path = plot_path, width = 25, height = 3, units = "in", dpi = 150)

# Save filtered object
rna_layers <- Layers(seu[["RNA"]])
keep_layers <- rna_layers[grepl("^counts.", rna_layers)]
seu[["RNA"]] <- subset(seu[["RNA"]], layers = keep_layers)
seu <- JoinLayers(seu)
saveRDS(seu, file.path("RDSfiles", "seu_011_singlet.RDS"))
