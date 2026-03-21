library(DoubletFinder)

for (sn in names(seu_list)) {
  seu <- seu_list[[sn]]
  ncells <- ncol(seu)
  
  message("Processing: ", sn, " (n=", ncells, ")")
  
  # Skip small samples
  if (ncells < 100) {
    message("  Skipping DoubletFinder: too few cells")
    seu$doublet_finder <- NA
    seu_list[[sn]] <- seu
    next
  }
  
  # Standard preprocessing
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  
  n_feats <- length(VariableFeatures(seu))
  max_pcs <- min(ncells, n_feats) - 1
  
  if (is.na(max_pcs) || max_pcs < 5) {
    message("  Skipping DoubletFinder: too few usable PCs")
    seu$doublet_finder <- NA
    seu_list[[sn]] <- seu
    next
  }
  
  pcs_use <- min(10, max_pcs)
  
  seu <- RunPCA(seu, npcs = pcs_use, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:pcs_use, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = 1:pcs_use, verbose = FALSE)
  seu <- FindClusters(seu, resolution = 0.1, verbose = FALSE)
  
  # pK identification
  sweep.list <- paramSweep(seu, PCs = 1:pcs_use)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  if (nrow(bcmvn) == 0) {
    message("  Skipping DoubletFinder: pK search failed")
    seu$doublet_finder <- NA
    seu_list[[sn]] <- seu
    next
  }
  
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric), ]
  optimal.pk <- as.numeric(as.character(bcmvn.max$pK))
  
  # Homotypic proportion
  annotations <- seu$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  
  # crude 10x estimate
  nExp.poi <- round((0.01 * ncells / 1000) * ncells)
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  nExp.poi.adj <- max(1, nExp.poi.adj)
  
  # run DoubletFinder
  seu <- doubletFinder(
    seu = seu,
    PCs = 1:pcs_use,
    pK = optimal.pk,
    nExp = nExp.poi.adj
  )
  
  # rename last metadata column
  colnames(seu@meta.data)[ncol(seu@meta.data)] <- "doublet_finder"
  
  seu_list[[sn]] <- seu
}

# rm(bcmvn, bcmvn.max, metadata, seu, sweep.list, sweep.stats, annotations, homotypic.prop, i, nExp.poi, nExp.poi.adj, optimal.pk)

# visualize doublet
# for(sn in names(seu_list)){
#   seu <- seu_list[[sn]]
#   
#   if (sum(is.na(seu$doublet_finder)) > 0) {
#     next
#   }
#   
#   dim <- DimPlot(seu, group.by = "doublet_finder") + NoAxes()
#   ggsave(paste0("DoubletFinder_", sn, ".png"), plot = dim,
#          path = plot_path,
#          width = 4, height = 3, units = "in", dpi = 150)
#   vln <- VlnPlot(seu_list[[sn]], features = c("nFeature_RNA", "nCount_RNA"), group.by = "doublet_finder")
#   ggsave(paste0("DoubletFinder_vln_", sn, ".png"), plot = vln,
#          path = plot_path,
#          width = 5, height = 3, units = "in", dpi = 150)
#   rm(dim, vln, sn)
# }
