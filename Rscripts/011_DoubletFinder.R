for (i in 1:length(seu_list)) {
  # Pre-process seurat object with standard seurat workflow
  seu <- NormalizeData(seu_list[[i]])
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:10)
  seu <- FindNeighbors(object = seu, dims = 1:10)              
  seu <- FindClusters(object = seu, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(seu, PCs = 1:10)
  sweep.stats <- summarizeSweep(sweep.list)
  bcmvn <- find.pK(sweep.stats)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round((0.01 * ncol(seu) / 1000) * ncol(seu)) ## doublet formation rate - tailor for your dataset, see below
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  seu <- doubletFinder(seu = seu, 
                       PCs = 1:10, 
                       pK = optimal.pk,
                       nExp = nExp.poi.adj)
  metadata <- seu@meta.data
  colnames(metadata)[ncol(metadata)] <- "doublet_finder"
  seu@meta.data <- metadata 
  
  # save
  seu_list[[i]] <- seu
}

rm(bcmvn, bcmvn.max, metadata, seu, sweep.list, sweep.stats, annotations, homotypic.prop, i, nExp.poi, nExp.poi.adj, optimal.pk)

# visualize doublet
for(i in 1:length(seu_list)){
  dim <- DimPlot(seu_list[[i]], group.by = "doublet_finder") + NoAxes()
  ggsave(paste0("DoubletFinder_", as.character(sample_name[i]), ".png"), plot = dim, 
         path = plot_path, 
         width = 4, height = 3, units = "in", dpi = 150)
  vln <- VlnPlot(seu_list[[i]], features = c("nFeature_RNA", "nCount_RNA"), group.by = "doublet_finder")
  ggsave(paste0("DoubletFinder_vln_", as.character(sample_name[i]), ".png"), plot = vln, 
         path = plot_path, 
         width = 5, height = 3, units = "in", dpi = 150)
  rm(dim, vln, i)
}