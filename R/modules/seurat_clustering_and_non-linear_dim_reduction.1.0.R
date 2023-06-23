############    ----------   SEURAT ANALYSIS    ---------    ############
############    -----   Clustering and non-linear   -----    ############
############    ------   dimensionality reduction   -----    ############

# ---> WIP

# Dependencies:
#   Seurat, Harmony and stringr
# Functions dependencies:
#   translate.ids, among others.
# Description:
#   Non-linear dimensionality reduction (UMAP and/or tSNE) and clustering.
# When 'do.harmony' is NULL, PCs from PCA are used as the reduction components for both processes. In the opposite case, harmony components are used.

# Arguments ------------------------->
# All argument values are taken from the main program.
#   seurat.obj, do.harmony.
# Value ----------------------------->
# Seurat object processed accordingly.

############    --------   Program definitions   --------    ############
# ---> This program's version-related details:
# Version: 1
# Version updates:
#   First version.
# Subversion: 0
#   First subversion.


dim.reduction.clustering <- function(){
  init.reduction <- if(is.null(do.harmony)) 'pca' else 'harmony'
  ### ------------------------ Cells clustering ------------------------ ###
  cat('### ------------------------ Cells clustering ------------------------ ###\n')
  cat('Finding clusters...\n')
  seurat.obj <- FindNeighbors(object=seurat.obj, dims=1:no.PCs, reduction=init.reduction)
  for(resolution in resolutions) seurat.obj <- FindClusters(object=seurat.obj, resolution=resolution)
  cat('All values provided were used.\n')

  ### ------------------ Non-linear dimensional reduction ------------------ ###
  cat('### ------------------ Non-linear dimensional reduction ------------------ ###\n')
  cmm.plots.path <- paste0(tmp.dim.reduction.path, '/UMAP_tSNE')
  create.dir(cmm.plots.path, 'Common plots')

  # ---> Cells to depict in plots.
  # If there's a number of cells defined in the seurat object larger than the maximum required (max.cells.no), we will depict only a random subset of them alongside all of them in different files.
  need.subset <- length(Cells(seurat.obj)) > max.cells.no
  if(need.subset) if(length(Cells(seurat.obj))*0.1 > max.cells.no) max.cells.no <- length(Cells(seurat.obj))*0.1
  if(need.subset) cells.subset <- sample(x=Cells(seurat.obj), size=max.cells.no, replace=FALSE) else cells.subset <- Cells(seurat.obj)

  # 1. UMAP --------------------------------------------------------------------
  cat('1. UMAP --------------------------------------------------------------------\n')
  cat('Running UMAP...\n\n')
  if('umap' %in% dim.reduction.methods){
    seurat.obj <- RunUMAP(seurat.obj, dims=1:pcs.to.comp, reduction=init.reduction, n.neighbors=15, min.dist=0.1, spread=1)
  }else{
    cat('UMAP will be skipped.\n')
  }

  # 2. tSNE --------------------------------------------------------------------
  cat('2. tSNE --------------------------------------------------------------------\n')
  cat('Running tSNE...\n\n')
  if('tsne' %in% dim.reduction.methods){
    seurat.obj <- RunTSNE(object=seurat.obj,  perplexity=100, dims=1:no.PCs, tsne.method="FIt-SNE", fast_tsne_path='/mnt/BioHome/ciro/bin/FIt-SNE2/bin/fast_tsne')
  }else{
    cat('tSNE will be skipped.\n')
  }

  # Output the dimensionality reduction plots (with each method) grouping by clusters according to each resolution value applied.
  for(resolution in resolutions){
    clusters.tag <- paste0(ifelse(test=gen.norm.method=='SCTransform', yes='SCT', no='RNA'), '_snn_res.', resolution)
    for(method in dim.reduction.methods){
      # All cells.
      tmp.file.name <- get.file.name(gen.path=cmm.plots.path, file.desc=paste0(toupper(method), '_Clusters'), extension='pdf', no.PCs=no.PCs, resolution=resolution)
      pdf(file=tmp.file.name)
      print(DimPlot(seurat.obj, reduction=method, dims=1:2, group.by=clusters.tag, label=TRUE))
      dev.off()
      if(need.subset){
        # Subset of cells.
        tmp.file.name <- get.file.name(gen.path=cmm.plots.path, file.desc=paste0(toupper(method), '_CellsSubset_Clusters'), extension='pdf', no.PCs=no.PCs, resolution=resolution)
        pdf(file=tmp.file.name)
        print(DimPlot(seurat.obj, reduction=method, dims=1:2, group.by=clusters.tag, label=TRUE, cells=cells.subset))
        dev.off()
      }
    }
  }

  return(seurat.obj)
}
