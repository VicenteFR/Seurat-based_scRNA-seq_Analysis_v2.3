############    ----------   SEURAT ANALYSIS    ---------    ############
############    - Clusters Definition and -Specific QC  -    ############

# Version: 1
# Version updates:
#   First stable version.
# Subversion: 3
# Subversion updates:
# ---> Subversion 2
#   We turned the creation of marker feature plots into a linear process (with lapply) instead of a parallel process (with mclapply).
#   Minor bug correction. Minimum number of cells to be depicted on plot is now properly calculated.
# ---> Subversion 3
#   Output violin plots depicting distribution of QC features across clusters with a better resolution.

# ---> Dependencies:
#   Seurat and stringr
library(gridExtra)
library(grid)
library(data.table)
# ---> Functions dependencies:
#   translate.ids, get.file.name, among others.
# ---> Description:
#   Clusters definition (Def) and Clusters-specific QC (QC) process as part of the overall Seurat analysis, that, up to this version, is composed by the next steps:
#     1. (Def) Visualization of markers expression across cells on dimensionality reduction plot, either on UMAP or tSNE (or both). This for both options, all cells or only a subset of cells.
#     2. (Def) Proportion of cells per cluster as barplots for each resolution.
#     3. (Def) Markers expresion distributions as violin plots, plotting with points and without.
#     4. (QC) QC features distributions across clusters for all resolutions.
#     5. (QC) QC features as heatmaps on dim. reduction plots (only for a subset of cells).
#   A lot of plots are output along the analysis.

# Arguments ------------------------->
# Markers as ENSEMBL IDs.
ensembl.markers <- c("ENSG00000010610", "ENSG00000167286", "ENSG00000198851", "ENSG00000160654", "ENSG00000153563", "ENSG00000172116", "ENSG00000110848", "ENSG00000102245", "ENSG00000160683", "ENSG00000113916", "ENSG00000049768", "ENSG00000030419", "ENSG00000107485", "ENSG00000100453", "ENSG00000163600", "ENSG00000111537", "ENSG00000136634", "ENSG00000138684", "ENSG00000134460", "ENSG00000113520", "ENSG00000113525", "ENSG00000083457", "ENSG00000180644", "ENSG00000143365", "ENSG00000020633", "ENSG00000073861", "ENSG00000186891", "ENSG00000008517")
# Markers as simple gene names.
gene.markers <- c('CD4', 'CD8A', 'CD8B', 'CD8B1', 'CD3G', 'CD3E', 'CD3D', 'TRB', 'TRAC', 'CXCR5', 'BCL6', 'FOXP3', 'IL2RA', 'TNFRSF18', 'ICOS', 'IFNG', 'TBX21', 'GATA3', 'RORC', 'PRF1', 'GZMB', 'IL4', 'IL5', 'IL21', 'IL10', 'CD69', 'RUNX3', 'ITGAE') # CD8B1 is not a gene for humans, but its equivalent in mouse gene language does exist as an orthologue of CD8B
gene.markers <- c(gene.markers, str_to_title(string=gene.markers))
# Other argument values are taken from the main program.
#   this.seurat.obj, resolutions, etc.
# Value ----------------------------->
# Seurat object post-QC analysis.

### --------------------------- Functions --------------------------- ###

### ------------------------- Main program -------------------------- ###

clusters.def.and.qc <- function(this.seurat.obj){
  cat('\n\n')
  ### ------------------------ Clusters definition ------------------------- ###
  cat('### ------------------------ Clusters definition ------------------------- ###\n')

  this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
  # ---> Cells to depict in plots.
  # If there's a number of cells defined in the seurat object larger than the maximum required (max.cells.no), we will depict only a random subset of them alongside all of them in different files.
  need.subset <- length(Cells(seurat.obj)) > max.cells.no
  if(need.subset) if(length(Cells(seurat.obj))*0.1 > max.cells.no) max.cells.no <- length(Cells(seurat.obj))*0.1
  if(need.subset) cells.subset <- sample(x=Cells(seurat.obj), size=max.cells.no, replace=FALSE) else cells.subset <- Cells(seurat.obj)

  cat('\n\n')
  # Markers definition -------------------------------------------------------->
  cat('# Markers definition -------------------------------------------------------->\n')

  # Define the kind of marker IDs that will be used as base markers.
  if(feature.id=='ensembl') markers.to.assess <- ensembl.markers else markers.to.assess <- gene.markers
  # If a markers file is provided, the markers defined there will also be considered.
  if(!is.null(markers.file)){
    if(file.exists(markers.file) & grep(x=markers.file, pattern='RData', ignore.case=TRUE)){
      cat('Markers file was provided so that the markers found there as lists will also be included in the analysis.\n')
      markers.env <- new.env()
      load(file=markers.file, envir=markers.env)
      markers.to.assess <- unique(c(markers.to.assess, unlist(lapply(X=ls(markers.env), FUN=function(tmp.list) unlist(get(x=tmp.list, envir=markers.env))))))
      rm(markers.env)
    }
    else{
      cat('Even though a markers file was provided, it didn\'nt seem to have the appropriate format. It should be an RData file and, of course, exist as provided with the absolute path.\n')
    }
  }
  # keep just the markers for which we have levels of expression.
  markers.to.assess <- markers.to.assess[markers.to.assess %in% rownames(this.seurat.obj@assays$RNA)]
  # Translate IDs if needed.
  if(feature.id=='ensembl') markers.to.assess <- translate.ids(ids=markers.to.assess) else names(markers.to.assess) <- markers.to.assess
  cat('Markers to assess:\n', paste(paste0(names(markers.to.assess)[1:length(markers.to.assess)-1], ', ', collapse=''), 'and', names(markers.to.assess)[length(markers.to.assess)], '\n\n'))

  cat('\n\n')
  # Markers visualization ----------------------------------------------------->
  cat("# Markers visualization ----------------------------------------------------->\n")

  markers.path <- paste0(tmp.dim.reduction.path, '/markers')
  create.dir(markers.path, 'Markers')

  # Split markers into different chunks.
  markers.to.assess.list <- split(markers.to.assess, ceiling(seq_along(markers.to.assess)/10))
  # Depict markers for each dimensional reduction method according to chunk.
  cat('Marker will be depicted for each dimensionality reduction map.\nWorking on that...\n\n')
  lapply(X=dim.reduction.methods, FUN=function(method){
    lapply(X=1:length(markers.to.assess.list), FUN=function(i){
    # Define chunk markers.
      tmp.markers.to.assess <- markers.to.assess.list[[i]]
    # Each marker in a separate pdf page.
      # All cells.
      tmp.file.name <- get.file.name(gen.path=markers.path, file.desc=paste0(toupper(method), 'MarkersFeaturePlotsChunk', i), extension='pdf', no.PCs=no.PCs)
      pdf(file=tmp.file.name)
      for(marker in names(tmp.markers.to.assess)){
        print(FeaturePlot(object=this.seurat.obj, features=markers.to.assess[marker], reduction=method) + scale_color_gradientn(colors=this.color.scale) + labs(title=marker, x='UMAP 1', y='UMAP 2'))
      }
      dev.off()
      # Subset of cells.
      if(need.subset){
        tmp.file.name <- get.file.name(gen.path=markers.path, file.desc=paste0(toupper(method), 'MarkersFeaturePlots_CellsSubset_Chunk', i), extension='pdf', no.PCs=no.PCs)
        pdf(file=tmp.file.name)
        for(marker in names(tmp.markers.to.assess)){
          print(FeaturePlot(object=this.seurat.obj, features=markers.to.assess[marker], reduction=method, cells=cells.subset) + scale_color_gradientn(colors=this.color.scale) + labs(title=marker, x='UMAP 1', y='UMAP 2'))
        }
        dev.off()
      }
    # All markers in the same page.
      # All cells.
      tmp.file.name <- get.file.name(gen.path=markers.path, file.desc=paste0(toupper(method), 'AllMarkersFeaturePlotsChunk', i), extension='pdf', no.PCs=no.PCs)
      pdf(file=tmp.file.name, width=20, height=10)
      print(FeaturePlot(object=this.seurat.obj, features=tmp.markers.to.assess, reduction=method, ncol=5, cols=c('#ffdf32', '#ff9a00', '#670000')))
      dev.off()
      # Subset of cells.
      if(need.subset){
        tmp.file.name <- get.file.name(gen.path=markers.path, file.desc=paste0(toupper(method), 'AllMarkersFeaturePlots_CellsSubset_Chunk', i), extension='pdf', no.PCs=no.PCs)
        pdf(file=tmp.file.name, width=20, height=10)
        print(FeaturePlot(object=this.seurat.obj, features=tmp.markers.to.assess, reduction=method, ncol=5, cells=cells.subset, cols=c('#ffdf32', '#ff9a00', '#670000')))
        dev.off()
      }
      # Mock return.
      return(NA)
    })
    # Mock return.
    return(NA)
  })
  # Output a descriptive list of the markers considered during the process.
  markers.list <- lapply(X=names(markers.to.assess.list), FUN=function(chunk.no){tmp.output <- data.frame(gene=names(markers.to.assess.list[[chunk.no]]), chunk=chunk.no); if(feature.id=='ensembl') tmp.output$ensembl.id <- markers.to.assess.list[[chunk.no]]; return(tmp.output)})
  markers.list <- rbindlist(markers.list)
  markers.list$index <- 1:nrow(markers.list)
  tmp.file.name <- get.file.name(gen.path=markers.path, file.desc='MarkersList', extension='csv', no.PCs=no.PCs)
  write.csv(file=tmp.file.name, x=markers.list, quote=FALSE, row.names=FALSE)

  cat('\n\n')
  # Cells proportions per cluster --------------------------------------------->
  cat('# Cells proportions per cluster --------------------------------------------->\n')

  props.per.cluster.path <- paste0(tmp.dim.reduction.path, '/proportions_per_cluster')
  create.dir(props.per.cluster.path, 'Proportions per cluster')

  cat('Also, a barplot depicting the percentage of cells per cluster (according to each resolution value applied) will be output.\n')
  mclapply(X=resolutions, mc.cores=total.workers, FUN=function(resolution){
    # A different cluster tag for each resolution
    clusters.tag <- paste0(ifelse(test=gen.norm.method=='SCTransform', yes='SCT', no='RNA'), '_snn_res.', resolution)
    cells.per.cluster <- table(this.seurat.obj@meta.data[, clusters.tag])
    cells.pct.per.cluster <- as.data.frame(cells.per.cluster*100 / sum(cells.per.cluster))
    colnames(cells.pct.per.cluster) <- c('Cluster', 'Percentage')
    cells.pct.per.cluster$Percentage <- round(x=cells.pct.per.cluster$Percentage, digits=2)
    # Ggplot.
    cell.pct.per.cluster.barplot <- ggplot(data=cells.pct.per.cluster) + aes(x=Cluster, y=Percentage, fill=Cluster) + geom_col(width=0.6) + labs(title='Cell percentage per cluster', caption=paste0('Total amount of cells, ', sum(cells.per.cluster))) +  geom_text(aes(label=Percentage), vjust=1.6, color="black", size=3.5) + theme_minimal() + theme(legend.position='none')
    # Output the ggplot.
    tmp.file.name <- get.file.name(gen.path=props.per.cluster.path, file.desc='CellPctPerClusterBarPlot', extension='pdf', resolution=resolution, no.PCs=no.PCs)
    pdf(file=tmp.file.name)
    print(cell.pct.per.cluster.barplot)
    dev.off()
    # Mock return.
    return(NA)
  })

  cat('\n\n')
  # Markers expression distributions ------------------------------------------>
  cat('# Markers expression distributions ------------------------------------------>\n')

  vln.plots.path <- paste0(tmp.dim.reduction.path, '/expression_distribution')
  create.dir(dir.path=vln.plots.path, path.desc='Violin Plots')

  cat('Outputting marker expression distributions (violin plots)...\n')
  mclapply(X=resolutions, mc.cores=total.workers, FUN=function(resolution){
    # Define resolution-specific reports path.
    tmp.vln.plots.path <- paste0(vln.plots.path, '/res_', resolution)
    create.dir(dir.path=tmp.vln.plots.path, path.desc=paste0('Violin Plots for resolution ', resolution))
    # Define clusters tag.
    clusters.tag <- paste0(ifelse(test=gen.norm.method=='SCTransform', yes='SCT', no='RNA'), '_snn_res.', resolution)
    mclapply(X=1:length(markers.to.assess.list), FUN=function(i){
    # Define chunk markers.
      tmp.markers.to.assess <- markers.to.assess.list[[i]]
      # With points.
      tmp.file.name <- get.file.name(gen.path=tmp.vln.plots.path, file.desc=paste0('MarkersExpDistributionChunk', i), extension='pdf', no.PCs=no.PCs)
      pdf(file=tmp.file.name, width=10)
      for(marker in names(tmp.markers.to.assess)){
        if(feature.id=='ensembl') tmp.marker <- tmp.markers.to.assess[marker] else tmp.marker <- marker
        tmp.ggplot <- vln.plot(seurat.obj=this.seurat.obj, feature=tmp.marker, slot='data', groups.tag=clusters.tag, na.rm=TRUE, color='burst.frequency', trim.val=TRUE, add.points=TRUE)
        print(tmp.ggplot + labs(title=marker, y='Normalized expression value', fill='Burst freq.') + theme_bw())
      }
      dev.off()
      # No points.
      tmp.file.name <- get.file.name(gen.path=tmp.vln.plots.path, file.desc=paste0('MarkersExpDistribution_NoPoints_Chunk', i), extension='pdf', no.PCs=no.PCs)
      pdf(file=tmp.file.name, width=10)
      for(marker in names(tmp.markers.to.assess)){
        if(feature.id=='ensembl') tmp.marker <- tmp.markers.to.assess[marker] else tmp.marker <- marker
        tmp.ggplot <- vln.plot(seurat.obj=this.seurat.obj, feature=tmp.marker, slot='data', groups.tag=clusters.tag, na.rm=TRUE, color='burst.frequency', trim.val=TRUE, add.points=FALSE)
        print(tmp.ggplot + labs(title=marker, y='Normalized expression value', fill='Burst freq.') + theme_bw())
      }
      dev.off()
      # Mock return.
      return(NA)
    })
  })
  # Output a descriptive list of the markers considered during the process.
  tmp.file.name <- get.file.name(gen.path=vln.plots.path, file.desc='MarkersList', extension='csv', no.PCs=no.PCs)
  write.csv(file=tmp.file.name, x=markers.list, quote=FALSE, row.names=FALSE)

  cat('\n\n')
  ### ------------------------ Cluster-specific QC ------------------------- ###
  cat("### ------------------------ Cluster-specific QC ------------------------- ###\n")
  cat("Analysis of the cluster-specific quality control.\n\n")

  # New directory for cluster-specific QC.
  cs.QC.path <- paste0(QC.path, '/cluster_specific')
  create.dir(dir.path=cs.QC.path, path.desc='Cluster-specific QC path')

  # These analysis will be splitted among clusters sets.
  cat('Next processes will be run for each set of clusters according to the different resolution values applied.\n')

  mclapply(X=resolutions, mc.cores=total.workers, FUN=function(resolution){
    cat('### ----------------------- CS QC for resolution', resolution, '------------------------ ###\n')
    tmp.cs.QC.path <- paste0(cs.QC.path, '/resolution_', resolution)
    create.dir(tmp.cs.QC.path, paste0('SC-QC path for resolution', resolution))
    clusters.tag <- paste0(ifelse(test=gen.norm.method=='SCTransform', yes='SCT', no='RNA'), '_snn_res.', resolution)
    # QC measures distributions ----------------------------------------------->
    cat('# QC measures distributions ----------------------------------------------->\n')
    tmp.file.name <- get.file.name(gen.path=tmp.cs.QC.path, file.desc='CSQCVlnPlots', extension='pdf', resolution=resolution)
    pdf(file=tmp.file.name, width=10)
    for(feature in names(feats.to.assess)){
      tmp.ggplot <- vln.plot(seurat.obj=this.seurat.obj, feature=feats.to.assess[[feature]], slot=NULL, groups.tag=clusters.tag, na.rm=TRUE, color='median', adjust.val=1.4, trim.val=TRUE, add.points=FALSE)
      tmp.ggplot <- tmp.ggplot + labs(x='Cluster', y=feature, fill='Median')
      print(tmp.ggplot + theme_bw())
    }
    dev.off()
    cat('Violin plots for number of genes, counts per cell and mitochondrial genes percentage saved for resolution', resolution, '.\n\n')
    # QC measures depicted on dim. reduction plots ---------------------------->
    cat('# QC measures depicted on dim. reduction plots ---------------------------->\n')
    for(method in dim.reduction.methods){
      tmp.file.name <- get.file.name(gen.path=tmp.cs.QC.path, file.desc=paste0('CSQCDimRed', toupper(method)), extension='pdf', resolution=resolution, no.PCs=no.PCs)
      pdf(file=tmp.file.name, width=12)
      for(feature in feats.to.assess){
        plot.1 <- DimPlot(this.seurat.obj, reduction=method, dims=1:2, group.by=clusters.tag, label=TRUE, cells=cells.subset)
        plot.2 <- FeaturePlot(object=this.seurat.obj, features=feature, dims=1:2, reduction=method, cells=cells.subset)
        print(CombinePlots(plots=list(plot.1, plot.2)))
      }
      dev.off()
    }
    cat('Dimensionality reduction plots depicting QC measures output.\n')
    # Mock return.
    return(NA)
  })
  cat("Cluster-specific QC completed!\n\n")

  # Mock return.
  return(NA)
}
