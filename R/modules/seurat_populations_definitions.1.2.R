############    ----------   SEURAT ANALYSIS    ---------    ############
############    ---   Population definition workflow   --    ############

# Dependencies:
#   Seurat
# Functions dependencies:
#   translate.ids, among others.
# Description:
# To describe.
# Arguments ------------------------->
# All argument values are taken from the main program.
#   seurat.obj, norm.method (which should be LogNormalize), FVFs.method, mean.cutoff, prop.cutoff, feats.for.dsa.
# To describe.
# Value ----------------------------->
# Processed seurat object.

# Version 1.1 to fix a bug for name-specific use.

do.pop.definition <- function(resolution){
  cat(paste('Process for resolution', resolution, '\n'))
  clusters.tag <- paste0(ifelse(test=gen.norm.method=='SCTransform', yes='SCT', no='RNA'), '_snn_res.', resolution)
  all.clusters <- levels(seurat.obj@meta.data[, clusters.tag])

  # Identify markers -----------------------------------------------------------
  # Find markers diferentially expressed for each cluster compared to the rest of the cells (Wilcoxon rank sum test).
  seurat.obj.markers <- FindAllMarkers(object=seurat.obj, only.pos=FALSE, min.pct=0.1, test.use='MAST', group.by=clusters.tag)
  if(nrow(seurat.obj.markers)==0) return(NULL)
  cat('Markers (DEGs) for each cluster have been identified.\nAll DEGs saved to a xsv file as well as a Heatmap depicting the top 10 of them will be output.\n')

  # Differential gene expression analysis --------------------------------------
  cat('# Differential gene expression analysis --------------------------------------\n')
  # Create directory for this clusters set and define the appropriate tag in meta data.
  tmp.DGEA.path <- paste0(DGEA.path, '/clusters_set_with_res_', resolution)
  create.dir(tmp.DGEA.path, paste0('DGEA path for clusters set with resolution ', resolution))

  # Top 10 Heatmap
  #top.10.seurat.obj.markers <- seurat.obj.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
  #tmp.file.name <- get.file.name(gen.path=tmp.DGEA.path, file.desc='HeatMapTop10DEGs', extension='pdf', resolution=resolution)
  #pdf(file=tmp.file.name, width=10)
  #print(DoHeatmap(object=seurat.obj, features=top.10.seurat.obj.markers$gene, group.by=clusters.tag) + NoLegend())
  #dev.off()

  # Translate gene name.
  if(feature.id=='ensembl') seurat.obj.markers$gene <- names(translate.ids(ids=seurat.obj.markers$gene))

  # All DEGs output to a csv file.
  tmp.file.name <- get.file.name(gen.path=tmp.DGEA.path, file.desc='DGEsTable', extension='csv', resolution=resolution)
  write.csv(x=seurat.obj.markers, file=tmp.file.name)

  if(do.dea){
    cat(paste0('This analysis will be saved for each cluster evaluated agains the rest of the cells and against each of the rest of the clusters individually.\n'))
    # Output them to a file.
    for(cluster in all.clusters){
      cat('Specific DGEA for cluster', cluster, '\n')
      cat(paste0('Cluster ', cluster, '...\n'))
      cluster.path <- paste0(tmp.DGEA.path, '/cluster', cluster)
      create.dir(dir.path=cluster.path, path.desc=paste0('Cluster ', cluster))
      # Save DGAE vs the rest of the cells.
      cluster.spc.rows <- seurat.obj.markers[, 'cluster'] == cluster
      cluster.markers <- seurat.obj.markers[cluster.spc.rows, -6]
      tmp.file.name <- get.file.name(gen.path=cluster.path, file.desc=paste0('DGEACluster', cluster, 'vsAllCells'), extension='csv', resolution=resolution)
      write.csv(x=cluster.markers, file=tmp.file.name)
    	# Also present a heatmap with the top 10 and the bottom 10 features specifically for the cells in this cluster.
    	top.10.cluster.markers <- as.character(head(cluster.markers[order(cluster.markers$avg_logFC), 'gene'], 10))
    	btm.10.cluster.markers <- as.character(head(cluster.markers[order(cluster.markers$avg_logFC, decreasing=TRUE), 'gene'], 10))
    	# Cells in the cluster in turn.
    	cluster.spc.cells <- seurat.obj@meta.data[, clusters.tag] == cluster
    	cluster.spc.cells <- colnames(seurat.obj@assays$RNA@scale.data)[cluster.spc.cells]
    	# Output heatmap.
    	tmp.file.name <- get.file.name(gen.path=cluster.path, file.desc=paste0('DGEAHeatMapCluster', cluster), extension='pdf', resolution=resolution)
    	pdf(file=tmp.file.name)
    	print(DoHeatmap(object=seurat.obj, features=c(top.10.cluster.markers, btm.10.cluster.markers), cells=cluster.spc.cells))
      dev.off()
    	# Now, against each other cluster.
      rest.clusters <- all.clusters[all.clusters != cluster]
      for(next.cluster in rest.clusters){
          cluster.markers <- FindMarkers(object=seurat.obj, ident.1=cluster, ident.2=next.cluster, only.pos=FALSE, min.pct=0.20)
          tmp.file.name <- get.file.name(gen.path=cluster.path, file.desc=paste0('DGEACluster', cluster, 'vs', next.cluster), extension='csv', resolution=resolution)
          write.csv(x=cluster.markers, file=tmp.file.name)
      }
      cat('Next...\n\n')
    }
    cat('DGEA finished for resolution', resolution, '\n\n')
  }else{
    cat('No DGEA for resolution', resolution, '\n\n')
  }

  # Ciro's significance plot ---------------------------------------------------
  cat('# Ciro\'s significance plot ---------------------------------------------------\n')
  # Define the delta value as the difference of the one cluster minus the rest of cells for the percentage of cells expressing a specific gene.
  seurat.obj.markers$Delta <- 100 * (seurat.obj.markers$pct.1 - seurat.obj.markers$pct.2)
  # Then, we define the significance for each DE gene so that it's negatively correlated with the p-value.
  seurat.obj.markers$Significance <- (-log10(seurat.obj.markers$p_val_adj))
  # If significance is way too high (approximated to an infinite value), this is set to the maximum significance captured + 1.
  seurat.obj.markers$Significance[seurat.obj.markers$Significance == Inf] <- max(seurat.obj.markers$Significance[seurat.obj.markers$Significance != Inf]) + 1
  # Then, just take the top genes as selected by the user.
  top.genes <- as.data.frame(rbindlist(lapply(X=unique(seurat.obj.markers$cluster), FUN=function(tmp.cluster){
    cluster.vals <- seurat.obj.markers[seurat.obj.markers$cluster==tmp.cluster, ]
    to.discard <- grep(pattern="^RPS|^RPL|^MT-", x=cluster.vals$gene)
    if(sum(to.discard)) cluster.vals <- cluster.vals[-to.discard, ] # Removing RPS and RPL genes
    cluster.vals <- cluster.vals[cluster.vals$avg_logFC > 0, ] # In order to keep only positive LFC markers.
    setorder(x=cluster.vals, p_val_adj) # Order according to p value ascending order.
    head(cluster.vals, no.top.sig.feats)
  })))
  # Grid definition.
  tmp.grid <- fitgrid(table(top.genes$cluster))
  # Info per cluster.
  cells.per.cluster <- table(seurat.obj@meta.data[, clusters.tag])
  cells.pct.per.cluster <- as.data.frame(cells.per.cluster*100 / sum(cells.per.cluster))
  colnames(cells.pct.per.cluster) <- c('cluster', 'percentage')
  cells.pct.per.cluster$cells <- cells.per.cluster
  cells.pct.per.cluster$percentage <- round(x=cells.pct.per.cluster$percentage, digits=2)
  cells.pct.per.cluster$clusters.info <- paste0('C', cells.pct.per.cluster$cluster, '-', cells.pct.per.cluster$cells, '-', cells.pct.per.cluster$percentage, '%')
  clusters.info <- cells.pct.per.cluster[, c('cluster', 'clusters.info')]
  # Merge clusters cells info with the cluster stats.
  top.genes <- merge(x=top.genes, y=clusters.info, by='cluster')
  # Plot.
  tmp.ggplot <- ggplot(data=top.genes, aes(x=gene, y=Delta, size=Significance, color=avg_logFC)) +
  labs(title='Significance plot', y="Delta value", x="Gene name", caption="Cluster-Number of cells-Percent of total") +
  geom_point() + facet_wrap(~clusters.info, ncol=tmp.grid[2], scale='free_x') +
  geom_hline(yintercept=0, linetype="dashed", color="red") + theme_minimal() + theme(axis.text.x=element_text(angle=45, hjust=1, size=4))
  # Output significance plot.
  tmp.file.name <- get.file.name(gen.path=significance.path, file.desc=paste0('SignificancePlotForResolution', resolution), extension='pdf')
  pdf(file=tmp.file.name, height=12, width=9)
  print(tmp.ggplot)
  dev.off()
  cat(paste('Significance plot output for resolution', resolution, '\n\n'))

  # Mock return.
  return(NA)
}
