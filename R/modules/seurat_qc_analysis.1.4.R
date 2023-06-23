############    ----------   SEURAT ANALYSIS    ---------    ############
############    ------------   QC Analysis   ------------    ############

# Version: 1
# Version updates:
#   First version.
# Subversion: 2
# Subversion updates:
# ---> Subversion 1
# Sub-subversion updates:
#   Simply changed the line that calculated the number of cells that would be kept after filtering, output to FilteringVSNoFiltering file. Changed operators to equal those used while actual filtering.
# ---> Subversion 2
# Sub-subversion updates:
#   The patterns to detect mitochondrial and ribosomal genes were modified so that mouse genes could also be detected.
# ---> Subversion 3
# Sub-subversion updates:
#   Forgot to specify.
# ---> Subversion 4
# Sub-subversion updates:
#   Big changes, but all related to the format of the QC plots. Importantly, the program will now split the QC features' distributions according to groups for major variables such as sequencing library ('library.id.tag') and FACS-sorting batch ('chrom.batch.tag'). Notably, this section is only applied if the variables are defined under the names specified within brackets before.

# ---> Dependencies:
#   Seurat and stringr
library(gridExtra)
library(grid)
# ---> Functions dependencies:
#   translate.ids, get.file.name, among others.
# ---> Description:
#   General Quality Control (QC) analysis as part of the overall Seurat analysis, that, up to this version, is composed by the next steps:
#     1. Remove confounding TCR-related genes if requested.
#     2. Caclulate percentage of mitochondrial and ribosomal counts per cell.
#     3. General QC assessment previous to filtering, which includes analysis of the features for QC description and their joint analysis.
#     4. If desired, filtering is applied according to QC thresholds input.
#     5. If applies, General QC assessment after filtering is carried out.
#     6. Post-filtering stats for both cell and gene measurements.
#   A lot of plots are output along the analysis.

# Arguments ------------------------->
feats.cols <- c('#b22222', '#319272', '#0000cd')
names(feats.cols) <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
# Other argument values are taken from the main program.
#   seurat.obj, QC thresholds (min.nCounts, max.nCounts, etc.)
# Value ----------------------------->
# Seurat object post-QC analysis.

### --------------------------- Functions --------------------------- ###
# 1 -------------------------------------------------------------------->
# Name: General QC assessment.
# Description:
# Outputs plots for a general QC assessment which basically depict a description of the features for quality assessment (violin plots with both normal and log scale) and a description of their relationships.
# Arguments:
#   filt.status - Char, either Pre or Post describing the filtering status.
# Value:
# NA

gen.qc.assessment <- function(filt.status, seurat.obj){
  # ---> Original plots.
  # Features for QC description.
  feats.dists <- lapply(X=names(feats.to.assess), FUN=function(feature.name){
    tmp.feature <- feats.to.assess[[feature.name]]
    seurat.obj@meta.data$tmp.col <- tmp.feature
    tmp.col <- feats.cols[tmp.feature]
    tmp.ggplot <- vln.plot(
        seurat.obj=seurat.obj,
        feature=tmp.feature, slot=NULL,
        groups.tag='tmp.col',
        trim.val=TRUE,
        color='color', this.color.scale=tmp.col,
    )
    # Add cutoffs
    tmp.tholds <- QC.tholds[[feature.name]]
    tmp.ggplot <- tmp.ggplot +
        geom_hline(yintercept=tmp.tholds['lower.thold'], col='blue', linetype='dashed')
    tmp.ggplot <- tmp.ggplot +
        geom_hline(yintercept=tmp.tholds['upper.thold'], col='red', linetype='dashed')
    # Final formatting
    tmp.ggplot <- tmp.ggplot +
        labs(title=feature.name) +
        theme_minimal() +
        theme(
            legend.pos='none',
            plot.title=element_text(face='bold', size=8),
            panel.grid=element_blank(),
            axis.line=element_line(size=0.8, color='black'),
            axis.text.x=element_blank(),
            axis.title=element_blank(),
            axis.text.y=element_text(size=unit(10, 'cm'))
        )
    seurat.obj@meta.data$tmp.col <- NULL
    return(tmp.ggplot)
  })
  names(feats.dists) <- names(feats.to.assess)
  # Joint analysis.
  feats.rels <- ggplot(data=seurat.obj@meta.data, aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) +
    geom_point(alpha=0.6, size=0.6) +
    scale_color_manual(scale_name='Mit. genes pct.', limits=c(0, 60)) +
    scale_color_gradient(name='Mit. genes pct.', low='#8b0000', high='#cccc00') + xlab('UMIs no. per cell') + ylab('Genes per cell') +
    labs(title='Features relationships') +
    theme_bw() +
    theme(plot.title=element_text(face='bold', size=8))
  # Add cutoffs
  tmp.tholds <- QC.tholds[['Genes per cell']]
  feats.rels <- feats.rels +
    geom_hline(yintercept=tmp.tholds['lower.thold'], col='blue', linetype='dashed')
  feats.rels <- feats.rels +
    geom_hline(yintercept=tmp.tholds['upper.thold'], col='red', linetype='dashed')
  tmp.tholds <- QC.tholds[['Seq. depth per cell']]
  feats.rels <- feats.rels +
    geom_vline(xintercept=tmp.tholds['lower.thold'], col='blue', linetype='dashed')
  feats.rels <- feats.rels +
    geom_vline(xintercept=tmp.tholds['upper.thold'], col='red', linetype='dashed')
  # Output.
  tmp.file.name <- paste0(QC.path, '/', filt.status, 'filteringFeatures.pdf')
  pdf(file=tmp.file.name)
  grid.arrange(feats.dists[[1]], feats.dists[[2]], feats.dists[[3]], layout_matrix=matrix(c(1,2,3), nrow=1), top=textGrob('Features for QC', gp=gpar(font='bold')))
  print(feats.rels)
  dev.off()

  # ---> Log-scaled plots and UMIs per cell relationships.
  # Log scale.
  feats.dists <- lapply(X=names(feats.to.assess), FUN=function(feature.name){
    tmp.feature <- feats.to.assess[[feature.name]]
    tmp.ggplot <- feats.dists[[feature.name]] +
        scale_y_log10()
    return(tmp.ggplot)
  })
  # Output.
  tmp.file.name <- paste0(QC.path, '/', filt.status, 'filteringFeatures_Log.pdf')
  pdf(file=tmp.file.name)
  grid.arrange(feats.dists[[1]], feats.dists[[2]], feats.dists[[3]], nrow=1, top=textGrob('Features for QC (log10 y scale)', gp=gpar(font='bold')))
  dev.off()

  # ---> Split according to tags of interest.
  tags.of.int <- c(`SeqLib`='library.id.tag', `Batch`='chrom.batch.tag', `FACS`='facs.sorting.batch.tag', `SeqRun`='seq.batch.tag')
  tags.of.int <- tags.of.int[tags.of.int %in% colnames(seurat.obj@meta.data)]
  if(length(tags.of.int)>0){
    for(tag.name in names(tags.of.int)){
        tag.of.int <- tags.of.int[tag.name]
        # Features for QC description.
        feats.dists <- lapply(X=names(feats.to.assess), FUN=function(feature.name){
            tmp.feature <- feats.to.assess[[feature.name]]
            tmp.ggplot <- vln.plot(
                seurat.obj=seurat.obj,
                feature=tmp.feature, slot=NULL,
                groups.tag=tag.of.int,
                trim.val=TRUE,
                color='mean'
            )
            # Add cutoffs
            tmp.tholds <- QC.tholds[[feature.name]]
            tmp.ggplot <- tmp.ggplot +
                geom_hline(yintercept=tmp.tholds['lower.thold'], col='blue', linetype='dashed')
            tmp.ggplot <- tmp.ggplot +
                geom_hline(yintercept=tmp.tholds['upper.thold'], col='red', linetype='dashed')
            # Final formatting
            tmp.ggplot <- tmp.ggplot +
                scale_y_log10() +
                labs(title=feature.name, x=tag.name, y=feature.name, fill='Mean') +
                theme_minimal() +
                theme(
                    legend.pos='bottom',
                    plot.title=element_text(face='bold', size=8),
                    panel.grid=element_blank(),
                    axis.line=element_line(size=0.8, color='black'),
                    axis.text.x=element_text(angle=45),
                    axis.text.y=element_text(size=unit(10, 'cm')),
                    legend.text=element_text(angle=45)
                )
            seurat.obj@meta.data$tmp.col <- NULL
            return(tmp.ggplot)
        })
        names(feats.dists) <- names(feats.to.assess)
        # Output.
        tmp.file.name <- paste0(QC.path, '/', filt.status, 'filteringFeatures_ForTag-', tag.name, '.pdf')
        pdf(file=tmp.file.name)
        for(tmp.val in names(feats.dists)) print(feats.dists[[tmp.val]])
        dev.off()
    }
  }

  return(NA)
}

### ------------------------- Main program -------------------------- ###

seurat.qc.analysis <- function(this.seurat.obj){
  # ---> Antigen receptor (BCR or TCR) genes.
  # Get rid of Ag receptor genes for any kind of chain; just in case this is a 10X 5' experiment. NOTE: This may not always be desired when analyzing cells from a single overall cell type (like when analyzing data only from T cells or B cells) but it comes in handy when having a mixture of different kinds of lymphocytes or lymphocytes with other cell types.
  if(is.five.prime){
    if(feature.id=='name'){
      non.agr.genes <- grep(pattern='^TR[BA]+[VDJC]{1}|^IG[HLK]{1}[VDJC]{1}', x=rownames(this.seurat.obj), perl=TRUE, value=TRUE, invert=TRUE)
    }else{
      non.agr.genes <- feature.info$ensembl[!grepl(x=feature.info$name, pattern='^TR[BA]+[VDJC]{1}|^IG[HLK]{1}[VDJC]{1}', perl=TRUE)]
    }
    this.seurat.obj <- subset(x=this.seurat.obj, features=non.agr.genes)
    cat('The dataset to be considered comes from a five prime version of 10x Genomics. Therefore, TCR genes will be disconted from downstream analyses as, otherwise, they may bias the final results.\n')
  }

  # ---> Mitochondrial and ribosomal genes percentage.
  if(feature.id=='name'){
    # When we've read the feature names since the start...
    this.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(object=this.seurat.obj, pattern = "^MT-|^mt-")
    this.seurat.obj[["percent.rb"]] <- PercentageFeatureSet(object=this.seurat.obj, pattern = "^RPS|^RPL|^Rps|^Rpl")
  }else{
      # Identify ribosomal and mitochondrial genes' ENSEMBL ID.
    mit.genes <- feature.info$ensembl[grepl(x=feature.info$name, pattern="^MT-|^mt-", perl=TRUE)]
    mit.genes <- mit.genes[mit.genes %in% rownames(this.seurat.obj)]
    rib.genes <- feature.info$ensembl[grepl(x=feature.info$name, pattern="^RPS|^RPL|^Rps|^Rpl", perl=TRUE)]
    rib.genes <- rib.genes[rib.genes %in% rownames(this.seurat.obj)]
    # Then, add the feature percentages.
    this.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(object=this.seurat.obj, features=mit.genes)
    this.seurat.obj[["percent.rb"]] <- PercentageFeatureSet(object=this.seurat.obj, features=rib.genes)
  }

  # ---> General QC assessment, prefiltering status.
  gen.qc.assessment(filt.status='Pre', seurat.obj=this.seurat.obj)
  cat('General QC assessment previous to filtering applied.\n\n')

  # ---> Filtering.
  # No matter if cells will be filtered out or not during this analysis, we keep track of the results of that plausible filtering step.
  # Current number of cells vs that after filtering process.
  no.cells.start <- length(Cells(this.seurat.obj))
  no.cells.kept <- length(Cells(subset(this.seurat.obj, subset=nCount_RNA >= min.nCounts & nCount_RNA <= max.nCounts & nFeature_RNA >= min.nFeatures & nFeature_RNA <= max.nFeatures & percent.mt <= max.MP)))
  # Report this comparison.
  tmp.df <- data.frame(no.filtering=no.cells.start, filtering=no.cells.kept)
  tmp.file.name <- get.file.name(gen.path=QC.path, file.desc='FilteringVSNoFiltering', extension='csv')
  write.csv(file=tmp.file.name, x=tmp.df, quote=FALSE)
  # Report also the QC criteria taken into account.
  tmp.df <- Reduce(x=QC.tholds, f=rbind)
  rownames(tmp.df) <- names(QC.tholds)
  tmp.file.name <- get.file.name(gen.path=QC.path, file.desc='QCCriteria', extension='csv')
  write.csv(file=tmp.file.name, x=tmp.df, quote=FALSE)
  # Then, if desired, we filter out the cells considering input parameters.
  if(filter.out){
    # Filter.
    this.seurat.obj <- subset(this.seurat.obj, subset=nCount_RNA >= min.nCounts & nCount_RNA <= max.nCounts & nFeature_RNA >= min.nFeatures & nFeature_RNA <= max.nFeatures & percent.mt <= max.MP)
    cat('Cell filtering applied using\n\tcell counts threshold:', min.nCounts, 'and', max.nCounts, '\n\tgene thresholds:', min.nFeatures, 'and', max.nFeatures, '\n\tMaximum mitochondrial gene percentage per cell:', max.MP, "\n\n")
    # ---> General QC assessment, postfiltering status.
    gen.qc.assessment(filt.status='Post', seurat.obj=this.seurat.obj)
    cat('General QC assessment after filtering applied.\n\n')
  }else{
    cat('There was no filtering of possible low quality cells even though the QC analyses were perfomed. You should see the QC per cluster and check if there is something really meaningful about the low quality cells.\n\n')
  }

  # ---> Post-filtering stats.
  # 1) Measures per cell.
  # a) Features captured per cell.
  feats.per.cell <- this.seurat.obj@meta.data[, 'nFeature_RNA']
  feats.per.cell.summ <- summary(feats.per.cell)
  # b) Sequencing depth (number of UMIs captured per cell).
  seq.depth.per.cell <- Matrix::colSums(x=this.seurat.obj@assays$RNA@counts)
  seq.depth.per.cell.summ <- summary(seq.depth.per.cell)
  # c) Mitochondrial gene percentage per cell.
  mit.pct.per.cell <- this.seurat.obj@meta.data[, "percent.mt"]
  mit.pct.per.cell.summ <- summary(mit.pct.per.cell)
  # 2) Measures per features.
  # a) Average expression
  feats.exp.mean <- Matrix::rowMeans(x=this.seurat.obj@assays$RNA@counts)
  feats.exp.mean.summ <- summary(feats.exp.mean)
  # b) Counts distribution.
  counts.per.gene <- Matrix::rowSums(x=this.seurat.obj@assays$RNA@counts)
  counts.per.gene.summ <- summary(counts.per.gene)
  # Output all summaries.
  tmp.file.name <- get.file.name(gen.path=QC.path, file.desc='/PostFilteringSummaries', extension='csv')
  write.summs(summs=c('feats.per.cell.summ', 'seq.depth.per.cell.summ', 'mit.pct.per.cell.summ', 'feats.exp.mean.summ', 'counts.per.gene.summ'), summs.names=c('Features per cell', 'Sequencing depth', 'Mitochondrial gene percentage',  'Features mean expression accross cells', 'Raw counts per gene'), file.name=tmp.file.name, tmp.env=environment())

  return(this.seurat.obj)
}
