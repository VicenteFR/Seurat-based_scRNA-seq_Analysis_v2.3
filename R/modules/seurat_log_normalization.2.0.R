############    ----------   SEURAT ANALYSIS    ---------    ############
############    ------   Normalization workflow   -------    ############

# ---> WIP

# Dependencies:
#   Seurat and stringr
# Functions dependencies:
#   translate.ids, among others.
# Description:
#   General normalization workflow as part of the seurat analysis. This is as follows:
#   1. Cell normalization: Normalization per cell (column of the seurat object).
#   2. Feature selection: Selection of most variable features (MVFs). Seurat provides many methods in order to standardize the original variance of each gene/feature in order to select only the most variable ones. feats.for.dsa variable indicates either the percentage of variance desired to be covered for the MVFs or the amount of top variable features to be selected (at least larger than 100 in this case).
#   Also, for a feature to be selected, this must be expressed at a valid level (threshold indicated in mean.cutoff) in a minimum number of cells (indicated in prop.cutoff).
#   3. Feature normalization: Scaling for the values of each gene (mean 0 and standard deviation equals 1).
#   A lot of plots are output along the analysis.

# Arguments ------------------------->
# All argument values are taken from the main program.
#   seurat.obj, norm.method (which should be LogNormalize), FVFs.method, mean.cutoff, prop.cutoff, feats.for.dsa.
# Value ----------------------------->
# List containing:
#   Idx. 1) Processed seurat object.
#   Idx. 2) Chracter vector with final variable features.

############    --------   Program definitions   --------    ############
# ---> This program's version-related details:
# Version: 2
# Version updates:
#   We open the possibility to regress out the signal from a list of multiple genes (n>10) as a signature score (calculated with the function 'AddModuleScore') instead of individual variables.
#   First version.
# Subversion: 0
# --> Subversion: 0
#   * Subversion updates: We open the possibility to regress out the signal from a list of multiple genes (n>100) as a signature score (calculated with the function 'AddModuleScore') instead of individual variables. Furthermore, if more variables should be regressed out as well and are defined independently in the metadata, they will be taken into account along the gene signature.
# Associations:
# -- (1)
#   Description: Same program, last subversion.
#   Version: 1.1
#   Path: /home/vfajardo/scripts/cell_type_classification/module_definition/seurat_log_normalization.1.1.R


do.normalization <- function(){
  # 1. Cell normalization ---------------------------------------------
  cat('1. Cell normalization ---------------------------------------------\n')
  seurat.obj <- NormalizeData(object=seurat.obj, normalization.method=norm.method) # Similar to CPM but with a different scale factor.
  cat('Measured data normalized based on ', norm.method, ' method.\n\n')

  # 2. Feature selection ---------------------------------------------
  cat('2. Feature selection ---------------------------------------------\n')
  # 2.1 Seurat identified variable features --------------------------
  cat('# 2.1 Seurat identified variable features --------------------------\n')
  seurat.obj <- FindVariableFeatures(object=seurat.obj, selection.method=FVFs.method) # Find most variable genes.
  cat('Most variable genes (identified with', FVFs.method, 'method) will be taken into account for downstream analyses.\n')
  top.var.feats <- head(VariableFeatures(seurat.obj), 25) # 10 highest variable genes.
  if(feature.id=='ensembl'){
    top.var.feats <- translate.ids(top.var.feats)
  }else{
    names(top.var.feats) <- top.var.feats
  }
  cat('Here are the top 10 most variable genes:\n', names(top.var.feats), '\n')

  # 2.2 Filtering out first variable features picked -----------------
  cat('# 2.2 Filtering out first variable features picked -----------------\n')
  # Add a mean cutoff (1CPM for 10X data).
  # Previouly done like this:
  #   has.min.mean.exp <- seurat.obj@assays$RNA@meta.features$vst.mean >= mean.cutoff
  feats.mean <- Matrix::rowMeans(seurat.obj@assays$RNA@data)
  seurat.obj@assays$RNA@meta.features$raw.mean <- feats.mean
  has.min.mean.exp <- seurat.obj@assays$RNA@meta.features$vst.mean > mean.cutoff
  seurat.obj@assays$RNA@meta.features$has.min.mean.exp <- has.min.mean.exp
  # Since previously we've already filtered the low quality cells out, we can define the appropriate amount of cells that a gene has to show expression (counts/UMIs>0) for it to be considered as candidate to variable feature or not according to the proportion cutoff for this data set given the one indicated by the user.
  ##### ---> Deprecated.
  #cells.min.no.cutoff <- round(x=prop.cutoff*ncol(seurat.obj), digits=0)
  #has.min.no.cells <- Matrix::rowSums(seurat.obj@assays$RNA@data > 0) > cells.min.no.cutoff
  #seurat.obj@assays$RNA@meta.features$has.min.no.cells <- has.min.no.cells

  # 2.3 Cummulative variance analysis --------------------------------
  # 2.4 Set of final variable features -------------------------------
  cat('# 2.3 Cummulative variance analysis --------------------------------\n')
  cat('# 2.4 Set of final variable features -------------------------------\n')
  cat('Analysis for both kinds of variance, raw and standardized.\nNotice that high levels out of the total variance are reached almost immediately for raw variance.\n')
  # Declare kinds of variance we'll look at.
  # DEVELOPER, NOTICE: Here the order of the kinds of variance we're analyzing is extremely important.
  # The standardized variance must be the last to be assessed since from that we'll take the real final variable features for downstream analysis.
  types.of.var <- c('vst.variance', 'vst.variance.standardized')
  names(types.of.var) <- c('RawVariance', 'StandardizedVariance') # Notice here the order important since the data frame for raw variance will still be used.
  # Take features-vst variance relationship
  for(type.of.var in names(types.of.var)){
    # ---------------------- kind of variance.
    cat('# ----------------------', type.of.var, '\n')
    variance.df <- seurat.obj@assays$RNA@meta.features[, ]
    # Filter out features not passing the cutoffs established.
    variance.df <- variance.df[variance.df$has.min.mean.exp, ]
    # In order to set a default name for the variance field.
    #variance.df <- data.frame(variance=variance.df[, types.of.var[type.of.var]], features=row.names(variance.df), row.names=row.names(variance.df), raw.mean=variance.df[, 'raw.mean'])
    variance.df$features <- row.names(variance.df)
    # Rank genes according to their variance...
    variance.df <- variance.df[order(variance.df[, types.of.var[type.of.var]], decreasing=TRUE), ]
    # Get their cummulative variance pct.
    variance.df$cumm.var.pct <- Reduce(x=variance.df[, types.of.var[type.of.var]]*100/sum(variance.df[, types.of.var[type.of.var]]), f=sum, accumulate=TRUE)
    variance.df$cumm.var.pct <- round(x=variance.df$cumm.var.pct, digits=2)
    # For plotting purposes, we add a ranking column to the data frame.
    variance.df$ranking <- 1:nrow(variance.df)
    # ---------> Take the final variable features.
    # If we're looking at the Standardized Variance, we'll keep the genes that will be used for downstream analysis according to the general option.
    # If it is a raw number to consider...
    if(feats.for.dsa > 100){
      variance.df$variable.feature <- variance.df$ranking<feats.for.dsa
    }else{
    # If it is a percentage.
      variance.df$variable.feature <- variance.df$cumm.var.pct<feats.for.dsa
    }
    # We take the variable features as final for downstream analysis if this iteration is triggered for standardized variance.
    if(type.of.var=='StandardizedVariance') final.variable.feats <- row.names(variance.df)[variance.df$variable.feature]
    # ---------> Cummulative variance plot.
    # Then, we create the plot. In it, we depict variable and non-variable features identified previously.
    no.var.feats.kept <- sum(variance.df$variable.feature)
    cummulative.var.plot <- ggplot(data=variance.df, aes(x=ranking, y=cumm.var.pct, color=variable.feature)) + geom_point() + labs(title='Cummulative variance', subtitle=paste0(type.of.var, ';', no.var.feats.kept, ' var. genes for downstream analyses.')) + xlab('Ranked genes') + ylab('Percentage of cummulative variance') + scale_color_manual(values=c('darkred', 'gold1')) + geom_vline(xintercept=no.var.feats.kept, linetype='dashed', color='red')
    # Output cummulative variance plot.
    tmp.file.name <- get.file.name(gen.path=feat.select.path, file.desc=paste0('Cummulative', type.of.var), extension='pdf')
    pdf(file=tmp.file.name)
    print(cummulative.var.plot)
    dev.off()
    tmp.file.name <- get.file.name(gen.path=feat.select.path, file.desc=paste0('TableCummulative', type.of.var), extension='csv')
    write.csv(file=tmp.file.name, x=variance.df)
    # Identify at which point 95% of variance is found.
    cat('Here are the amounts of genes necessary to reach certain variance percentages.\n\nVariance percent\tAmount of genes:\n')
    for(tmp.var in c(25, 50, 75, 90, 95)){
      tmp.no.feats <- head(variance.df[variance.df$cumm.var.pct > tmp.var, 'ranking'], 1)
      cat(as.character(tmp.var), '\t', as.character(tmp.no.feats), '\n')
    }
    cat('\n')
  }
  cat('\n\n')

  # 2.5 Mean by variance analysis ------------------------------------
  cat('# 2.5 Mean by variance analysis ------------------------------------\n')
  # Add a column to meta features names 'final.variable' as a summary depicting if the gene was selected as variable by seurat and it passess both cutoffs, for minimum proportion of cells and for minimum mean expression.
  is.final.variable.feature <- row.names(seurat.obj@assays$RNA@meta.features) %in% final.variable.feats
  # For plotting purposes, give them the value of 'Variable' or 'Non variable'
  seurat.obj@assays$RNA@meta.features$final.variable <- ifelse(test=is.final.variable.feature, yes='Variable', no='Non variable')
  # Ggplot.
  mean.by.var.plot <- ggplot(data=seurat.obj@assays$RNA@meta.features, aes(x=log10(vst.mean), y=vst.variance.standardized)) +
  geom_point(aes(color=final.variable), size=0.5) +
  labs(title='Genes\' mean by variance plot', x='log10 Average Expression', y='Standardized Variance') + scale_color_manual(values=c('black', 'firebrick')) + geom_vline(xintercept=log10(mean.cutoff), color='khaki1', size=1.5)
  tmp.file.name <- get.file.name(gen.path=feat.select.path, file.desc='MeanByVariancePlot', extension='pdf')
  # Check if top variable features passed the expression cutoff.
  top.var.feats.kept <- top.var.feats[top.var.feats %in% final.variable.feats]
  pdf(file=tmp.file.name)
  print(LabelPoints(plot=mean.by.var.plot, points=top.var.feats.kept, labels=names(top.var.feats.kept), repel=TRUE))
  dev.off()
  cat('Plot depicting the relation log10 expression mean by variance has been output. Cutoff is shown as well as the varibale features picked are marked.')

  # 2.6 Final variable features picked -------------------------------
  cat('# 2.6 Final variable features picked -------------------------------\n')
  # Set final variable features.
  VariableFeatures(seurat.obj) <- final.variable.feats
  # Identify what identified variable features we are throwing away due to cutoffs.
  variable.features.thrown <- rownames(seurat.obj@assays$RNA@meta.features)[seurat.obj@assays$RNA@meta.features$vst.variable & !seurat.obj@assays$RNA@meta.features$has.min.mean.exp]
  #base::setdiff(x=VariableFeatures(seurat.obj), y=final.variable.feats)
  cat('From the variable features identified before that have passed the established cutoffs, only those asked will be used (total', as.character(length(final.variable.feats)), ').\n')
  if(feature.id=='ensembl') variable.features.thrown <- translate.ids(ids=variable.features.thrown) else names(variable.features.thrown) <- variable.features.thrown
  # Then, output them.
  tmp.file.name <- get.file.name(gen.path=feat.select.path, file.desc='VariableFeaturesThrownAway', extension='txt')
  write(file=tmp.file.name, x=names(variable.features.thrown))
  tmp.file.name <- get.file.name(gen.path=feat.select.path, file.desc='VariableFeaturesKept', extension='txt')
  write(file=tmp.file.name, x=as.character(names(final.variable.feats)))
  cat('They\'ve also been output to a file, along with the list of variable features thrown away because of the cutoff. You should chek them, too.\n\n')

  # 3. Gene normalization ---------------------------------------------
  cat('3. Gene normalization ---------------------------------------------\n')
  # ---> Cell cycle scoring.
  if(is.null(regress.cc)){
    cat('No cell cycle scoring to be regressed out.\n')
  }else{
    switch(EXPR=regress.cc,
    'human'={
      cat('Cell cycle scoring for human signature to be regressed out.\n')
      # @ Get scoring if requested for human data.
      seurat.obj <- CellCycleScoring(object=seurat.obj, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=FALSE)
      vars.to.regress <- c(vars.to.regress, c("S.Score", "G2M.Score"))
    },
    'mouse'={
      cat('Cell cycle scoring for mouse signature to be regressed out.\n')
      # @ Get scoring if requested for mouse data.
      #       Dependency needed specifically for this task.
      library(biomaRt)
      #       Obtain mouse gene seignature.
      ensembl.human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
      s.genes <- getBM(attributes='mmusculus_homolog_associated_gene_name',
          filters='hgnc_symbol', # HUGO Gene Nomenclature Comittee Symbol.
          values=cc.genes$s.genes,
          mart=ensembl.human)[, 1]
      g2m.genes <- getBM(attributes='mmusculus_homolog_associated_gene_name',
          filters='hgnc_symbol', # HUGO Gene Nomenclature Comittee Symbol.
          values=cc.genes$g2m.genes,
          mart=ensembl.human)[, 1]
      #       Obtain cell cycle scoring.
      seurat.obj <- CellCycleScoring(object=seurat.obj, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
      vars.to.regress <- c(vars.to.regress)
      rm(ensembl.human)
    },
    stop('No valid option about cell cycle signature info regression. See valid options for flag \'RegressCC\'.\n')
    )
  }

  # ---> Further module scoring
  if(length(vars.to.regress)>100){
    # Reserve the variables saved in the metadata that may be of interest.
    reserved.vars <- vars.to.regress[vars.to.regress %in% colnames(seurat.obj@meta.data)]
    vars.to.regress <- setdiff(vars.to.regress, reserved.vars)
    # Check whether the variables not defined in the metadata are in the seurat object's features. We must have at least 10 features to calculate a signature score.
    vars.to.regress <- vars.to.regress[vars.to.regress %in% row.names(seurat.obj)]
    tmp.check <- length(vars.to.regress) > 10
    if(!tmp.check){
      tmp.warn <- paste0('Various features were requested to be regressed out, many of which may not be defined in the dataset and the ones that are indeed defined (listed below) are not enough (n<10) to calculate a proper signature score to be regressed out. Therefore, though they are defined in the dataset, they will be ignored overall.\n', paste0(vars.to.regress, collapse=', '), '\n')
      warning(tmp.warn)
      vars.to.regress <- if(length(reserved.vars)>0) reserved.vars else NULL
    }else{
      tmp.mssg <- paste0('The following variables, originally requested to be regressed out as individual variables, will be used to infer a single gene signature and such a signature will be regressed out when normalizing the data.\n', paste0(vars.to.regress, collapse=', '), '\n')
      cat(tmp.mssg)
      seurat.obj <- AddModuleScore(object=seurat.obj, features=list(vars.to.regress), name='signature.to.regress')
      colnames(seurat.obj@meta.data)[length(colnames(seurat.obj@meta.data))] <- 'signature.to.regress'
      vars.to.regress <- c(reserved.vars, 'signature.to.regress')
    }
  }else{
    vars.to.regress <- grep(x=vars.to.regress, pattern='-', invert=TRUE, value=TRUE)
  }

  # ---> Scaling and regression.
  # z-scores for all genes.
  seurat.obj <- ScaleData(object=seurat.obj, vars.to.regress=vars.to.regress, block.size=2000, verbose=TRUE)

  # ---> Final ouptut.
  seurat.obj <- list(seurat.obj, final.variable.feats)
  names(seurat.obj) <- c('seurat.obj', 'var.feats')
  return(seurat.obj)
}
