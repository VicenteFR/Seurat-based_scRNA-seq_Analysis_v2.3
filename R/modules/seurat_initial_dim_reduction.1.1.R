############    ----------   SEURAT ANALYSIS    ---------    ############
############    ------   Normalization workflow   -------    ############
############    -------   Initial dim. reduction   ------    ############

# Version: 1
# Version updates:
#   First version.
# Subversion: 1
# Subversion updates:
# ---> Subversion 1
# Sub-subversion updates:
#   Harmony broke when using the latest versions of its dependencies up to June 2023. The problem has to do with one of Seurat's functions , Key. Thus, I had to modify this function directly within the program. It was an easy fix, by defining the key directly instead of using Seurat's function. See line 82 of the v1.1 of this script. This was actually done based on recommendations froms: https://github.com/immunogenomics/harmony/issues/187

# Dependencies:
#   Seurat, Harmony and stringr
# Functions dependencies:
#   translate.ids, among others.
# Description:
#   Initial dimensionality reduction thrugh these processes:
#   1. PCA: Standard PC calculation and further strategies to aid pick the best number of PCs for downstream analysis.
#   2. Harmony-based dimensionality reduction: Batch-effect removal and calculation of components based on Harmony.
# When 'do.harmony' is NULL, the second process is not performed. When not NULL, do.harmony indicates the tags defined in the seurat object that list the batches desired to be removed.

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
#   Subversion updates:
#   First subversion.


RunHarmony <- function (object, group.by.vars, reduction = "pca", dims.use = NULL,
    theta = NULL, lambda = NULL, sigma = 0.1, nclust = NULL,
    tau = 0, block.size = 0.05, max.iter.harmony = 10, max.iter.cluster = 20,
    epsilon.cluster = 1e-05, epsilon.harmony = 1e-04, plot_convergence = FALSE,
    verbose = TRUE, reference_values = NULL, reduction.save = "harmony",
    assay.use = NULL, project.dim = TRUE, ...)
{ 
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Running Harmony on a Seurat object requires Seurat")
    }
    assay.use <- Seurat::DefaultAssay(object)
    if (reduction == "pca" && !reduction %in% Seurat::Reductions(object = object)
) {
        if (isTRUE(x = verbose)) {
            message("Harmony needs PCA. Trying to run PCA now.")
        }
        object <- tryCatch(expr = Seurat::RunPCA(object = object,
            assay = assay.use, verbose = verbose, reduction.name = reduction),
            error = function(...) {
                stop("Harmony needs PCA. Tried to run PCA and failed.")
            })
    }
    if (!reduction %in% Seurat::Reductions(object = object)) {
        stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Seurat::Embeddings(object, reduction = reduction)
    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(embedding))
    }
    dims_avail <- seq_len(ncol(embedding))
    if (!all(dims.use %in% dims_avail)) {
        stop("trying to use more dimensions than computed. Rereun dimension reduction\n         with more dimensions or run Harmony with fewer dimensions")
    }
    if (length(dims.use) == 1) {
        stop("only specified one dimension in dims.use")
    }
    metavars_df <- Seurat::FetchData(object, group.by.vars, cells = Seurat::Cells(x = object[[reduction]]))
    harmonyEmbed <- HarmonyMatrix(embedding, metavars_df, group.by.vars,
        FALSE, 0, theta, lambda, sigma, nclust, tau, block.size,
        max.iter.harmony, max.iter.cluster, epsilon.cluster,
        epsilon.harmony, plot_convergence, FALSE, verbose, reference_values)
    # reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
    reduction.key <- "harmony_" # Fixed code line.
    rownames(harmonyEmbed) <- rownames(embedding)
    colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))
    object[[reduction.save]] <- Seurat::CreateDimReducObject(embeddings = harmonyEmbed,
        stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
        assay = Seurat::DefaultAssay(object = object[[reduction]]),
        key = reduction.key)
    if (project.dim) {
        object <- Seurat::ProjectDim(object, reduction = reduction.save,
            overwrite = TRUE, verbose = FALSE)
    }
    return(object)
}


init.dim.reduction <- function(){
  # 1. PCA ---------------------------------------------------
  cat('1. PCA ---------------------------------------------------\n')
  PCA.reports.path <- paste0(dim.reduction.path, '/PCA_reports')
  create.dir(PCA.reports.path, 'PCA reports')
  # Applying PCA
  seurat.obj <- RunPCA(object=seurat.obj, features=VariableFeatures(seurat.obj), npcs=pcs.to.comp, nfeatures.print=15)
  # Saving PCA results.
  tmp.file.name <- get.file.name(gen.path=PCA.reports.path, file.desc='PCALoadingAndPC1vsPC2', extension='pdf')
  pdf(file=tmp.file.name, width=15)
  plot.1 <- VizDimLoadings(seurat.obj, dims=1:2, reduction="pca")
  plot.2 <- DimPlot(seurat.obj, reduction = "pca") + theme(legend.position='element_blank')
  print(CombinePlots(plots=list(plot.1, plot.2)))
  dev.off()
  cat('PC1 and PC2 loading and general info saved to PCA reports path.\n\n')

  # 2. Picking PCs ------------------------------------------------
  cat('2. Picking PCs ------------------------------------------------\n')
  # Save general plots thought to aid in PCs picking decision.
  cat('Saving plots useful for PCs picking decision...\n')
  #   1. Dimensions Heatmap of PC values for the variable features among different cells.
  tmp.file.name <- get.file.name(gen.path=PCA.reports.path, file.desc='DimHeatMap', extension='pdf')
  pdf(file=tmp.file.name)
  for(i in 1:nos.PCs[1]){
  	DimHeatmap(object=seurat.obj, dims=i, cells=500, balanced=TRUE)
  }
  dev.off()
  #   2. Jack Straw Plot
  #seurat.obj <- JackStraw(seurat.obj, num.replicate=100)
  #seurat.obj <- ScoreJackStraw(seurat.obj, dims=1:20)#JackStrawPlot(seurat.obj, dims = 1:20)
  #   3. Elbow Plot
  tmp.file.name <- get.file.name(gen.path=PCA.reports.path, file.desc='ElbowPlot', extension='pdf')
  pdf(file=tmp.file.name, width=10)
  ElbowPlot(seurat.obj, ndims=pcs.to.comp)
  dev.off()
  #   4. Explaining features gain
  # 50 most important features according to weight per PC.
  exp.feats.per.pc <- lapply(X=1:nos.PCs[1], FUN=function(tmp.dim){
    feat.weights <- abs(seurat.obj@reductions$pca@feature.loadings[, tmp.dim])
    feat.weights <- sort(x=feat.weights, decreasing=TRUE)
    tmp.output <- names(feat.weights[1:nos.PCs[1]])
    if(feature.id=='ensembl'){tmp.output <- translate.ids(ids=tmp.output); tmp.output <- names(tmp.output)}
    return(tmp.output)
  })
  # Identify only new features per PC.
  for(tmp.dim in 2:nos.PCs[1]){
    # Take explaining features from this dimension and the last one.
    exp.feats.on <- exp.feats.per.pc[[tmp.dim]]
    exp.feats.prev <- unlist(exp.feats.per.pc[1:(tmp.dim-1)])
    # Remove explaining features already seen in previous PCs.
    exp.feats.on <- exp.feats.on[!exp.feats.on %in% exp.feats.prev]
    exp.feats.per.pc[[tmp.dim]] <- exp.feats.on
  }
  # Sum up.
  exp.feats.gain <- lapply(X=exp.feats.per.pc, FUN=function(exp.feats){
    feats.summ <- if(length(exp.feats)==0) 'None' else if(length(exp.feats)==1) exp.feats else paste0(paste0(exp.feats[1:(length(exp.feats)-1)], collapse=', '), ' and ', exp.feats[length(exp.feats)])
    tmp.output <- data.frame(new.genes.no=length(exp.feats), new.genes.names=feats.summ)
    return(tmp.output)
  })
  exp.feats.gain <- Reduce(x=exp.feats.gain, f=rbind)
  rownames(exp.feats.gain) <- paste0('PC', 1:nos.PCs[1])
  # Output.
  tmp.file.name <- get.file.name(gen.path=PCA.reports.path, file.desc='ExplainingFeatsGain', extension='csv')
  write.csv(file=tmp.file.name, x=exp.feats.gain, quote=TRUE)
  cat('DimHeatMap, Jack Straw and elbow plots saved! You should check them all in order to make a decision about how many PCs would be suitable for downstream analyses!\n\n')
  cat('For this running,', paste0(nos.PCs, collapse=", "), 'PCs will be considered for downstream analyses. If multiple values of number of PCs are required, downstream analyses will be run separately for all of them.\n\n')
  if(!is.null(do.harmony)){
    # 3. Harmony -----------------------------------------------
    cat('3. Harmony -----------------------------------------------\n')
    cat(paste0('Tags to be considered for batch-effect removal:\n', paste0(do.harmony, collapse='\n'), '\n'))
    new.lvls <- as.character(sort(unique(seurat.obj@meta.data[, do.harmony])))
    seurat.obj@meta.data[, 'harmony.var'] <- factor(x=as.character(seurat.obj@meta.data[, do.harmony]), levels=new.lvls)
    seurat.obj <- RunHarmony(object=seurat.obj, group.by.vars='harmony.var')
  }
  return(seurat.obj)
}
