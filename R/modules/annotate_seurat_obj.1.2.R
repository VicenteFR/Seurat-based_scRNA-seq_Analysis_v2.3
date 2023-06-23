############    ----------   SEURAT ANALYSIS    ---------    ############
############    ---   Get a seurat object annotated   ---    ############

# Version: 1
# Version updates:
#   First stable version.
# Subversion: 2
# Subversion updates:
# ---> Subversion 2
# Sub-subversion updates:
#   We allow for the user to decide whether it's ok having more batches defined in the aggr table provided than the actual number of batches defined for the seurat object.
#   Minor bug fix: When hashtag data column is provided in the aggr table but no valid file is indicated (i.e., no file with value different than NA), the program deals with this issue by setting all of the values of column 'hashtag.tag' as NA.

# Dependencies:
#   Seurat and stringr
# Functions dependencies:
#   None.
# Description:
#    Given an aggregation table with annotations regarding each library (i.e., this is intended for seurat object created based on aggregations), this function will provide the appropriate annotations for each cell in the metadata section. Possible annotations are as follows:
#   HTO tags:
#     Description: Hashtag class per cell based on a file for each independent library, as reported by any method applied over hashtag data.
#     Value in aggr table: Absolute path to file describing class per cell.
#   HTO tags:
#     Description: Hashtag class per cell based on a file for each independent library.
#     Value in aggr table: Absolute path to file describing class per cell.
#   Chromium batch tag:
#     Description: Library ID.
#     Value in aggr table: Library ID.
#   Lane tag:
#     Description: Other info conveyed for the library ID. Could be a single value or different tags separated by semicolon.
#     Value in aggr table: Tag value or values separated by semicolon.
#   Chromium batch tag:
#     Description: Library ID.
#     Value in aggr table: Library ID.
#   Sequencing batch tag:
#     Description: Sequencing batch.
#     Value in aggr table: Sequencing batch.
#   Overlay tag:
#     Description: Value conveyed based on the combination of both tags, chromium batch and lane tag.
#     Value in aggr table: None. Added in case the base tags are provided.

# Arguments ------------------------->
# All argument values are taken from the main program.
#   seurat.obj, norm.method (which should be LogNormalize), FVFs.method, mean.cutoff, prop.cutoff, feats.for.dsa.
# Value ----------------------------->
# List containing:
#   Idx. 1) Processed seurat object.
#   Idx. 2) Chracter vector with final variable features.

annotate.seurat.obj <- function(seurat.obj, aggr.table, lane.tag.id, chrom.batch.tag.id, seq.batch.tag.id, hto.tag.id, overlay.tag.id, mult.ann.per.lane, exceded.batches=FALSE){
  ### ---------------------- Data preprocessing ----------------------- ###
  # ---> Samples suffix per cell.
  # Get sample-specific suffix of cell IDs.
  cell.sffxs <- as.integer(str_extract(string=Cells(seurat.obj), pattern='\\d+$'))
  # ---> Batches post-aggregation.
  # Find amount of 10x batches aggregated based on cell barcode suffixes and make sure it's the same amount of batch definitions provided.
  tenx.batches <- as.integer(unique(str_extract(string=Cells(seurat.obj), pattern='\\d+$')))
  if(length(tenx.batches)!=nrow(aggr.table) & !exceded.batches) stop('Not same amount of batches identified in Seurat object as the ones provided in table.')
  # ---> ClassificationPerCell files from aggregation table.
  # Just in case there's HTO data provided.
  if('hto.tag' %in% colnames(aggr.table)){
    hto.tag.vals.list <- lapply(X=aggr.table$hto.tag, FUN=function(hto.tag.file){
      if(is.na(hto.tag.file)){
        hto.tag.vals <- NULL
      }else{
        hto.tag.vals <- read.csv(file=hto.tag.file, stringsAsFactors=FALSE)
        colnames(hto.tag.vals) <- c('barcode', 'hashtag')
      }
      return(hto.tag.vals)
    })
    # Names designated as batch order.
    names(hto.tag.vals.list) <- as.character(1:length(hto.tag.vals.list))
    # Remove NULL values (batches with no HTO info.)
    hto.tag.vals.list <- hto.tag.vals.list[!sapply(X=hto.tag.vals.list, FUN=is.null)]
    if(length(hto.tag.vals.list)==0){
      hto.tag.vals <- NULL
    }else{
      # Unlist sets of HTO tag values and assing identifier according to batch. Then, names will be as in seurat object post-aggregation, composed by barcode and batch.
      # First, get complete cell IDs.
      cell.ids <- unlist(lapply(X=names(hto.tag.vals.list), FUN=function(batch.sfx){
        cell.id <- str_extract(string=hto.tag.vals.list[[batch.sfx]]$barcode, pattern='[ACTG]+')
        cell.id <- paste0(cell.id, '-',batch.sfx)
        return(cell.id)
      }))
      # Then get hashtag vaues unlisted and assign names.
      hto.tag.vals <- unlist(lapply(X=hto.tag.vals.list, FUN=function(hto.tag.vals) return(hto.tag.vals$hashtag)))
      names(hto.tag.vals) <- cell.ids
    }
  }
  ### ---------------------------- Lane tag --------------------------- ###
  if('lane.tag' %in% colnames(aggr.table)){
    # If there are multiple annotations.
    if(mult.ann.per.lane){
      # Tag values.
      lane.tags <- as.data.frame(str_split(string=aggr.table$lane.tag, pattern=';', simplify=TRUE))
      # Correct for NA vals.
      lane.tags[lane.tags=='NA'] <- NA
      # Tag IDs.
      lane.tag.id <- str_split(string=lane.tag.id, pattern=';', simplify=TRUE)[1, ]
      # Add one if it is not properly defined (best guess).
      if(length(lane.tag.id)!=ncol(lane.tags)){
        lane.tag.id <- paste(lane.tag.id[1], 1:ncol(lane.tags), sep='.')
        warning('No appropriate lane tag IDs were defined according to the number of values in aggr table. Next, number of values and number of IDs: ', paste(ncol(lane.tags), length(lane.tag.id), sep='/'))
      }
      # Then, set format.
      colnames(lane.tags) <- lane.tag.id
      lane.tag.vals <- lapply(X=1:nrow(lane.tags), FUN=function(tag.val) return(lane.tags[tag.val, ]))
      # Find lane tag value per cell.
      lane.val.per.cell <- rbindlist(lane.tag.vals[cell.sffxs])
      # Add lane tag value per cell.
      seurat.obj@meta.data <- cbind(seurat.obj@meta.data, lane.val.per.cell)
    }else{
    # If there's a single value.
      lane.tag.vals <- lapply(X=aggr.table$lane.tag, FUN=function(tag.val) return(tag.val))
      # Find lane tag value per cell.
      lane.val.per.cell <- unlist(lane.tag.vals[cell.sffxs])
      # Add lane tag value per cell.
      seurat.obj[[lane.tag.id]] <- lane.val.per.cell
    }
  }
  ### ---------------------------- HTO tag ---------------------------- ###
  if('hto.tag' %in% colnames(aggr.table)){
    if(is.null(hto.tag.vals)){
      seurat.obj[[hto.tag.id]] <- NA
    }else{
      classes.df <- data.frame(barcode.id=names(hto.tag.vals), class=hto.tag.vals, stringsAsFactors=FALSE)
      cells.df <- data.frame(barcode.id=Cells(seurat.obj), stringsAsFactors=FALSE)
      hto.tag.vals.per.cell <- merge(x=cells.df, y=classes.df, by='barcode.id', all.y=FALSE, all.x=TRUE, sort=FALSE)
      rownames(hto.tag.vals.per.cell) <- hto.tag.vals.per.cell$barcode.id
      hto.tag.vals.per.cell <- hto.tag.vals.per.cell[cells.df$barcode.id, ]
      # Add HTO tag value per cell.
      seurat.obj[[hto.tag.id]] <- hto.tag.vals.per.cell$class
    }
  }
  ### ----------------------- Chromium batch tag ---------------------- ###
  if('chrom.batch.tag' %in% colnames(aggr.table)){
    # Find batch tag value per cell.
    chrom.batch.tag.vals <- lapply(X=aggr.table$chrom.batch.tag, FUN=function(tag.val) return(tag.val))
    chrom.batch.val.per.cell <- unlist(chrom.batch.tag.vals[cell.sffxs])
    # Add batch tag value per cell.
    seurat.obj[[chrom.batch.tag.id]] <- chrom.batch.val.per.cell
  }
  ### ---------------------- Sequencing batch tag --------------------- ###
  if('seq.batch.tag' %in% colnames(aggr.table)){
    # Find batch tag value per cell.
    seq.batch.tag.vals <- lapply(X=aggr.table$seq.batch.tag, FUN=function(tag.val) return(tag.val))
    seq.batch.val.per.cell <- unlist(seq.batch.tag.vals[cell.sffxs])
    # Add batch tag value per cell.
    seurat.obj[[seq.batch.tag.id]] <- seq.batch.val.per.cell
  }
  ### ------------------------- Overlay tag --------------------------- ###
  # Tag provided by the combination of both, HTO and chromium batch tags.
  if(all(c('chrom.batch.tag', 'hto.tag') %in% colnames(aggr.table))){
    overlay.vals.per.cell <- paste(seurat.obj@meta.data[, hto.tag.id], seurat.obj@meta.data[, chrom.batch.tag.id], sep='-')
    # Find NA values from any independent column and report the same for overlay.
    overlay.vals.per.cell[is.na(seurat.obj@meta.data[, hto.tag.id]) | is.na(seurat.obj@meta.data[, chrom.batch.tag.id])] <- NA
    # Add overlay tag per cell.
    seurat.obj[[overlay.tag.id]] <- overlay.vals.per.cell
  }
  # ---> Return
  return(seurat.obj)
}
