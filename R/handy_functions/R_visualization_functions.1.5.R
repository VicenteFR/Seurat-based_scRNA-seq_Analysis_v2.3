### ----------------------------------------------------------------- ###
### --------------- Functions for data visualization ---------------- ###
### ----------------------------------------------------------------- ###

# By: Vicente Fajardo.

# 1 -------------------------------------------------------------------->
# Name: Violin plot.
# Updates:
# @     Compared to the function in version 1.3, this one:
# 1. Can take log2(CPM + 1) values for the feature of interest apart of the rest already available.
# 2. Allows for the filtering of cells based on multiple tags and, tehrefore, based on values of multiple tags.
# @     Compared to the function in version 1.4, this one:
# 1. Allows for the merging of the groups of interest into a unique group which is argument 'groups.label'. Additionally, this argument also controls whether such groups should merged (when an appropriate label is provided) or not (when indicated as NULL).
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If required (groups.of.int not NULL), keep the original identity of the groups indicated only, setting the others' identity to 'rest'
# 3. If indicated (filter.tags set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of one or multiple tags (filter.tags). See details below as to how to indicate filtring based on multiple tags.
# 4. Set a color scale for each group (either the final identities -"color"- or a stat from mean, median or percentage of expressing cells).
# 5. Generate a ggplot, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# feature - Feature for which to output density. Should be defined either in the metadata (slot could either be NULL for straight and unambiguous picking or simply not be defined in any data slot) or as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. NULL indicates 'data', 'cpm' indicates log2(CPM + 1) relative counts, while the rest are the usual ones: counts, data or scale.data. No other value is allowed. Default: data.
# groups.tag - Tag defined in the seurat onject metadata whose groups will be the basis for the x axis of the product. No default.
# na.rm - Logical, indicates whether NA values should be discarded for the tag listing groups. Default: TRUE.
# groups.of.int - Character listing the set of tag groups whose identity should be kept as originally defined, otherwise they'll be set to 'rest'. This step is skipped is it's set to NULL. Default: NULL.
# groups.label - Character. When a non-NULL values provided, indicates that the groups of interest provided through argument 'groups.of.int' should be in turn merged into a single group and the value given through this argument will be used as a label for such a single group. Otherwise, a NULL value indicates that the groups of interest shouldn't be merged. Default: NULL.
# filter.tags - Chracter vector, which should list the tag or tags in the seurat object metadata whose groups will be used to filter out cells. When filtering should be applied based on multiple tags, their groups to be considered to keep or exclude cells must be passed as chracter vectors embedded within a list; then, this character vector must have the same length as that of the list provided through the argument 'groups.to.filter'. No cell filtering is applied when it's set to NULL. Default: NULL.
# groups.to.filter - List of vectors for the set of tag groups that'll be kept ot removed according to the argument 'keep' value. Must have the smae length as the vector provided in 'filter.tags'. This step is skipped is it's set to NULL. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
# color - Value for color scale of violin, either the final identities -"color"- or a stat from mean -"mean"-, median -"median"-, burst frequency (percentage of expressing cells) -"burst.frequency"- or burst size (mean of expressing cells, where a cell is considered as such when expression > threshold, (see argument size.thold)) -"burst.size"-. Default: 'color'.
# size.thold - Threshold over which cells are going to be considered as 'expressing or not'. Ignored when color different than 'burst.size'. Default: 0
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# adjust.val - Value for argument 'adjust' of geom_violin. Default: 0.8. See details at https://ggplot2.tidyverse.org/reference/geom_violin.html
# trim.val - Value for argument 'trim' of geom_violin. Default: FALSE. See details at https://ggplot2.tidyverse.org/reference/geom_violin.html
# plot.limits - Value for argument adjust of scale_fill_gradientn. Only valid if argument color's value is different than 'color'. Default: NULL. See details at: https://ggplot2.tidyverse.org/reference/scale_gradient.html
# add.points - Logical, indicates whether jittered points should be depicted or not in the violin plot. If TRUE, only 20% of the data points (picked randomly) will be depicted. Default: FALSE.
# this.color.scale - Character. Color scale or color values to use.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# Function:

vln.plot <- function(seurat.obj, feature, slot='data', groups.tag, na.rm=TRUE, groups.of.int=NULL, groups.label=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, color='color', size.thold=0, file.name=NULL, adjust.val=0.8, trim.val=FALSE, plot.limits=NULL, add.points=FALSE, this.color.scale=NULL){
  # ---> Arguments and checks.
  # Argument defining violin plots' fill color.
  if(!color %in% c('median', 'mean', 'burst.frequency', 'burst.size', 'color')) stop(paste0('Valid value for "color" must be provided among:\n\t', paste0(c('median', 'mean', 'burst.frequency', 'burst.size', 'color'), collapse=', '), '\n'))
  if(is.null(this.color.scale) & color!='color') this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # provided by Ciro.
  # Define CPMs-related arguments.
  if(!is.null(slot)){
    if(slot == 'cpm'){
      slot <- 'counts'
      cpm <- TRUE
    }else{
      cpm <- FALSE
    }
  }else{
    cpm <- FALSE
  }
  # Filter tags and groups to filter dimensions must agree with each other.
  if(!is.null(filter.tags)){
    if(length(filter.tags)>1){
      if(!is.list(groups.to.filter) | length(filter.tags)!=length(groups.to.filter)) stop(paste0('When multiple tags are provided for cell filtering, the argument groups.to.filter must be a list with the same dimension length.\n'))
    }else{
      if(!is.list(groups.to.filter)){
        if(!is.vector(groups.to.filter)) stop('Argument "groups.to.filter" must be either a list (with the same length as that for the argument "filter.tags") or a vector.\n') else groups.to.filter <- list(groups.to.filter)
      }
    }
    names(groups.to.filter) <- paste0('filter.tag.', 1:length(filter.tags))
  }
  # ---> Fetch data.
  # Look for feature.
  # First, look for it in the metadata if slot is set to NULL.
  if(is.null(slot)){
    if(!(feature %in% colnames(seurat.obj@meta.data))) stop('Slot set to NULL, but feature required was not found in seurat object\'s metadata. Were you trying to look for it any slot?\n')
    feature.vals <- data.table(feature=seurat.obj@meta.data[, c(feature)])
  }else{
    # If set to a value different than NULL, make sure it's a valid slot and try to find the feature there.
    if(!(slot %in% c('counts', 'data', 'scale.data'))) stop('Not valid slot. Should be anyone among "counts", "data", "scale.data".\n')
    # Proceed if all good, though process branches based on CPM possible requirement.
    if(cpm){
      if(!feature %in% row.names(seurat.obj)) stop('Feature could not be found in the count slor, which is required to calculate CPMs as requested. May want to try a different slot?\n')
      feature.vals <- RelativeCounts(data=GetAssayData(object=seurat.obj, slot=slot, assay='RNA'), scale=1e6, verbose=FALSE)
      feature.vals <- log2(feature.vals[feature, ] + 1)
      feature.vals <- as.data.table(feature.vals)
      colnames(feature.vals) <- 'feature'
    }else{
      feature.vals <- try(expr=as.data.table(FetchData(object=seurat.obj, vars=feature, slot=slot)), silent=TRUE)
      # Only if it was not possible to find the feature in the data slot appointed, try again in the metadata.
      if(any(class(feature.vals)=='try-error')){
        if(!(feature %in% colnames(seurat.obj@meta.data))) stop('Feature found nowhere. Please enter a valid feature. Alternatively, package \'data.table\' may not be imported.\n')
        feature.vals <- data.table(feature=seurat.obj@meta.data[, c(feature)])
      }else{
        # If not, set the appropriate column name for format consistency.
        colnames(feature.vals) <- 'feature'
      }
    }
  }
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!(groups.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for groups of interest not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tags)){
    tmp.assess <- filter.tags %in% colnames(seurat.obj@meta.data)
    if(!all(tmp.assess)){
      tmp.assess <- filter.tags[!tmp.assess]
      stop(paste0('Tags for filtering set to a value different than NULL, but some of them not defined in seurat object\'s metadata, as follows.\n', paste0(tmp.assess, collapse='\n'), '\n'))
    }
    tmp.data <- data.table(seurat.obj@meta.data[, filter.tags])
    colnames(tmp.data) <- paste0('filter.tag.', 1:ncol(tmp.data))
    tmp.data <- cbind(data.table(groups.tag=seurat.obj@meta.data[, groups.tag]), tmp.data)
  }else{
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag])
  }
  # Join altogether.
  tmp.data <- cbind(tmp.data, feature.vals)
  # ---> Groups tag.
  # Filter out non-desired cells according to groups.tag
  if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
  # Set rest group for the groups tag if necessary.
  if(!is.null(groups.of.int)){
    tmp.data[!groups.tag %in% groups.of.int, groups.tag:='Rest']
  # Similarly, merge groups of interest if their label is provided.
    if(!is.null(groups.label)){
      tmp.data[groups.tag %in% groups.of.int, groups.tag:=as.character(groups.label)]
      # Set factors leaving 'Rest' label right at the end of the levels' order.
      new.lvls <- c(as.character(groups.label), 'Rest')
      tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=new.lvls)
    }
  }
  # ---> Filter out non-desired cells.
  if(!is.null(filter.tags)){
    # Filter all cells with not available data regarding all column-to-filter information.
    if(na.rm){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        tmp.data <- tmp.data[!is.na(get(tmp.col))]
      }
    }
    # Then, when requested, keep or remove input values.
    if(!is.null(groups.to.filter)){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        if(keep) tmp.data <- tmp.data[get(tmp.col) %in% groups.to.filter[[tmp.col]]] else tmp.data <- tmp.data[!get(tmp.col) %in% groups.to.filter[[tmp.col]]]
      }
    }
  }
  # ---> Apply feature thereshold.
  if(!is.null(feature.thold)){
    if(!is.numeric(feature.thold)) stop('Not valid feature threshold provided. Must be a numeric value.\n')
    tmp.data <- tmp.data[feature > feature.thold]
  }
  # ---> Color scale.
  if(color=='color'){ # Color per group.
    tmp.data[, color:=groups.tag]
  }else{ # Stats
    switch(EXPR=color,
      "median"={
        tmp.dt <- tmp.data[, .(color=median(feature, na.rm=TRUE)), by=groups.tag]
      },
      "mean"={
        tmp.dt <- tmp.data[, .(color=mean(feature, na.rm=TRUE)), by=groups.tag]
      },
      "burst.frequency"={ # Percent of expresing cells.
        tmp.dt <- tmp.data[, .(color=.SD[feature > 0, .N]*100/.N), by=groups.tag]
      },
      "burst.size"={ # Mean of expressing cells.
        tmp.dt <- tmp.data[feature>size.thold, .(color=mean(feature, na.rm=TRUE)), by=groups.tag]
      }
    )
    tmp.data <- merge(x=tmp.data, y=tmp.dt, by='groups.tag')
  }
  # ---> Points
  # 20% of the points per group are going to be added.
  if(add.points){
    rows.to.keep <- sample(x=1:nrow(tmp.data), size=nrow(tmp.data)*0.05)
    tmp.data.sample <- tmp.data[rows.to.keep, ]
  }
  # ---> Violin plot.
  # Get ggplot.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=groups.tag, y=feature, fill=color)) + geom_violin(alpha=0.7, trim=trim.val, adjust=adjust.val) + geom_boxplot(width=0.05, alpha=0.7, outlier.shape=NA) + labs(x=paste0('Groups (', groups.tag, ')'), y=feature)
  # Adjust color scale.
  if(!is.null(this.color.scale) & color=='color') tmp.ggplot <- tmp.ggplot + scale_fill_manual(values=this.color.scale)
  if(color!='color') tmp.ggplot <- tmp.ggplot + scale_fill_gradientn(colors=this.color.scale, limits=plot.limits)
  # Add points if necessary.
  if(add.points) tmp.ggplot <- tmp.ggplot + geom_jitter(data=tmp.data.sample, color='black', size=0.2, alpha=0.6, width=0.2)
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}


# 2 -------------------------------------------------------------------->
# Name: Dot plot.
# Of note: This function's idea was based on Seurat's output dot plot, but it was developed independently up to the previous version (see past function's file version). The size scale function implementation was copied from Seurat's function.
# Updates:
# @     Compared to the function in version 1.4, this one:
# 1. If required, applies a row-wise normalization by z-scoring.
# 2. Can set boundaries for both scales, color and size, through arguments col.min/col.max and size.min/size.max.
# 3. If required, sets a hot and cold color scale taking into account the minimum and maximum color scale values by centering at zero. This won't work if there's only negative or only positive values. (This needs further revising).
# 4. Sets a scale by taking into consideration either the size or the radius of the points in the plot. Original seurat's function does apply a radius scaling and it seems to be the best option for this kind of plots.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If required (groups.of.int not NULL), keep the original identity of the groups indicated only, setting the others' identity to 'rest'
# 3. If indicated (filter.tag set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of a tag (filter.tag).
# 4. Calculate expression mean and burst frequency for each gene across groups.
# 5. Normalize the mean expression values across the groups of interest (i.e., per gene or in a row-wise manner).
# 6. Set color and size scales by setting boundaries and plugging the appropriate arguments.
# 7. Generate a dot plot (based on seurat's style) as a ggplot, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
#     @ About data.
#         Data specificities.
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# features - Features for which to output expression values' maetrics. All should be defined as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. Valid values are the usual ones: counts, data or scale.data. Default: data.
# do.norm - Logical, indicates whether a row-wise normalization should be applied through z-scoring. Default: FALSE.
# ensembl - Logical, defines whether input features' vector (defined through argument 'fetaures') lists ENSEMBL IDs (then expecting the gene names to be defined as the vector's names). Default: FALSE.
#         Data filtering and grouping.
# groups.tag - Tag defined in the seurat object metadata whose groups will be the basis for the x axis of the product. No default.
# groups.of.int - Character listing the set of tag groups whose identity should be kept as originally defined, otherwise they'll be set to 'rest'. This step is skipped is it's set to NULL. Default: NULL.
# groups.order - Order to show the groups in the x axis. When an order is provided, all groups defined in the object (post filtering, if any) are expected to be in the order vector and the other way around. If none provided (i.e., NULL value), the function will attempt to sort them based on their names (in an alphanumeric order). Default: NULL
# filter.tag - Tag defined in the seurat object metadata whose groups will be used to filter out cells. No cell filtering is applied when this argument's set to NULL. Default: NULL.
# groups.to.filter - Character listing the set of tag groups that'll be kept ot removed according to the argument 'keep' value. This step is skipped is it's set to NULL. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Default: TRUE
# na.rm - Logical, indicates whether NA values should be discarded for the tag listing groups. Default: TRUE.
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
#     @ About plot and output.
#         Color scale.
# this.color.scale - Character, indicates color scale to use to represent average expression. If set to NULL, a color scale that's broadly used in Vijay's lab will be used instead. If a length-1 vector that equals 'hot.and.cold', 'hot.n.cold' or 'cold.and.hot' is input, a hot and cold color scale is infered by taking into account the range of average expression values seen in the data (post normalization if required).
# col.min - Numeric, lower boundary for color scale. Must be smaller than col.max (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
# col.max - Numeric, upper boundary for color scale. Must be larger than col.min (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
#         Color scale.
# scale.by - Character, indicates whether a radius scaling or a standard size scaling should be applied. Accepted values are 'radius' and 'size', respectively. Default: 'radius'.
# dot.scale - Character, second value in the argument 'range' of 'scale_radius' or 'scale_size'. Default: 6.
# size.min - Numeric, lower boundary for size scale. Must be smaller than size.max (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
# size.max - Numeric, upper boundary for size scale. Must be larger than col.min (if provided) to input an appropriate color range. If NULL, no boundary is set. Default: NULL.
#         Output.
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# Function:

dot.plot <- function(
  # @   About data.
  # Data specificities.
  seurat.obj, features, slot='data', do.norm=FALSE, ensembl=FALSE,
  # Data filtering and grouping.
  groups.tag, groups.order=NULL, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=TRUE, na.rm=TRUE, feature.thold=NULL,
  # @   About plot and output.
  # Color scale
  this.color.scale=NULL, col.min=NULL, col.max=NULL,
  # Size scale.
  scale.by='radius', dot.scale=6, size.min=NA, size.max=NA,
  # File.
  file.name=NULL
){
  # ---> Arguments.
  # Set color scale.
  if(is.null(this.color.scale)){
    this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
  }else{
    if(this.color.scale=='hot.and.cold' | this.color.scale=='cold.and.hot' | this.color.scale=='cold.n.hot') this.color.scale <- NA
  }
  # Set size scale.
  scale.func <- switch(
    EXPR=scale.by,
    size=scale_size,
    radius=scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  # ---> Fetch data.
  # Look for all features. Make sure it's been saved to a valid slot and try to find the features there.
  if(!(slot %in% c('counts', 'data', 'scale.data'))) stop('Not valid slot. Should be anyone among "counts", "data", "scale.data".\n')
  feature.vals <- try(expr=as.data.table(FetchData(object=seurat.obj, vars=features, slot=slot)), silent=TRUE)
  # Kill if it was not possible to find the feature in the data slot appointed.
  if(any(class(feature.vals)=='try-error')) stop('Feature found nowhere. Please enter a valid feature.\n')
  colnames(feature.vals) <- if(ensembl) names(features) else features
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!(groups.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for groups of interest not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tag)) if(!(filter.tag %in% colnames(seurat.obj@meta.data))) stop('Tag for filtering set to a value different than NULL, but not defined in seurat object\'s metadata.\n')
  if(!is.null(filter.tag)){
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag], filter.tag=seurat.obj@meta.data[, filter.tag])
  }else{
    tmp.data <- data.table(groups.tag=seurat.obj@meta.data[, groups.tag])
  }
  # Join altogether.
  tmp.data <- cbind(tmp.data, feature.vals)
  # ---> Groups tag.
  # Filter out non-desired cells according to groups.tag
  if(na.rm) tmp.data <- tmp.data[!is.na(groups.tag)]
  # Set rest group for the groups tag if necessary.
  if(!is.null(groups.of.int)){
    tmp.data[!groups.tag %in% groups.of.int, groups.tag:='rest']
  }
  # ---> Filter out non-desired cells.
  # According to filter tag.
  if(!is.null(filter.tag) & !is.null(groups.to.filter)){
    if(keep) tmp.data <- tmp.data[filter.tag %in% groups.to.filter] else tmp.data <- tmp.data[!filter.tag %in% groups.to.filter]
    tmp.data[, filter.tag:=NULL]
  }
  # ---> Tidy data up!
  tmp.data <- as.data.table(gather(data=tmp.data, key='feature', value='exp.val', -groups.tag))
  # ---> Apply feature thereshold.
  if(!is.null(feature.thold)){
    if(!is.numeric(feature.thold)) stop('Not valid feature threshold provided. Must be a numeric value.\n')
    tmp.data <- tmp.data[exp.val > feature.thold]
  }
  # ---> Values for color and size scales.
  tmp.data <- tmp.data[, .(exp.mean=mean(exp.val, na.rm=TRUE), burst.freq=.SD[exp.val > 0, .N]*100/.N), by=.(groups.tag, feature)]
  # ---> Row-wise normalization.
  # Apply only if required.
  if(do.norm){
    tmp.data.2 <- as.data.frame(spread(data=tmp.data[, .(groups.tag, feature, exp.mean)], key=groups.tag, value=exp.mean))
    row.names(tmp.data.2) <- tmp.data.2$feature; tmp.data.2$feature <- NULL
    tmp.data.2 <- as.data.frame(scale(x=tmp.data.2))
    tmp.data.2$feature <- row.names(tmp.data.2); tmp.data.2 <- as.data.table(tmp.data.2)
    tmp.data.2 <- gather(data=tmp.data.2, key='groups.tag', value='exp.mean', -feature)
    tmp.data <- merge(x=tmp.data[, .(groups.tag, feature, burst.freq)], y=tmp.data.2, by=c('feature', 'groups.tag'))
  }
  # ---> Set boundaries for scales and get final color definitions.
  # Color scale.
  if(!is.null(col.min) & !is.null(col.max)) if(col.min > col.max) stop('Boundaries for color scale were provided, but the minimum value to be set is larger than the maximum one. Provide a valid range if any is required.\n')
  if(!is.null(col.min)) tmp.data[exp.mean < col.min, exp.mean:=col.min]
  if(!is.null(col.max)) tmp.data[exp.mean > col.max, exp.mean:=col.max]
  if(is.na(this.color.scale)){
    scale.lims <- tmp.data[, c(max=max(exp.mean), min=min(exp.mean))]
    is.max <- names(scale.lims)[abs(scale.lims)==max(abs(scale.lims))]
    max.col.count <- 1001
    min.col.count <- floor(min(abs(scale.lims))/max(abs(scale.lims)) * max.col.count)
    max.col.scale <- heatmaply::cool_warm(max.col.count)
    min.col.scale<- heatmaply::cool_warm(min.col.count)
    if(is.max=='max'){
      this.color.scale <- c(min.col.scale[1:floor(min.col.count/2)], max.col.scale[floor(max.col.count/2):length(max.col.scale)])
    }else{
      this.color.scale <- c(max.col.scale[1:floor(max.col.count/2)], min.col.scale[floor(min.col.count/2):length(min.col.scale)])
    }
  }
  # Size scale.
  if(!is.na(size.min) & !is.na(size.max)) if(size.min > col.max) stop('Boundaries for size scale were provided, but the minimum value to be set is larger than the maximum one. Provide a valid range if any is required.\n')
  if(!is.na(size.min)) tmp.data[burst.freq<size.min, burst.freq:=size.min]
  if(!is.na(size.max)) tmp.data[burst.freq>size.max, burst.freq:=size.max]
  # ---> Set factors.
  # For groups.
  if(is.null(groups.order)){
    new.lvls <- mixedsort(x=unique(as.character(tmp.data$groups.tag)))
    tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=new.lvls)
  }else{
    if(!(all(groups.order %in% unique(as.character(tmp.data$groups.tag))) & all(unique(as.character(tmp.data$groups.tag)) %in% groups.order))) stop('When an order is provided, all groups defined in the object (post filtering, if any) are expected to be in the order vector and the other way around.\n')
    tmp.data$groups.tag <- factor(x=as.character(tmp.data$groups.tag), levels=as.character(groups.order))
  }
  # For features (same order as the one provided)
  new.lvls <- if(ensembl) names(features) else features
  tmp.data$feature <- factor(x=tmp.data$feature, levels=rev(new.lvls))
  # ---> Dot plot.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=groups.tag, y=feature, size=burst.freq, col=exp.mean)) +
    geom_point() +
    scale_color_gradientn(colors=this.color.scale) +
    scale.func(range=c(0, dot.scale), limits=c(size.min, size.max)) +
    labs(x='Groups', y='Gene', col='Average expression', size='Percent expressed') +
    theme(panel.background=element_blank(), panel.border=element_blank(), axis.line=element_line())
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}


# 3 -------------------------------------------------------------------->
# Name: Heatmap with DEA results sorted according to sets' values and sizes.
# Description:
# Given a table (provided either as a data frame or a data table) listing the metrics ussually gotten when applying multiple pairwise comparisons through DEA (mainly, adjusted p-value, LFC and mean values for the groups that were compared), this function will provide a tidied version of such table listing the DEGs in an ordered way that depends on both the reference elements of each comparison (set value) and the size of groups of elements that share a DEG (set size). When a heatmap is plotted based on this order (provided an order of the elements), it must have a scalonated pattern. In order to get such results, the function will:
# 1. Tidy the data up keeping relevan columns and calculate the size set each DEG belongs to.
# 2. Set order, first according to set size, then to LFC and finally to elements.
# Arguments ------------------------->
# dea.results - DEA results table. Data frame or data table, must contain all DEA metrics as well as some identifier for each DEG. No default.
# feat.id - Character, name of the column listing DEG's ID. Default: 'gene.id'
# element - Character, name of the column listing each element value (i.e., the value for which the DEA was applied (either as reference or not)). Make sure this column's name is different than 'element', since the function will not work otherwise. No default.
# vals.order - Character vector containing the order for the element values desired to show on the heatmap. No default.
# p.adj - Character, name of the column listing DEG's adjusted p values. Default: 'p.adj'
# lfc - Character, name of the column listing DEG's LFC. Default: 'lfc'
# mean.tag - Character, pattern to identify columns listing DEA group means. No default.
# element.eq.means - Logical, indicates whether element values are the same values for which a mean is provided. If TRUE, the function will try to get rid of the information taht's repeated along the process. Default: FALSE.
# file.name - To implement. If not NULL, the function should output the heatmap, too.
# Value ----------------------------->
# Datatable listing DEGs main metrics ordered according to both element values and set sizes.
# Function:

dea.heatmap <- function(dea.results, feat.id='gene.id', element, vals.order, p.adj='p.adj', lfc='lfc', mean.tag, element.eq.means=FALSE, file.name=NULL){
  # ---> Preflights.
  if(!is.data.table(dea.results)) dea.results <- as.data.table(dea.results)
  # ---> Define stuff
  # Stats.
  mean.cols <- colnames(dea.results)[str_detect(string=colnames(dea.results), pattern=mean.tag)]
  stat.cols <- c(p.adj, lfc, mean.cols)
  col.vals <- str_replace(string=mean.cols, pattern=mean.tag, replacemen='')
  # ---> Tidy data up!
  stats.vals <- apply(X=as.matrix(dea.results[, ..stat.cols]), MARGIN=1, FUN=paste, collapse=';')
  dea.results[, stats:=stats.vals]
  tmp.data <- dea.results[, .(gene.id=get(feat.id), element=get(element), stats)]
  tmp.data <- spread(data=tmp.data, key='element', value='stats', fill=paste('0;0', paste0(rep(x='0', times=length(mean.cols)), collapse=';'), sep=';'))
  for(element.val in vals.order){
    tmp.data <- separate(data=tmp.data, col=element.val, into=paste(c('p.val', 'lfc', paste0('mean.', col.vals)), element.val, sep='.'), sep=';', convert=TRUE)
  }
  # Discard repeated information when the mean groups (values) are the same as the element values.
  if(element.eq.means){
    mean.cols <- paste('mean', col.vals, col.vals, sep='.')
    if(!all(mean.cols %in% colnames(tmp.data))) stop('It was requested to remove repeated information, though mean column values don\'t seem to be the same as the element values (could be that not all element values have a mean column). Please set element.eq.means to FALSE and see if your results are still favorable. Otherwise, check your input.\n')
    other.cols <- colnames(tmp.data)[!str_detect(string=colnames(tmp.data), pattern='mean.')]
    cols.to.keep <- c(other.cols, mean.cols)
    tmp.data <- tmp.data[, ..cols.to.keep]
    mean.cols <- str_replace(string=mean.cols, pattern=paste0('.', col.vals, '$'), replacement='')
    cols.to.rep <- c(other.cols, mean.cols)
    colnames(tmp.data) <- cols.to.rep
  }
  # Identify for what cell types a gene is a DEG as well as the total amount of cell types it is DEG for.
  lfc.cols <- grep(x=colnames(tmp.data), pattern='lfc', value=TRUE)
  p.cols <- grep(x=colnames(tmp.data), pattern='p.val', value=TRUE)
  rows.to.keep.1 <- apply(X=tmp.data[, ..lfc.cols], MARGIN=1, FUN=function(row) return(row>0.25))
  rows.to.keep.2 <- apply(X=tmp.data[, ..p.cols], MARGIN=1, FUN=function(row) return(row<0.05))
  rows.to.keep.mat <- rows.to.keep.1&rows.to.keep.2
  row.names(rows.to.keep.mat) <- str_replace(string=row.names(rows.to.keep.mat), pattern='lfc\\.', replacement='')
  # While doing so, filter for genes that are significantly differentially expressed for at least one cell type (apply to both objects).
  rows.to.keep <- apply(X=rows.to.keep.mat, MARGIN=2, FUN=any)
  tmp.data <- tmp.data[rows.to.keep]
  rows.to.keep.mat <- rows.to.keep.mat[, rows.to.keep]
  # Add sets sizes and values.
  degs.sets.info <- t(apply(X=rows.to.keep.mat, MARGIN=2, FUN=function(col){
    tmp.vals <- row.names(rows.to.keep.mat)[col]
    tmp.size <- length(tmp.vals)
    tmp.vals <- paste0(sort(tmp.vals), collapse=';')
    return(c(tmp.vals, tmp.size))
  }))
  colnames(degs.sets.info) <- c('set.values', 'set.size')
  tmp.data <- cbind(tmp.data, degs.sets.info)
  # ---> Set order.
  # Get sets' sizes and values to apply iterations.
  set.sizes <- 1:(dea.results[, length(unique(get(element)))])
  set.vals <- dea.results[, unique(get(element))]
  if(!all(sort(set.vals)==sort(vals.order))) stop('Not all element values were provided in argument \'vars.order\'')
  set.vals <- vals.order
  vals.order <- 1:length(set.vals)
  names(vals.order) <- set.vals
  # Get DEGs' order.
  set.degs <- lapply(X=set.sizes, FUN=function(tmp.size){
    # Get DEGs set for the size of a given iteration.
    set.degs <- tmp.data[set.size==tmp.size]
    # Then, identify possible values for the whole size set.
    sets.vals <- str_split(string=set.degs[, set.values], pattern=';', simplify=TRUE)
    # Order size set according to values.
    set.degs <- lapply(X=set.vals, FUN=function(tmp.val){
      # Idenitfy values that have already been taken into account in previos iterations.
      to.exclude <- set.sizes[set.vals == tmp.val]-1
      if(to.exclude>0) to.exclude <- 1:to.exclude
      to.exclude <- set.vals[to.exclude]
      # Get DEGs size subset for the value of a given iteration, disregarding the DEGs belonging to any value to exclude.
      rows.to.keep <- apply(X=sets.vals, MARGIN=1, FUN=function(row) return(tmp.val %in% row & !any(to.exclude %in% row)))
      set.degs <- set.degs[rows.to.keep]
      sets.vals <- sets.vals[rows.to.keep, ]
      # Check we've got info to work with.
      if(all(!rows.to.keep)) return(NA)
      # Identify LFC columns to consider to get a mean LFC value.
      iter.set.vals <- set.degs[, unique(set.values)]
      for(tmp.val.set in iter.set.vals){
        tmp.lfc.cols <- paste0('lfc.', str_split(string=tmp.val.set, pattern=';', simplify=TRUE))
        set.degs[, tmp.avg.lfc:=rowMeans(set.degs[, ..tmp.lfc.cols])]
        set.degs[set.values==tmp.val.set, avg.lfc:=tmp.avg.lfc]
      }
      set.degs[, tmp.avg.lfc:=NULL]
      # Check we've got dimensions to work over.
      if(is.null(dim(sets.vals))) return(set.degs)
      # Set order according to mean LFC while keeping the order for object having set's values, too.
      set.degs[, idx:=1:.N]
      setorderv(x=set.degs, cols='avg.lfc', order=-1)
      sets.vals <- sets.vals[set.degs[, idx], ]
      set.degs[, idx:=NULL]
      # Set order according to set values.
      sets.vals <- as.data.table(t(apply(X=sets.vals, MARGIN=1, FUN=function(row) return(sort(vals.order[row])))))
      row.names(sets.vals) <- as.character(1:nrow(sets.vals))
      colnames(sets.vals) <- as.character(1:ncol(sets.vals)); cols.to.order <- as.character(1:ncol(sets.vals)); sets.vals[, idx:=1:.N]
      for(tmp.col in rev(cols.to.order)) setorderv(x=sets.vals, col=tmp.col, order=1)
      set.degs <- set.degs[sets.vals[, idx]]
      # return
      return(set.degs)
    })
    set.degs <- set.degs[!is.na(set.degs)]
    set.degs <- rbindlist(l=set.degs, use.names=TRUE, idcol=FALSE)
    return(set.degs)
  })
  set.degs <- rbindlist(l=set.degs, use.names=TRUE, idcol=FALSE)
  # ---> Signature depicted by a heatmap.
  # To implement.
  # Return sorted data table.
  return(set.degs)
}


# 4 -------------------------------------------------------------------->
# Name: Get density values.
# Taken from Ciro's clever functions.
# Description:
# Given the x values and the y values (compraising a grid containing 'n' cells -see parameter 'n'-) for a set of 2D points, this function will calculate the density values of such a grid and then will define the 2D density values for each individual point.
# Arguments ------------------------->
# x - x values for the set of points.
# y - y values for the set of points.
# ... - Other arguments taken by function 'kde2d' from the package 'MASS'. Among them, check parameter 'n' which basically defines the number of cells to be contained in the grid (see https://rdrr.io/cran/MASS/man/kde2d.html).
# Value ----------------------------->
# Matrix with 2D density values (columns) for each input point (rows).
# Function:

get.densities <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


# 5 -------------------------------------------------------------------->
# Name: Co-expression plot.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If indicated (filter.tag set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of a tag (filter.tag).
# 3. Define expression-related clasess for each cell by checking whether expression for each feature is larger than a set threshold (express.thold). Possible classes for each cell are: Double positive (DP, expresses both features over the threshold); Double negative (DN, expresses none of the features over the threshold); and single positive for either feature (SP.X or SP.Y when expressing the feature depicted in the x coordinate or the y coordinate of the scatter plot, respectively).
# 4. Calculate the fraction of each expression-related class to depict on the final plot.
# 5. Calculate the 2D kernel density for each cell.
# 6. Generate a ggplot object depicting a scatter plot with density as color, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# feature.x & feature.y - Features for which to depict their expression values in scatter plot (x and y axes, respectively). Both should be defined as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. NULL indicates 'data', 'cpm' indicates log2(CPM + 1) relative counts, while the rest are the usual ones: counts, data or scale.data. No other value is allowed. Default: data.
# filter.tags - Tags defined in the seurat onject metadata whose groups will be used to filter out cells. No cell filtering is applied when it's set to NULL. Default: NULL.
# na.rm - Logical, indicates whether NA values defined the filter.tag column should be discarded for the tag listing groups. Default: TRUE.
# groups.to.filter - Character or list, listing the set of tag groups that'll be kept ot removed according to the argument 'keep' value. This step is skipped is it's set to NULL. When more than one tag are provided for filtering, this argument must be a list with the same dimensions as that argument. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Recycled for all filtering tags if several. Default: TRUE
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
# express.thold - Threshold used to define expression-related clasess as described above. Default: 1
# use.dp = Logical, indicates whether only the x and y values of cells classified as double positives (DP) should be used to calculate 2d kernel density values. If TRUE, the density value for the cells belonging to any other class is set to 0. Default: FALSE.
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# this.color.scale - Character. Color scale or color values to use. To be implemented, i.e., is totally ignored for now.
# Dependencies ---------------------->
# Function: get.densities to calculate 2D kernel densities for each cell.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# To do ----------------------------->
# Consider a specific color scale.
# Expand 'keep' argument as was done for the rest of the arguments it's associated to (i.e., 1 to 1 correspondance for all those arguments).
# Function:

coexpress.plot <- function(seurat.obj, feature.x, feature.y, slot='data', filter.tags=NULL, na.rm=TRUE, groups.to.filter=NULL, keep=TRUE, feature.thold=NULL, express.thold=1, use.dp=FALSE, file.name=NULL, this.color.scale=NULL){
  # ---> Arguments and checks.
  if(is.null(this.color.scale)) this.color.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000') # To implement.
  if(is.null(slot)) slot <- 'data'
  both.feats <- c(feature.x, feature.y)
  # Define CPMs-related arguments..
  if(slot == 'cpm'){
    slot <- 'counts'
    cpm <- TRUE
  }else{
    cpm <- FALSE
  }
  # Filter tags and groups to filter dimensions must agree with each other.
  if(!is.null(filter.tags)){
    if(length(filter.tags)>1){
      if(!is.list(groups.to.filter) | length(filter.tags)!=length(groups.to.filter)) stop(paste0('When multiple tags are provided for cell filtering, the argument groups.to.filter must be a list with the same dimension length.\n'))
    }else{
      if(!is.list(groups.to.filter)){
        if(!is.vector(groups.to.filter)) stop('Argument "groups.to.filter" must be either a list (with the same length as that for the argument "filter.tags") or a vector.\n') else groups.to.filter <- list(groups.to.filter)
      }
    }
    names(groups.to.filter) <- paste0('filter.tag.', 1:length(filter.tags))
  }
  # Feature thresholds (for cells filtering). Get a 2-length vector if a single, non-NULL value provided.
  if(!is.null(feature.thold)){
    if(length(feature.thold)==0 & length(feature.thold)>2) stop('Feature threshold must have length > 0 and < 3.\n')
    if(length(feature.thold)==1) feature.thold <- c(feature.thold, feature.thold)
  }
  # Expression-related threshold. Must be a single value (positive scalar).
  if(length(express.thold)!=1 | express.thold<=0) stop('Expression threshold (argument \'express.thold\') value must be a single value (positive scalar).\n')
  # ---> Fetch data.
  # Look for feature.
  # Make sure we've got a valid slot and try to find the feature there.
  if(!(slot %in% c('counts', 'data', 'scale.data'))) stop('Not valid slot. Should be anyone among "counts", "data", "scale.data".\n')
  # Then, make sure we've got both features defined in whatever data slot that was defined.
  tmp.assess <- both.feats %in% rownames(GetAssayData(object=seurat.obj, slot=slot, assay='RNA'))
  if(!all(tmp.assess)){
    tmp.assess <- both.feats[!tmp.assess]
    stop(paste0('Not all features defined in slot requested, ', slot, ', as follows:\n', paste0(tmp.assess, collapse='\n'), '\n'))
  }
  # Proceed if all good.
  if(cpm){
    feature.vals <- RelativeCounts(data=GetAssayData(object=seurat.obj, slot=slot, assay='RNA'), scale=1e6, verbose=FALSE)
    feature.vals <- log2(t(as.matrix(feature.vals[both.feats, ])) + 1)
    feature.vals <- as.data.table(feature.vals)
  }else{
    feature.vals <- as.data.table(FetchData(object=seurat.obj, vars=both.feats, slot=slot))
  }
  colnames(feature.vals) <- c('feature.x', 'feature.y')
  # Fetch the data for the rest of the tags after checking everything exists.
  if(!is.null(filter.tags)){
    tmp.assess <- filter.tags %in% colnames(seurat.obj@meta.data)
    if(!all(tmp.assess)){
      tmp.assess <- filter.tags[!tmp.assess]
      stop(paste0('Tags for filtering set to a value different than NULL, but some of them not defined in seurat object\'s metadata, as follows.\n', paste0(tmp.assess, collapse='\n'), '\n'))
    }
    tmp.data <- data.table(seurat.obj@meta.data[, filter.tags])
    colnames(tmp.data) <- paste0('filter.tag.', 1:ncol(tmp.data))
  }
  # If necessary, join altogether.
  if(!is.null(filter.tags)) tmp.data <- cbind(tmp.data, feature.vals) else tmp.data <- feature.vals
  # ---> Filter out non-desired cells.
  if(!is.null(filter.tags)){
    # Filter all cells with not available data regarding all column-to-filter information.
    if(na.rm){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        tmp.data <- tmp.data[!is.na(get(tmp.col))]
      }
    }
    # Then, when requested, keep or remove input values.
    if(!is.null(groups.to.filter)){
      for(idx in 1:length(filter.tags)){
        tmp.col <- paste0('filter.tag.', idx)
        if(keep) tmp.data <- tmp.data[get(tmp.col) %in% groups.to.filter[[tmp.col]]] else tmp.data <- tmp.data[!get(tmp.col) %in% groups.to.filter[[tmp.col]]]
      }
    }
  }
  # ---> Apply feature threshold.
  if(!is.null(feature.thold)){
    if(!is.numeric(feature.thold)) stop('Not valid feature threshold provided. Must be a numeric value.\n')
    tmp.data <- tmp.data[feature.x > feature.thold[1]]
    tmp.data <- tmp.data[feature.y > feature.thold[2]]
  }
  # ---> Define expression-related classes.
  # @  Define classes in general data,
  # Double positives.
  tmp.data[feature.x>express.thold & feature.y>express.thold, express.class:='DP']
  # Single positives (either for x or y).
  tmp.data[feature.x>express.thold & feature.y<=express.thold, express.class:='SP.X']
  tmp.data[feature.x<=express.thold & feature.y>express.thold, express.class:='SP.Y']
  # Double negatives.
  tmp.data[feature.x<=express.thold & feature.y<=express.thold, express.class:='DN']
  # @  Calculate fraction of each class.
  class.counts <- tmp.data[, .(class.fraction=.N/tmp.data[, .N]), by=express.class]
  # @  Get caption on classess' fractions.
  tmp.caption <- class.counts[, paste(express.class, round(x=class.fraction, digits=2), sep=' - ', collapse='; ')]
  # ---> Calculate densities.
  # Depends on whether all cells should be used or only double positive cells.
  if(use.dp){
    dens.no <- min(c(tmp.data[express.class=='DP', .N], 100))
    tmp.data[express.class=='DP', density:=get.densities(x=tmp.data[express.class=='DP', feature.x], y=tmp.data[express.class=='DP', feature.y], n=dens.no)]
    tmp.data[is.na(density), density:=0]
  }else{
    dens.no <- min(c(tmp.data[, .N], 100))
    tmp.data[, density:=get.densities(x=tmp.data[, feature.x], y=tmp.data[, feature.y], n=dens.no)]
  }
  # ---> Scatter plot depicting density, a.k.a. Co-expression plot.
  # Get ggplot.
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=feature.x, y=feature.y)) + geom_density2d(colour = "#c0c5ce") + geom_point(aes(col=density), size=0.2) + viridis::scale_color_viridis(option = "magma") + labs(x=feature.x, y=feature.y, col='Density', caption=tmp.caption)
  # ---> Return.
  # If reports path is set to null, return ggplot.
  if(is.null(file.name)){
    return(tmp.ggplot)
  }else{
    pdf(file=file.name)
    print(tmp.ggplot)
    dev.off()
    return(NA)
  }
}

# ------------ IN PROCESS

# 5 -------------------------------------------------------------------->
# Name: Co-expression plot.
# Description:
# Given a seurat object, this function will:
# 1. Fetch the appropriate data columns from it.
# 2. If indicated (filter.tag set to a -valid- value other than NULL), filter out or forward (keep set to F or T, respectively) a set of groups (groups.to.filter) of a tag (filter.tag).
# 3. Define expression-related clasess for each cell by checking whether expression for each feature is larger than a set threshold (express.thold). Possible classes for each cell are: Double positive (DP, expresses both features over the threshold); Double negative (DN, expresses none of the features over the threshold); and single positive for either feature (SP.X or SP.Y when expressing the feature depicted in the x coordinate or the y coordinate of the scatter plot, respectively).
# 4. Calculate the fraction of each expression-related class to depict on the final plot.
# 5. Calculate the 2D kernel density for each cell.
# 6. Generate a ggplot object depicting a scatter plot with density as color, which will be either returned (when file.name set to NULL) or output (when a valid file.name provided).
# Arguments ------------------------->
# seurat.obj - seurat object. The tags indicated should already be defined. No default.
# feature.x & feature.y - Features for which to depict their expression values in scatter plot (x and y axes, respectively). Both should be defined as part of a valid data slot (choose slot accordingly if any prefered -'data' as default-). No default.
# slot - Valid data slot to fetch data from. NULL indicates 'data', 'cpm' indicates log2(CPM + 1) relative counts, while the rest are the usual ones: counts, data or scale.data. No other value is allowed. Default: data.
# filter.tags - Tags defined in the seurat onject metadata whose groups will be used to filter out cells. No cell filtering is applied when it's set to NULL. Default: NULL.
# na.rm - Logical, indicates whether NA values defined the filter.tag column should be discarded for the tag listing groups. Default: TRUE.
# groups.to.filter - Character or list, listing the set of tag groups that'll be kept ot removed according to the argument 'keep' value. This step is skipped is it's set to NULL. When more than one tag are provided for filtering, this argument must be a list with the same dimensions as that argument. Default: NULL.
# keep - Logical indicating if groups defined in 'groups.to.filter' should be filter out or forward when set to FALSE and TRUE, respectively. Recycled for all filtering tags if several. Default: TRUE
# feature.thold - Numeric, must be a number put as a threshold to filter forward only certain cells.
# express.thold - Threshold used to define expression-related clasess as described above. Default: 1
# use.dp = Logical, indicates whether only the x and y values of cells classified as double positives (DP) should be used to calculate 2d kernel density values. If TRUE, the density value for the cells belonging to any other class is set to 0. Default: FALSE.
# file.name - Though NULL indicates that the ggplot should be returned instead, this value indicates the name of the file to output plot to. Default: NULL.
# this.color.scale - Character. Color scale or color values to use. To be implemented, i.e., is totally ignored for now.
# Dependencies ---------------------->
# Function: get.densities to calculate 2D kernel densities for each cell.
# Value ----------------------------->
# Either a ggplot object or NULL value and the plot is output as a file (when file.name is different than NULL)
# To do ----------------------------->
# Consider a specific color scale.
# Expand 'keep' argument as was done for the rest of the arguments it's associated to (i.e., 1 to 1 correspondance for all those arguments).
# Function:

# upset.plot <- function(data, is.seurat=FALSE, ){
#   # ---> Prepare general input data.
#   # In case we start from a seurat object, take the meta data as a data table.
#   if(is.seurat) data <- as.data.table(data@meta.data)
#   # Get tidy data while keeping only clonotypes seen in clonally expanded cells.
#   tmp.data.1 <- data[!is.na(clonotype.tag), .(count=.N), by=.(clonotype=clonotype.tag, cluster=RNA_snn_res.0.2)]
#   tmp.data.2 <- data[!is.na(clonotype.tag), .(count=.N), by=.(clonotype=clonotype.tag)][count>1]
#   tmp.data <- tmp.data.1[clonotype%chin%tmp.data.2[, clonotype]]
#   rm(tmp.data.1); rm(tmp.data.2)
#   tmp.data <- spread(data=tmp.data, key=cluster, value=count, fill=0, sep='.')
#   tmp.data <- as.data.frame(tmp.data, stringsAsFactors=FALSE)
#   # Get sets.
#   upset.sets <- paste0('cluster.', data[, as.character(unique(get(tf.pt.clusts.lab)))])
#   all(upset.sets %in% colnames(tmp.data))
#   tmp.data <- as.matrix(tmp.data[, upset.sets])
#   tmp.idxs <- tmp.data>1
#   tmp.data[tmp.idxs] <- 1
#   tmp.data[!tmp.idxs] <- 0
#   tmp.data <- as.data.frame(tmp.data)
#   for(tmp.col in upset.sets) tmp.data[, tmp.col] <- as.integer(tmp.data[, tmp.col])
#   # Output UpSet plot.
#   tmp.file.name <- paste0(tmp.reports.path, '/Upset_ClonotypeSharingAmongClusters.pdf')
#   pdf(file=tmp.file.name)
#   print(UpSetR::upset(
#     data=tmp.data,
#     sets=upset.sets,
#     mainbar.y.label='Clonotype sharing',
#     sets.x.label='Clonotypes per cluster',
#     main.bar.color='#0080ff',
#     sets.bar.color='#0d0d0d',
#     text.scale=c(1, 1, 1, 1, 1, 1), # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set
#     point.size=2,
#     line.size=0.7
#   ))
#   dev.off()
# }
