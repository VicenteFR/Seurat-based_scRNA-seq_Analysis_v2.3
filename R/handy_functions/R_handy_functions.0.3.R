### --------------------------- Libraries --------------------------- ###
library(dplyr)
#library(Seurat)

### ------------------------ Other variables ------------------------ ###
binary.col.scale <- c('yellow'='#cccc00', 'blue'='#0074e6')

# --------------------------------------- 1 --------------------------------------- #
# Name: Load Packages Set ---------------------------------------------------------->
# Description: Given a set of packages, the function will load them to R session and report any missing libraries.
# Calling: load_pacakages(packages, )
# Ciro's equivalent function: load_packs()
# Arguments:
# packages - Vector listing set of packages.
# lib - Directory to libraries.
# v - Boolean indicating verbose mode for this function.

# Function
load_packages <- function(packages, lib=.libPaths(), v=TRUE){
  options(warn=-1) # To suppress warnings in the function environment.
  if(v) cat('*******************************************************************\n')
  .libPaths(new = lib[1]) # Set library path as that of the environment from where function was called.
  if(v) cat('Path to Libraries Directory:', .libPaths()[1], '\n\n')
  if(v){cat('Loading', length(packages), 'package(s) as listed below:\n\t'); cat(packages, sep='\n\t')}
  if(v) prog.bar <- txtProgressBar(min = 0, max = length(packages), width = NA,  style = 3)
  req.out <- vector() # Output for the require function run for each package below.
  for(i in 1:length(packages)){
    if(v) setTxtProgressBar(prog.bar, i)
    req.out <- append(req.out, suppressMessages(require(packages[i], character.only = T))) # Even though this suppresses messages, it will still output errors and warnings.
  }; if(v) cat('\n'); if(v) close(prog.bar)
  req.out <- as.logical(req.out) # Even though the output of require is a logical value, the process through suppressMessages turns it into a character value, so here it's converted again into logical.
  if(v) cat(sum(req.out), 'package(s) loaded.\n')
  if(sum(req.out) < length(packages)){ # To check missing packages
    packages <- packages[!req.out] # packages variable becomes a vector listing the missing ones.
    if(v){ cat('Missing packages:\n'); cat(packages, sep = '\n\t'); cat('\nTry to install them.\n') }
  }
  else{cat('That is, they were all loaded appropriately!\n')}
  if(v) cat('*******************************************************************\n')
}

# 1 -------------------------------------------------------------------->
# Name: Get file name
# Description:
# Given the arguments required, this outputs a character vector as the general path to the file describing its content as well as the important arguments used in the program and others as required.
# Arguments ------------------------->
# gen.pat - Absolute path where the file will be saved.
# file.dess - Description of the file's content that, for clarity, is expected to be formatted withoud blank spaces and with upper case letters indicating them instead.
# extension - Chracter vector indicating the extension of the file.
# added.args - List of character vectors indicating the names of the variables saving other arguments you desire to add to the file name description.
# Function:

get.file.name <- function(gen.path, file.desc, extension, added.args=NULL){
  tmp.file.name <- paste0(gen.path, '/', file.desc)
  if(!is.null(added.args)){
    added.args <- paste(sapply(X=added.args, FUN=function(x){
      paste0('_', x, '_', as.character(get(x)))
    }), collapse='')
    tmp.file.name <- paste0(tmp.file.name, added.args, '.', extension)
    return(tmp.file.name)
  }
  tmp.file.name <- paste0(tmp.file.name, '.', extension)
  return(tmp.file.name)
}

# 2 -------------------------------------------------------------------->
# Name: Get features metadata
# Description:
# Returns a simple dataframe for a given surat object's features metadata, with features' means and variances across cells as the dataframe columns named feat.mean and feat.var, respectively.
# Arguments ------------------------->
# seurat.obj - Seurat object to retrieve features' metadata.
# Function:

get.feats.metadata <- function(seurat.obj){
  feat.mean <- apply(X=seurat.obj@assays$RNA@data, MARGIN=1, FUN=mean)
  feat.var <- apply(X=seurat.obj@assays$RNA@data, MARGIN=1, FUN=var)
  return(data.frame(feat.mean, feat.var))
}

# 3 -------------------------------------------------------------------->
# Name: Analysis of the relation between features' mean and variance.
# Description:
# Given a dataframe describing the simple metadata of a set of features, this program saves to a file (in the general path indicated for the function) a scatter plot of the relation mean vs variance and returns a vector listing the amount of lowly-expressed genes discarded by two arbitrary approaches (these also depicted in the scatter plot) as described below:
# Approach A: Genes whose expression level is below 1TPM or its equivalent given the normalization factor with LogNormalize.
# Approach B: Those with a level of expression below the equivalent to the first quantile calculated after removing non-expressed genes.
# Arguments ------------------------->
# feats.metadata - Dataframe storing the normal simple metadata of a set of features.
# norm.factor - The same used to normalize data with method LogNormalize.
# data.set - The name of the dataset the features are coming from.
# gen.path - Absolute path where the scatter plot will be saved.
# Function:

mean.by.var.analysis <- function(feats.metadata, norm.factor, data.set, gen.path){
  total.no.feats <- nrow(feats.metadata)
  lt.zero <- feats.metadata$feat.mean > 0
  TPM.threshold <- norm.factor/10e6
  lt.TPM.threshold <- feats.metadata$feat.mean > TPM.threshold
  # Turn variables into their log2 values for plotting purposes.
  # Approach A thresholds to plot.
  TPM.threshold <- log2(TPM.threshold) # Will be the x threshold.
  A.y.threshold <- log2(min(feats.metadata$feat.var[lt.TPM.threshold]))
  # Features discarded by approach A conditions.
  discarded.by.A <- total.no.feats - sum(lt.TPM.threshold)
  # Approach B thresholds to plot.
  first.quant <- quantile(x=feats.metadata$feat.mean[lt.zero], probs=0.25) # Will be the x threshold.
  lt.1st.quant <- feats.metadata$feat.mean > first.quant
  B.y.threshold <- log2(min(feats.metadata$feat.var[lt.1st.quant]))
  # Features discarded by approach B conditions.
  discarded.by.B <- total.no.feats - sum(lt.1st.quant)
  # Building plot.
  tmp.caption <- paste0('Total no. feats:', total.no.feats, ';discarded by A:', discarded.by.A, ';by B:', discarded.by.B)
  mean.by.var.plot <- ggplot(feats.metadata, aes(x=log2(feat.mean), y=log2(feat.var))) +
  geom_point(size=0.5) +
  labs(title='Genes\' mean by variance plot', subtitle=data.set, x='log2 Mean', y='log2 Variance', caption=tmp.caption) +
  annotate('rect', xmin=-Inf, xmax=TPM.threshold, ymin=-Inf, ymax=A.y.threshold, alpha=0.2) + annotate('text', x=TPM.threshold+2, y=A.y.threshold, label='Approach A') +
  annotate('rect', xmin=-Inf, xmax=log2(first.quant), ymin=-Inf, ymax=B.y.threshold, alpha=0.4) + annotate('text', x=log2(first.quant)+2, y=B.y.threshold, label='Approach B')
  # Save to a file the plot.
  tmp.file.name <- paste0(gen.path, '/MeanByVariancePlotFor', data.set, '.pdf')
  pdf(file=tmp.file.name)
  print(mean.by.var.plot)
  dev.off()
  # Finally, return the amount of genes discarded by each approach.
  output <- c(discarded.by.A, discarded.by.B, first.quant)
  names(output) <- c('Approach A', 'Approach B', '1st quantile')
  return(output)
}

# 4 -------------------------------------------------------------------->
# Name: Cummulative variance
# Description:
# Given a simple features metadata dataframe, this function will sort all features in decreasing order according to variance and will output it along with the cummulative sum and its percentage out of the total variance. A threshold value can be passed so that features below it are not considered at all. Moreover, it also saves to a general path a scatter plot and a table depicting the cummulative variance for the dataset under consideration.
# Arguments ------------------------->
# feats.metadata - A simple features metadata dataframe (i.e., with feat.mean and feat.var as columns among other possible ones).
# threshold - Threshold to set on features variance.
# gen.path - Absolute path to save both cummulative variance plot and table.
# data.set - Name of the data set for which the process will be computed.
# Function:

cummulative.variance <- function(feats.metadata, threshold=0, gen.path, data.set){
  lt.threshold <- feats.metadata$feat.mean > threshold
  feats.var <- sort(feats.metadata$feat.var[lt.threshold], decreasing=TRUE)
  total.var <- sum(feats.var)
  cumm.var <- cumsum(feats.var)
  pct.cumm.var <- (cumm.var*100)/total.var
  variance.summ <- data.frame(feats.var, cumm.var, pct.cumm.var)
  rownames(variance.summ) <- rownames(feats.var)
  # Create scatter plot.
  cummulative.variance.plot <- ggplot(data.frame(index=1:length(cumm.var), cumm.var), aes(x=index, y=pct.cumm.var)) + geom_point(size=1) + labs(title='Cummulative plot', subtitle=data.set, x='Sorted variable genes', y='% cummulative variance', caption='Notice that features have been filtrated.')
  tmp.file.name <- paste0(gen.path, '/CummulativeVariancePlot4', data.set, '.pdf')
  pdf(file=tmp.file.name)
  print(cummulative.variance.plot)
  dev.off()
  # Return cummulative variance summary after saving it to a cvs file.
  tmp.file.name <- paste0(gen.path, '/CummulativeVarianceTable4', data.set, '.csv')
  write.csv(x=variance.summ, file=tmp.file.name)
  return(variance.summ)
}

# 5 -------------------------------------------------------------------->
# Name: are.cells.in.this.population
# Description:
# Given a table of tags' expression values, this function checks for the expression values (either positive or negative) of a set of these markers. Moreover, it can identifiy if certain sets show expression while others no. In the end, based on marker expression, this function identifies cells' populations.
# Arguments ------------------------->
# tags.exp.values - Data frame from which check the expression of the markers of interest set as arguments below. No Default.
# markers.sets - List of different char vectors representing different sets markers with the restriction that either 1) all the members of a vector show the expression signs indicated in the sets.exp.signs argument or 2) at least one does. This is stated in the any.or.all argument. No Default.
# sets.exp.signs - List of char vectors indicating the expression sign (+ or -) for their corresponding set of markers according to their order.
# any.or.all - Char vector (value 'any' or 'all') where index i indicates if 1) any or 2) all genes should be expressed for the set in the index i of the list markers.sets. Default 1.
#logical.rels - vector of length equivalent to that of the amount of sets minus 1 to evaluate all sets logical relationships (values must be 'and' and 'or'). Non valid values will be deleted and then, if length(logical.rels) is not the specified or equals to 1, function will stop. Logical values will be avaluated according to the order of the markers sets list. No default. Also, notice that if just one markers set is provided, then, this argument value will just be ignored (i.e., will be turned into NULL).
# log.or.tag - Value either 'log' (logical) or 'tag' defining if a cell should be defined as TRUE/FALSE (logical) for a given population or should be defined by the value given to population.tag.
# output.opt - Integer with either value 1, 2 or 3 indicating if the vector result (population+ or population- for each cell) should be returned as an independent vector, as the vector pasted to the original tags.exp.values dataframe or as it attached to a seurat.obj metadata dataframe, respectively. If 3, seurat.obj argument cannot be NULL. Default 1.
# seurat.obj - Seurat object to add the vector result to. If output.opt is set to 3, this argument cannot be NULL or function will output an error. Default NULL.
# population.tag - Char indicating the name to assign to the population if this is added to the tags dataframe or the seurat object. If NULL, it is equal to the concatenation of all sets with their logical relationships. Default NULL.
# Example:
# are.cells.in.this.population(tags.exp.values=tags.exp.values, markers.sets=list('FOXP3', c('GZMB', 'PRF1')), sets.exp.signs=list('-', c('+', '+')), any.or.all=c('all', 'any'),logical.rels='and', output.opt=3, seurat.obj=seurat.obj, population.tag='CD4CTLs')
# Applied in this way to a set of CD4 T cells, this function would identify those non-Treg-like CD4 CTLs.
# Function:

are.cells.in.this.population <- function(tags.exp.values, markers.sets, sets.exp.signs, any.or.all='all',logical.rels=NULL, log.or.tag='tag', output.opt=1, seurat.obj=NULL, population.tag=NULL){
  # Check that the lengths of sets of markers and signs lists are equivalent.
  if(!(length(markers.sets)==length(sets.exp.signs))) stop('You should give the same amount of expression values signs sets as that of markers sets.')
  # Output error if no seurat object is provided when it is wanted to be part of the output.
  if(output.opt==3 & is.null(seurat.obj)) stop('output.opt set to 3 but no seurat object provided as argument.')
  # Same but for any.or.all argument
  if(!(any.or.all=='any' | any.or.all=='all')) stop('Non valid value for any.or.all argument')
  # If just one set of markers is introduced and valid logical relationships are passed as arguments, the function will simply turn it into NULL. Alert user as a warning.
  if(length(markers.sets)==1 & !is.null(logical.rels)) warning('Logical relationships passed as argument when there is only one markers set. These values will ignored (i.e., logical.rels argument value will be turned into NULL).')
  # Function to apply iteratively for lists of sets and their expression values.
  to.do.mapply <- function(markers.set, set.exp.signs, any.or.all){
    markers.tags <- paste0(markers.set, '.tag')
    tags.exp.values <- tags.exp.values[, markers.tags]
    tags.exp.pattern <- paste0(markers.set, set.exp.signs)
    if(is.null(dim(tags.exp.values))){
      return(tags.exp.values==tags.exp.pattern)
    }
    else{
      set.eval <- t(apply(X=tags.exp.values, MARGIN=1, FUN=function(row){return(row==tags.exp.pattern)}))
      if(any.or.all=='any') return(apply(X=set.eval, MARGIN=1, FUN=any)) else return(apply(X=set.eval, MARGIN=1, FUN=all))
    }
  }
  sets.evals <- mapply(FUN=to.do.mapply, markers.sets, sets.exp.signs, any.or.all)
  # If logical relationships are enough to compare the different markers sets, the function will be applied normally. If it's an atomized value (length=1), it will be recycled unless it's NULL, when there should be no value. Any other value will just be rejected.
  if(length(logical.rels)>0 & length(markers.sets)==1) logical.rels <- NULL
  if(!(length(logical.rels) == ncol(sets.evals) - 1 | length(logical.rels) == 1)) stop('Non-valid logical relationships passed for logical.rels argument.')
  # Then, logical values are passed to the right syntax.
  logical.rels <- gsub(pattern='and', replacement='&', x=logical.rels, fixed=TRUE)
  logical.rels <- gsub(pattern='or', replacement='|', x=logical.rels, fixed=TRUE)
  # When there's only one logical relationship, this will be repeated times the number of comparisons needed (i.e., the number of sets minus one). When it's NULL (i.e., when there's just one set), logical.rels will be empty (which is good for the fuction goal).
  if(length(logical.rels) == 1) logical.rels <- rep(x=logical.rels, times=ncol(sets.evals)-1)
  # For the porpuses of using paste appropriately.
  logical.rels <- c(logical.rels, '')
  # Finally, the sets evaluations are pasted with their logical relationships in order to furtehr evaluate them and figure out if a cell is part of the population chracterized by the markers under consideration.
  are.cells.in.this.population <- apply(X=sets.evals, MARGIN=1, FUN=function(cell.evals){
    to.eval.set.rels <- paste0(cell.evals, logical.rels, collapse='')
    return(eval(expr=parse(text=to.eval.set.rels)))
  })
  # Create a tag (i.e, a name) for the population if it's not defined already.
  if(is.null(population.tag)) population.tag <- paste0(unlist(markers.sets), unlist(sets.exp.signs), collapse='')
  # Then, based on the choice of the user (or 'tag' by default), either we leave the values of the cells in the population as logical values ('log') or we replace the logical values with the same tag given for the population (non-'tag' and 'tag').
  if(log.or.tag=='tag'){ are.cells.in.this.population <- gsub(x=are.cells.in.this.population, pattern='TRUE', replacement=population.tag, fixed=TRUE); are.cells.in.this.population <- gsub(x=are.cells.in.this.population, pattern='FALSE', replacement=paste0('non-', population.tag))}
  # In the end, the result is output according to the option.
  if(output.opt==1){
    return(are.cells.in.this.population)
  }
  if(output.opt==2){
    tags.exp.values[, population.tag] <- are.cells.in.this.population
    return(tags.exp.values)
  }
  else{
    seurat.obj@meta.data[, population.tag] <- are.cells.in.this.population
    return(seurat.obj)
  }
}

# 6 -------------------------------------------------------------------->
# Name: find.tags.per.cell
# Description:
# Given the expression values per population tag (i.e., columns in metadata defined by tags pop.tag and 'non-'pop.tag) for a set of cells and also given a set of population tags of interest, this function returns a char vector specifying per cell if it 1) belongs to a unique population (defined by the population tag), 2) belongs to more than one population (defined by the many.tags tag) or 3) belongs to none population tag of interest.
# Arguments ------------------------->
# tags.exp.values - Data frame from which check the expression of the population tags of interest set as arguments below. No Default.
# tags.to.eval - A set of population tags (char vector) of interest to evaluate each cell finding if it belongs to it or not. If any of the population tags is not defined in the dataframe tags.exp.values, the funciton will stop.
# none.tag - String defining the tag for cells that belong to none population of interest.
# many.tags - String defining the tag for cells that are part of more that one population of interest.
# prio - Logical indicating if tags should be prioritized according to their order in tags.to.eval vector (thereby, not returning an overlapping population ever) or not.

find.tags.per.cell <- function(tags.exp.values, tags.to.eval, none.tag='None', many.tags=NULL, prio=FALSE){
  # Make sure all tags have already been defined.
  all.tags <- colnames(tags.exp.values)
  if(sum(!(tags.to.eval%in%all.tags))>0) stop('Not all tags defined in the dataframe passed to tags.exp.values.')
  # First, in case it wasn't done before passing tags.exp.values as an argument, subset the dataframe to keep tags of interest.
  tags.exp.values <- tags.exp.values[, tags.to.eval]
  # Then, look for what populations are defined per cell.
  tags.per.cell <- apply(X=tags.exp.values, MARGIN=1, FUN=function(cell.tags){
    # Which tags are defined.
    tags.found <- grep(pattern='non-', x=cell.tags, invert=TRUE)
    # Return tags accoridng to which were found.
    # If many tags were found.
    if(length(tags.found)>1){
      if(prio) return(tags.to.eval[tags.found[1]])
      if(is.null(many.tags)) many.tags <- paste(tags.to.eval[tags.found], collapse='-')
      return(many.tags)
    }
    # If none was found or just a non-overlapping population was found.
    if(length(tags.found)==0) return(none.tag) else return(tags.to.eval[tags.found])
  })
  return(tags.per.cell)
}

# 7 -------------------------------------------------------------------->
# Name: read.seurat.input.data
# Description:
# Given a file/path name, this function figures out what its format is to work as a data object appropriate to become a seurat object. It looks for three formats
# a) TSV with extensions .tsv and .txt
# b) CSV with extension csv
# c) H5 with extension h5
# The default is for it to be taken as a path to cellranger output where matrix, features and cells data is saved.
# The function returns the file read and saved to an object.
# Arguments ------------------------->
# file - Character value to be assessed for any of the different formats described.
# feature.id - Character with possible values 'name' or 'ensembl'. Used just in case a cellranger output is provided (either from features.csv file or from h5 matrix). Default, name.

read.10X.data <- function(file, feature.id='name'){
  if(!is.null(dim(file))) stop('file argument must be a character vector storing in its first index the name of the file to be read.')
  file <- file[1]
  # ---> Check feature.id option value.
  # Valid values: 'name', 'ensembl'
  if(!(feature.id=='name'|feature.id=='ensembl')) stop('Non valid entry for feature.id argument.')
  # ---> Read data.
  #   ---> Normal table
  # Tab-separated values.
  if(grepl(x=file, pattern='\\.txt$', ignore.case=TRUE, perl=TRUE) | grepl(x=file, pattern='\\.tsv', ignore.case=TRUE, perl=TRUE)) return(read.table(file=file))
  # Comma-separated values.
  if(grepl(x=file, pattern='\\.csv$', ignore.case=TRUE, perl=TRUE)) return(read.csv(file=file, row.names=1))
  #   ---> Cellranger output.
  # Case h5 format.
  if(grepl(x=file, pattern='.h5$', ignore.case=TRUE, perl=TRUE)){
    feature.id <- feature.id=='name'
    return(Seurat::Read10X_h5(filename=file, use.name=feature.id))
  }
  # Directory with 10X output (default).
  return(Seurat::Read10X(data.dir=file, gene.column=ifelse(test=feature.id=='name', yes=2, no=1)))
}

# 7 -------------------------------------------------------------------->
# Name: sample.fav.colors
# Description:
# Not really a function. This just saves my favorite colors so that you can get a sample of them.
# Arguments ------------------------->
# size - For sample function.

favorite.colors <- c('firebrick', 'cornsilk3', 'aquamarine4', 'darkgoldenrod1', 'cyan4', 'mediumvioletred', 'lightsalmon3', 'khaki', 'lightblue4', 'darkseagreen', 'dodgerblue4', 'dimgrey', 'coral3', 'aliceblue', 'azure1', 'bisque3', 'brown', 'gray40', 'midnightblue', 'orchid4')
sample.fav.colors <- function(size){
  return(sample(x=favorite.colors, size=size, replace=TRUE))
}

# 8 -------------------------------------------------------------------->
# Name: Create directory
# Description:
# Given an absolute path to a directory, this function checks if it exists and if not, it creates it.
# Useful for bigger pipelines reports.
# Arguments ------------------------->
# dir.path - Absolute path of a diretcory to create/check if exists.
# path.desc - A brief description/name of the absolute path.
# Function:

create.dir <- function(dir.path, path.desc=NULL){
	if(dir.exists(dir.path)){
		cat(path.desc, 'directory already exists.\n')
	}else{
		dir.create(dir.path)
    if(is.null(path.desc)) path.desc <- basename(dir.path)
		cat(path.desc, 'directory has been created. Check warnings if any.\n')
	}
}

# 9 -------------------------------------------------------------------->
# Name: read.gen.R.obj
# Description:
# Given an absolute path to a file, this function identifies its format either from the RDS or RData extensions (case insensitive) and reads appropriately.
# Useful for bigger pipelines reports.
# Arguments ------------------------->
# obj.file - Absolute path to R data file.
# obj.name - Character vector or NULL values indicating either the name of the object or none, respectively.
# v - Logical indicating if the process should process verbosely.
# Function:

read.gen.R.obj <- function(obj.file, obj.name=NULL, v=FALSE){
  env <- new.env()
  if(v) cat('Reading object ')
  if(grepl('rds$', obj.file, ignore.case=TRUE)){
    cat("RDS\n")
    env$lol <- readRDS(obj.file)
  }else{ cat("RDATA\n"); load(obj.file, envir = env) }
  loadedObjects <- objects(env, all = TRUE)
  if(is.null(obj.name)){
    if(length(loadedObjects) > 1) stop(length(loadedObjects), ' objects in data, specify one')
    obj.name <- loadedObjects
  }else if(is.numeric(obj.name)){ # mmda para hacerlo bonito
    if(obj.name > length(loadedObjects)){ # if there fewer objects than the asked number
      if(grepl('1$', as.character(obj.name))) tvar <- paste0(obj.name, 'st')
      if(grepl('2$', as.character(obj.name))) tvar <- paste0(obj.name, 'nd')
      if(grepl('3$', as.character(obj.name))) tvar <- paste0(obj.name, 'rd')
      if(sum(c('12', '13', '11') %in% as.character(obj.name))) tvar <- paste0(obj.name, 'th')
      if(!exists('tvar')) tvar <- paste0(obj.name, 'th')
      stop(length(loadedObjects), ' objects in data, ', tvar, ' not found')
    }
    obj.name <- loadedObjects[obj.name]
  }
  if(v) cat('Taking:', obj.name, '\n')
  if(!sum(obj.name %in% loadedObjects)) stop('Object not found')
  env[[obj.name]]
}

# 10 ------------------------------------------------------------------->
# Name: output.piechart
# Description:
# To add
# Useful for bigger pipelines reports.
# Arguments ------------------------->
# tmp.df - A data frame depicting the relations between groups and their associated frequencies, with columns (notice they have to be called as shown):
#   Group: Names assigned to the groups (to fill in the pie chart).
#   value: Frequency assigned to each group
# file.name - Absolute path to output file. If none is given, the function returns the ggplot.
# text.color - Colot for the text showing proportions. Default: white
# tmp.cols - Colors to fill with the pie chart sections. Default: NULL so that they're taken from 'favorite colors'
# Function:

output.piechart <- function(tmp.df, file.name=NULL, text.color="white", tmp.cols=NULL, piechart.title=NULL){
  # Get total frequency for the groups considered.
  total.freq <- sum(tmp.df$value)
  # Add the proportion for each group
  tmp.df$value.prop <- round(x=tmp.df$value*100/total.freq, digits=3)
  # Add the label y position in pie chart.
  tmp.df <- tmp.df %>% arrange(desc(Group)) %>% mutate(lab.pos=cumsum(value) - 0.5*value)
  # Create a blank theme
  blank_theme <- theme_minimal()+
    theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
    )
  # Get pie chart
  if(is.null(tmp.cols)) tmp.cols <- sample.fav.colors(nrow(tmp.df))
  tmp.plot <- ggplot(data=tmp.df, aes(x="", y=value, fill=Group)) + geom_bar(width=1, stat="identity") +
    coord_polar("y", start=0) + blank_theme +
      theme(axis.text.x=element_blank()) +
      geom_text(aes(y=lab.pos, label=value.prop), color=text.color) +
      scale_fill_manual(values=tmp.cols)
  if(!is.null(piechart.title)) tmp.plot <- tmp.plot + ggtitle(label=piechart.title)
  # Output to a file
  if(!is.null(tmp.file.name)){
    pdf(file=tmp.file.name)
    print(tmp.plot)
    dev.off()
  }else{
    return(tmp.plot)
  }
}

# 11 ------------------------------------------------------------------->
# Name: Write summaries
# Description:
# For different summary objects to be output to a file.
# Arguments ------------------------->
# summs - Char with the summary objects names to be output.
# summs.names - Char with the same length storing the alternative names of the objects to output.
# file.name - Output file.

write.summs <- function(summs, summs.names=NULL, file.name, tmp.env=environment()){
  write(x='Description,Minimum,1st Qu.,Median,Mean,3rd Qu.,Maximum', file=file.name)
  if(is.null(summs.names)){
    for(summ in summs){
      tmp.output <- matrix(data=c(summ, as.character(get(summ))), nrow=1)
      write(x=tmp.output[1, ], file=file.name, sep=',', ncolumns=7, append=TRUE)
    }
  }
  else{
    if(length(summs)!=length(summs.names)){warning('Different amounts of summary objects than that of names.'); return('Check warning!')}
    for(i in 1:length(summs.names)){
      tmp.output <- matrix(data=c(summs.names[i], as.character(get(x=summs[i], envir=tmp.env))), nrow=1)
      #write(x=summs.names[i], file=file.name, append=TRUE)
      write(x=tmp.output[1, ], file=file.name, sep=',', ncolumns=7, append=TRUE)
    }
  }
}

# 12 ------------------------------------------------------------------->
# Name: Get CDR3 sequences from a non-tidy table.
# Description:
# Given the format provided by 10x for single cell TCR data, CDR3 sequences for different chains are usually provided in a single column separated by a semicolon (;). This function takes that table and returns the CDR3 sequences in a character vector.
# Arguments ------------------------->
# cdr3s.table - Non-tidy table storing CDR3 sequences.
# cdr3.col - Name of the column depicting the CDR3 sequences.
# chain.ptn - Pattern to fetch a specific type of chain ('TRA' and 'TRB' for TRA or TRB, respectively, duh) or both (pattern 'TR[A|B]').
# Value ----------------------------->
# Character vector with CDR3 sequences.

get.cdr3s.from.table <- function(cdr3s.table, cdr3.col='cdr3s_nt', chain.ptn='TRB'){
  cdr3s <- unlist(str_split(string=cdr3s.table[, cdr3.col], pattern=';'))
  cdr3s <- grep(x=cdr3s, pattern=chain.ptn, value=TRUE)
  cdr3s <- str_replace(string=cdr3s, pattern='TR[AB]:', replacement='')
  return(cdr3s)
}

# 13 ------------------------------------------------------------------->
# Name: Get CDR3 sequences along info from a non-tidy table.
# Description:
# Just as function 12, but provides metadata (other columns) along each CDR3 sequence.
# Consider that there could be multiple columns for the same clonotype describing different CDR3 sequences (whichever chains are requested).
# Arguments ------------------------->
# cdr3s.table - Non-tidy table storing CDR3 sequences.
# cdr3.col - Name of the column depicting the CDR3 sequences.
# chain.ptn - Pattern to fetch a specific type of chain ('TRA' and 'TRB' for TRA or TRB, respectively, duh) or both (pattern 'TR[A|B]').
# cols.to.fetch - Character describing which columns (as in the original vdj clonotypes table) should be fetched along clonotype CDR3 sequences.
# Value ----------------------------->
# Character vector with CDR3 sequences.

get.cdr3s.info.from.table <- function(cdr3s.table, cdr3.col='cdr3s_nt', chain.ptn='TRB', cols.to.fetch=c('clonotype_id', 'frequency')){
  # Extract the desired pattern for each table row.
  cdr3s <- str_split(string=cdr3s.table[, cdr3.col], pattern=';')
  cdr3s <- lapply(X=cdr3s, FUN=function(tmp.cdr3s) grep(x=tmp.cdr3s, pattern=chain.ptn, value=TRUE))
  cdr3s <- lapply(X=cdr3s, FUN=function(tmp.cdr3s) str_replace(string=tmp.cdr3s, pattern='TR[AB]:', replacement=''))
  # Extract the desired metadata.
  cdr3s.metadata <- cdr3s.table[, cols.to.fetch]
  # Add the same data to each CDR3 sequence.
  cdr3s.info <- lapply(X=1:length(cdr3s), FUN=function(clon.idx){
    clon.cdr3s <- cdr3s[[clon.idx]]
    clon.metadata <- cdr3s.metadata[clon.idx, ]
    clon.cdr3s.info <- lapply(X=clon.cdr3s, FUN=function(tmp.cdr3){
      return(cbind(tmp.cdr3, clon.metadata))
    })
    clon.cdr3s.info <- Reduce(x=clon.cdr3s.info, f=rbind)
    return(clon.cdr3s.info)
  })
  # In case there was no CDR3 sequence for the chain requested (yes, it can happen since a clonotype could be defined by a single TRA or a single TRB), we get rid of NULL values.
  idx.is.null <- sapply(cdr3s.info, FUN=is.null)
  cdr3s.info <- cdr3s.info[!idx.is.null]
  # Finally, bind everything.
  cdr3s.info <- Reduce(x=cdr3s.info, f=rbind)
  colnames(cdr3s.info)[1] <- cdr3.col
  # Important to change the factors as characters.
  cdr3s.info[, 1] <- as.character(cdr3s.info[, 1])
  return(cdr3s.info)
}

# 14 ------------------------------------------------------------------->
# Name: Get Colors From Seurat Object Tag
# Description:
# Given a tag defined in a seurat object, provide colors it would provide for their values.
# Arguments ------------------------->
# seurat.obj - Seurat object with UMAP evaluated.
# tag.of.int - Tag of interest. Must be defined for the seurat object.
# get.names - Logical indicating wether the output vector should be named or not (just when tag of interest is clusters tag). Rationale: If it is named, ggplot scale (color) functions will try to match this to group values; then, make sure to give appropriate input.
# Value ----------------------------->
# Character vector with tag colors.

get.tag.cols <- function(seurat.obj, tag.of.int, get.names=FALSE){
  # Check tag's defined.
  if(!tag.of.int %in% colnames(seurat.obj@meta.data)) stop('Tag not defined.\n')
  # Get colors.
  tmp.ggplot <- UMAPPlot(object=seurat.obj, group.by=tag.of.int) # Generate the tSNE plot, but save it as an object
  build.obj <- ggplot2::ggplot_build(tmp.ggplot) # Use ggplot_build to deconstruct the ggplot object
  build.data <- build.obj$data[[1]] # Pull the data used for the plot
  build.data <-  build.data[order(build.data$group), ] # Order the plot data by group
  tags.cols <- unique(build.data$colour) # Get a vector of unique colors
  if(get.names) names(tags.cols) <- unique(build.data$group)-1 # Add the groups to the vector of colors as names. Minus 1 because clusters are 0-based.
  # Return.
  return(tags.cols)
}

# 15 ------------------------------------------------------------------->
# Name: Get CPM counts from a seurat object info.
# Description:
# Given a seurat object, this function will extract the raw counts from the desired assay and convert them to CPM, returning them as such if desired or returning the same seurat object with a new assay, 'CPM', with the result saved to slot 'data'.
# Arguments ------------------------->
# seurat.obj.query - seurat object.
# counts.assay - Character, assay for counts to use for normalization. Default, 'RNA'
# add.mat - Logical indicating if CPM values should be added to the seurat object as a new assay and returned that way or simply return the CPM values.
# Value ----------------------------->
# Depends. If add.mat=TRUE, same seurat object with CPM values within data slot of a new assay, 'CPM'; else, dgcMatrix with CPM values.

get.cpm.from.obj <- function(seurat.obj.query, counts.assay='RNA', add.mat=FALSE){
  # Get CPM matrix.
  cpm.matrix <- Matrix::t(Matrix::t(GetAssayData(object=seurat.obj.query, slot='counts', assay=counts.assay))*1e6/Matrix::colSums(GetAssayData(object=seurat.obj.query, slot='counts', assay='RNA')))
  # if desired, add to a new assay.
  if(add.mat){
    seurat.obj.query[['CPM']] <- CreateAssayObject(data=cpm.matrix)
    return(seurat.obj.query)
  }
  # else
  return(cpm.matrix)
}

# 16 ------------------------------------------------------------------->
# Name: Apply DGEA based on MAST.
# Dependencies:
# MAST and data.table
# Description:
# Adapted from Seurat package and some code of myself. This function applies MAST-based DEA from a seurat object given a tag of interest annotated in meta data. The tag of interest should be value-binary in the sense that should only have two possible value; then, order of tag values should be provided through main.condition and second.condition arguments.
# This function assumes that data slot countains counts normalized through LogNormalize.
# Arguments ------------------------->
# seurat.obj - seurat object with RNA assay defines and previously log-normalized.
# groups.tag - Tag, already annotated in meta data, listing the groups of interest to be compared.
# main.condition - Reference group for differential analysis.
# second.condition - Group to compare vs reference.
# Value ----------------------------->
# Data table, DEA table result.

do.mast.dea <- function(seurat.obj, tag.of.int, main.condition, second.condition){
  ### ---------------------- Data preprocessing ----------------------- ###
  # Set variables as should be passed to original Seurat MASTDETest function.
  # ----> Seurat object.
  # Subset seurat object to keep just the cells from the conditions to compare.
  cells.1 <- Cells(seurat.obj)[seurat.obj@meta.data[, tag.of.int]==main.condition & !is.na(seurat.obj@meta.data[, tag.of.int])]
  cells.2 <- Cells(seurat.obj)[seurat.obj@meta.data[, tag.of.int]==second.condition & !is.na(seurat.obj@meta.data[, tag.of.int])]
  seurat.obj <- subset(x=seurat.obj, cells=c(cells.1, cells.2))
  # ---> Counts (log-normalized) matrix.
  data.use <- GetAssayData(object=seurat.obj, slot='data', assay='RNA')
  # ---> CDR.
  # Cellular detection rate
  latent.vars <- Matrix::colSums(GetAssayData(object=seurat.obj, slot='data', assay='RNA')>0)
  latent.vars <- scale(data.frame(cdr=latent.vars))
  ### ------------------- Single Cell Assay Object -------------------- ###
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars.names <- c("condition", colnames(x = latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(x = latent.vars)
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(exprsArray=as.matrix(x=data.use), cData=latent.vars, fData=fdat)
  ### --------------- Differential Expression Analysis ---------------- ###
  # ---> Set conditions.
  conditions <- factor(x=SummarizedExperiment::colData(sca)$group)
  conditions <- relevel(x=conditions, ref="Group2")
  SummarizedExperiment::colData(sca)$condition <- conditions
  # ---> Formula for DA.
  fmla <- as.formula(object=paste0(" ~ ", paste(latent.vars.names, collapse = "+")))
  # ---> DA
  zlmCond <- MAST::zlm(formula=fmla, sca=sca, parallel=TRUE)
  # ---> Organize data.
  vs.condition <- 'conditionGroup1'
  summaryCond <- summary(object=zlmCond, doLRT=vs.condition)
  summary.data.table <- summaryCond$datatable
  hurdle.component <- summary.data.table[contrast==vs.condition & component=='H',.(primerid, `Pr(>Chisq)`)]
  lfc.component <- summary.data.table[contrast==vs.condition & component=='logFC', .(primerid, coef, ci.hi, ci.lo)]
  dea.results <- merge(hurdle.component, lfc.component, by='primerid')
  # Add FDR-corrected p-value.
  dea.results[,fdr:=p.adjust(p=`Pr(>Chisq)`, method='fdr')]
  # Return results.
  return(dea.results)
}

# 17 ------------------------------------------------------------------->
# Name: Get annotations for a seurat obejct.
# Version:
#   2.0
# Updates:
#   * Consider NA values witin any aggr table index for the HTO column.
#   * Take into account the possibility of having multiple annotations (separated by semicolon) within the lane.tag column.
# Updates since: functions file version 0.3
# Dependencies:
# Seurat and stringr
# Description:
#   Given an aggregation table with annotations regarding each library (i.e., this is intended for seurat object created based on aggregations), this function will provide the appropriate annotations for each cell in the metadata section.
# Arguments ------------------------->
# seurat.obj - seurat object (raw or any kind)
# Value ----------------------------->
# Annotated seurat object.

annotate.seurat.obj <- function(seurat.obj, aggr.table, lane.tag.id, chrom.batch.tag.id, seq.batch.tag.id, hto.tag.id, overlay.tag.id, mult.ann.per.lane){
  ### ---------------------- Data preprocessing ----------------------- ###
  # ---> Samples suffix per cell.
  # Get sample-specific suffix of cell IDs.
  cell.sffxs <- as.integer(str_extract(string=Cells(seurat.obj), pattern='\\d+$'))
  # ---> Batches post-aggregation.
  # Find amount of 10x batches aggregated based on cell barcode suffixes and make sure it's the same amount of batch definitions provided.
  tenx.batches <- as.integer(unique(str_extract(string=Cells(seurat.obj), pattern='\\d+$')))
  if(length(tenx.batches)!=nrow(aggr.table)) stop('Not same amount of batches identified in Seurat object as the ones provided in table.')
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
    # Unlist sets of HTO tag values and assing identifier according to batch. Then, names will be as in seurat object post-aggregation, composed by barcode and batch.
    # First, get complete cell IDs.
    cell.ids <- unlist(lapply(X=names(hto.tag.vals.list), FUN=function(batch.sfx){
      cell.id <- paste0(hto.tag.vals.list[[batch.sfx]]$barcode, '-',batch.sfx)
      return(cell.id)
    }))
    # Then get hashtag vaues unlisted and assign names.
    hto.tag.vals <- unlist(lapply(X=hto.tag.vals.list, FUN=function(hto.tag.vals) return(hto.tag.vals$hashtag)))
    names(hto.tag.vals) <- cell.ids
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
    # Find HTO tag value per cell.
    hto.tag.vals.per.cell <- hto.tag.vals[Cells(seurat.obj)]
    # Add HTO tag value per cell.
    seurat.obj[[hto.tag.id]] <- hto.tag.vals.per.cell
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

# 18 ------------------------------------------------------------------->
# Name: Subset seurat objects
# Dependencies:
# Seurat and stringr
# Description:
#   Given a seurat object, this program will subset this seurat object according to two kind of rules: 1) Rules per featurs and 2) rules per tags.
# Arguments ------------------------->
# seurat.obj - seurat object.
# tags.criteria - Data frame with criteria to discard info from tags in meta data.
# feats.criteria - Data frame with criteria to discard info from features described in meta data.
# Value ----------------------------->
# Subsetted seurat object.

subset.seurat.obj <- function(seurat.obj, tags.criteria, feats.criteria){
  # ---> Subset seurat object.
  # ---> Assess criteria per feature/tag
  # By tag criteria.
  if(!is.null(tags.criteria)){
    tags.evals <- sapply(X=colnames(tags.criteria), FUN=function(tmp.tag){
      tmp.criteria <- tags.criteria[, tmp.tag]
      to.do <- rownames(tags.criteria)[!is.na(tmp.criteria)]
      tmp.criteria <- tmp.criteria[!is.na(tmp.criteria)]
      tmp.criteria <- str_split(string=tmp.criteria, pattern=";")[[1]]
      # Check if NA value is part of the criteria.
      na.flag <- any(grepl(x=tmp.criteria, pattern='^[\'"]NA[\'"]$', perl=TRUE))
      tmp.criteria <- tmp.criteria[!grepl(x=tmp.criteria, pattern='^[\'"]NA[\'"]$', perl=TRUE)]
      # ---> Get evaluations per cell per tag value that need to be evaluated.
      # Check there's at least one value to evaluate.
      if(length(tmp.criteria)>0){ to.keep <- sapply(X=tmp.criteria, FUN=function(criteria.term) seurat.obj@meta.data[, tmp.tag]==criteria.term); normal.evals <- TRUE} else normal.evals <- FALSE
      # Evaluate NA value if required.
      if(na.flag){
        na.evals <- is.na(seurat.obj@meta.data[, tmp.tag])
        if(!normal.evals) to.keep <- na.evals else to.keep <- cbind(to.keep, na.evals)
      }
      if(to.do=='discard'){
        to.keep <- !to.keep
        # In case many tag values were evaluated.
        if(!is.null(dim(to.keep))) to.keep <- apply(X=to.keep, MARGIN=1, FUN=all, na.rm=TRUE)
      }else{
        # In case many tag values were evaluated.
        if(!is.null(dim(to.keep))) to.keep <- apply(X=to.keep, MARGIN=1, FUN=any, na.rm=TRUE)
      }
      return(to.keep)
    })
    if(!is.null(dim(tags.evals))) tags.evals <- apply(X=tags.evals, MARGIN=1, FUN=all)
  }else{
    tags.evals <- NULL
  }
  # By features criteria.
  if(!is.null(feats.criteria)){
    feats.evals <- sapply(X=colnames(feats.criteria), FUN=function(tmp.feat){
      upper.thold <- feats.criteria['upper', tmp.feat]
      lower.thold <- feats.criteria['lower', tmp.feat]
      if(is.na(upper.thold)) upper.thold <- max(seurat.obj@meta.data[, tmp.feat])
      if(is.na(lower.thold)) lower.thold <- min(seurat.obj@meta.data[, tmp.feat])
      to.keep <- seurat.obj@meta.data[, tmp.feat]>=lower.thold & seurat.obj@meta.data[, tmp.feat]<=upper.thold
      return(to.keep)
    })
    if(!is.null(dim(feats.evals))) feats.evals <- apply(X=feats.evals, MARGIN=1, FUN=all)
  }else{
    feats.evals <- NULL
  }
  # ---> Assess all evaluation over cells.
  if(is.null(feats.evals)|is.null(tags.evals)){
    if(is.null(feats.evals)) all.evals <- tags.evals else all.evals <- feats.evals
  }else{
    all.evals <- cbind(feats.evals, tags.evals)
    all.evals <- apply(X=all.evals, MARGIN=1, FUN=all)
  }
  names(all.evals) <- rownames(seurat.obj@meta.data)
  # ---> Subset seurat object accordingly.
  cells.to.keep <- names(all.evals)[all.evals]
  seurat.obj <- subset(x=seurat.obj, cells=cells.to.keep)
  # ---> Return subsetted seurat object.
  return(seurat.obj)
}

# 19 ------------------------------------------------------------------->
# Name: Translate ID into name.
# Description:
# Given a set (vector) of ENSEMBL IDs, this function will return the same vector with their common gene names as vector names.
# ~~~* Requires the info from the data frame called feature.names (see code within seurat analysis script versions 1.4.2 or higher).
# Arguments ------------------------->
# ids - Vector with ENSEMBL IDs as values.
# Values ---------------------------->
# * Same vector with common gene names as vector names.
# Function:

translate.ids <- function(ids){
  names(ids) <- sapply(X=ids, FUN=function(tmp.id){
    tmp.names <- unique(feature.info$name[feature.info$ensembl==tmp.id])
    if(length(tmp.names)>1) feat.name <- paste0(tmp.names, collapse=';') else feat.name <- tmp.names
    return(feat.name)
  })
  return(ids)
}


### TO CHECK (EITHER BECAUSE THEY'RE COMING FROM A DIFFERENT SOURCE OR BECAUSE I HAVEN'T FINISHED IT)

# ? -------------------------------------------------------------------->
# Name: Draw general Euler diagram
# Description ----------------------->
# Given a list with NAMES items sets (where they can be unique or not), this function will output the appropriate Venn diagram.
# Arguments ------------------------->
# lib.info - Data frame storing clootypes info (MIGEC output) from a single library.
# depths - Integer vector. Minimum amount of samples a clonotype must be captured by to be marked as TRUE. It can be different depth values as log as they are less or equal the amount of samples in the library.
# get.clons - If requested, the output value will be a list storing 1) the clonotypes nt sequences assessed (those found in the library's samples) and 2) the matrix depicting the presence/absence and overlaps for each depth for all clonotypes.
draw.gen.venn <- function(sets.list){
  if(is.null(names(sets.list))) stop('Items sets must be named.')
  if(length(sets.list) > 5) stop('Sets number limits to 5.')
  tmp.cols <- sample.fav.colors(size=length(sets.list))
  venn.d <- NULL
  if(length(sets.list) == 1) venn.d <- draw.single.venn(area=length(unique(sets.list[[1]])), category=names(sets.list)[1], fill=tmp.cols)
  if(length(sets.list) == 2) venn.d <- draw.pairwise.venn(area1=length(unique(sets.list[[1]])), area2=length(unique(sets.list[[2]])), cross.area=length(calculate.overlap(sets.list)[[3]]), category=names(sets.list), fill=tmp.cols)
  if(length(sets.list) == 3){
    ovrlp.12 <- calculate.overlap(sets.list[[c(1,2)]])
    ovrlp.13 <- calculate.overlap(sets.list[[c(1,3)]])
    ovrlp.23 <- calculate.overlap(sets.list[[c(2,3)]])
    venn.d <- draw.triple.venn(area1=length(unique(sets.list[[1]])), area2=length(unique(sets.list[[2]])), area3=length(unique(sets.list[[3]])), n12=length(ovrlp.12[[3]]), n13=length(ovrlp.13[[3]]), n23=length(ovrlp.23[[3]]), n123=length(unique(Reduce(x=sets.list, f=intersect))), category=names(sets.list), fill=tmp.cols)
  }
  if(length(sets.list) == 4){
    ovrlp.12 <- calculate.overlap(sets.list[[c(1,2)]])
    ovrlp.13 <- calculate.overlap(sets.list[[c(1,3)]])
    ovrlp.14 <- calculate.overlap(sets.list[[c(1,4)]])
    ovrlp.23 <- calculate.overlap(sets.list[[c(2,3)]])
    ovrlp.24 <- calculate.overlap(sets.list[[c(2,4)]])
    ovrlp.34 <- calculate.overlap(sets.list[[c(3,4)]])
    venn.d <- draw.quad.venn(area1=length(unique(sets.list[[1]])), area2=length(unique(sets.list[[2]])), area3=length(unique(sets.list[[3]])), area4=length(unique(sets.list[[4]])), n12=length(ovrlp.12[[3]]), n13=length(ovrlp.13[[3]]), n14=length(ovrlp.14[[3]]), n23=length(ovrlp.23[[3]]), n24=length(ovrlp.24[[3]]), n34=length(ovrlp.34[[3]]), n123=length(unique(Reduce(x=sets.list[[c(1,2,3)]], f=intersect))), n124=length(unique(Reduce(x=sets.list[[c(1,2,4)]], f=intersect))), n134=length(unique(Reduce(x=sets.list[[c(1,3,4)]], f=intersect))), n1234=length(unique(Reduce(x=sets.list, f=intersect))), category=names(sets.list), fill=tmp.cols)
  }
  if(!is.null(venn.d)){return(venn.d)}else{return('There was an error.')}
}

#?
# how to fit a number of objects in a 2D dimention
fitgrid <- function(x, nCol = 'x', v = FALSE){
  if(is.numeric(x) && length(x) == 1) x <- 1:x
  vecL <- length(x)
  if(v) cat('Length:', vecL, '\n')
  if (is.null(x = nCol)) {
    if(v) cat('Adjust for nCol:', nCol, '\n')
    nCol <- 2
    if (length(vecL) == 1) nCol <- 1
    if (length(vecL) > 6) nCol <- 3
    if (length(vecL) > 9) nCol <- 4
    nRow <- floor(x = length(vecL) / nCol - 1e-5) + 1
  }else{
    if(vecL%in%1:2) nCol <- 1 else nCol <- round(sqrt(vecL)+.5)
    if(vecL==1) nRow <- 1 else nRow <- round(vecL/nCol)
    if(nRow*nCol < vecL) nRow <- nRow + 1
    if(!sqrt(vecL)%%1) nRow <- nCol <- sqrt(vecL)
  }
  if(v) cat('nRow:', nRow, '\nnCol', nCol, '\n')
  return(c(nRow, nCol, ifelse(vecL <= 2,  5, 4)))
}


# 11 ------------------------------------------------------------------->
# Name: Add feature tag
# Description:
# Given a set of features, this function reports for any cell if that feature is expressed or not according to a threshold value.
# Arguments ------------------------->
# To describe.
add_gene_tag <- function(
  lgenes,
  annot,
  mat,
  thresh = 0,
  tag = c('tag', '+', '-'),
  prefix = FALSE,
  v = FALSE
){
  if(length(thresh) != length(lgenes)){ # check if there is a threshold for each gene
    if(length(thresh) == 1) thresh <- rep(thresh, length(lgenes))
    if(!is.null(names(thresh))) thresh <- thresh[lgenes]
    thresh <- thresh[1:length(lgenes)]
    if(sum(is.na(thresh))) warning('No threshold specified for ', commas(lgenes[is.na(thresh)]), '; using 0\'s')
    thresh[is.na(thresh)] <- 0
  }
  names(thresh) <- lgenes
  newcolms <- cbindList(lapply(lgenes, function(thisgene){
    if(v) cat('----\nGene:', thisgene, '- threshold:', thresh[thisgene], '\n')
    if(prefix){
      tags <- paste0(tag, "_", casefold(thisgene, upper = TRUE))
    }else{
      tags <- paste0(tag[1], "_", casefold(thisgene, upper = TRUE))
      tags <- c(tags, paste0(casefold(thisgene, upper = TRUE), tag[-1]))
    }
    annot[, tags[1]] <- ifelse(as.vector(t(mat[thisgene, rownames(annot)] > thresh[thisgene])), tags[2], tags[3])
    tvar <- rownames(annot[annot[, tags[1]] == tags[3], ]) # check values
    tvar <- unlist(mat[thisgene, tvar])
    # if(v){ cat("Negative summary:\n"); print(summary(tvar)) }
    tvar <- rownames(annot[annot[, tags[1]] == tags[2], ]); tvar <- unlist(mat[thisgene, tvar])
    # if(v){ cat("Positive summary:\n"); print(summary(tvar)) }
    if(v){ cat('Proportions'); print(table(annot[, tags[1]])) }
    return(annot[, tags[1], drop = FALSE])
  }))
  newcolms
}
