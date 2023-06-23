############    ----------   SEURAT ANALYSIS    ---------    ############
############    ------   Transform current tags   -------    ############

# Version: 0
# Version updates:
#   Not a stable version yet.
# Subversion: 3
# Subversion updates:
# ---> Subversion 3
# Sub-subversion updates:
#   Metadata columns may now be removed before being recalculated (this may be necessary when previous annotations are not right for certain).


# Dependencies:
#   Seurat, data.table and stringr
# Functions dependencies:
#   None.
# Description ----------------------->
#   The idea and main code body for this script comes from: /home/vfajardo/scripts/seurat_analysis/tag_specific_analysis.2.8.R (and probably, comes originally from previous versions of that script).
#   Given the arguments described by the columns in the table provided that is listing rules to create new tags, this script will perform the task of creating new tags by taking advantage of the information distributed in multiple tags to create new tags (e.g., to keep information only for a subset of cells with a specific set of values from a given tag). In order to perform this task, the script works through the execution of two different functions as described below.
#   Of note, this script cheks that all tags provided in the 'merge' column of the rules file do exist.


### --------------------------- Functions --------------------------- ###

# 1 -------------------------------------------------------------------->
# Name: Add new tags.

# Description ----------------------->
# This is the function that should be called straight from the main program.
# Based on the set of rules provided in the data frame 'new.tags.rules' and the execution of the function 'get.new.tag', this function obtains the new tags and adds them straight to the metadata of the input seurat object, returning such a modified object. This is achieved in three main steps:
#   1. Checking of the rules object to make sure it's a valid one.
#   2. Definition of the new tags (one at a time) through the execution of the function 'get.new.tag'.
#   3. Addition of the new information to the input object's metadata.

# Arguments ------------------------->
# seurat.obj - Input seurat object whose metadata should contain the tags to be used to create the new tags as defined in the rules object.
# new.tags.rules - Data frame, whose specificities are defined next.
#   Describes comprehensively the rules to merge two or more different tags already defined in the input seurat object, thereby creating a new tag. The data frame should have next columns (in any order):
#     1) tag.names: For each row, describes the name of the new tag to be produced according to the specifications given in the rest of the columns. This is the only case where NAs are discarded before concatenating all values.
#     2) merge: For each row, sequence of already defined tags to be -somehow- combined, which should be seprated by and only by semicolon (e.g., 'virus.tag;donor.id.tag;batch.tag'). If rule equals 'join' or 'add' (see below), then the order given will be the one to use for the output.
#     3) rule: For each row, one of three options as follows:
#         a) join: Tag values are joined/concatenated (in the order provided) and separated by the pattern provided.
#         b) remove: Final tag values for the new tag will be produced by the appending of the columns specified in merge (separated by a dash and given in the order provided in that column), such that the tag values provided in 'pattern' will be discarded and, therefore, the final value for the cells having them will be set to NA. The tag values defined in 'pattern' and to be removed (i.e., to be set to NA) must be defined as valid values for any of the tags defined in 'merge', however, if any these 'pattern' values are defined for multiple tags defined in 'merge' (i.e., they're redundant), the program will fail outputting an error.
#           For example, if you provided as merge argument 'virus.tag;donor.id.tag;batch.tag' and you set pattern to 'Flu;P07', assuming that they're valid values for 'virus.tag' and 'donor.id.tag', respectively, then the final values for the new column will be composed by the appending of 'virus.tag-donor.id-batch.tag' (notice that '-' is 'joining' the classes) and any cell 1) with value 'Flu' and/or 'P07' for the tags 'virus.tag' and 'donor.id.tag', respectively or 2) with value set originally to NA for any of the original tags will be marked by NA in the final column/new tag.
#         c) keep: Same action as for 'remove', but insted of setting NA as final value to the cells being classified as the values in the 'pattern' argument, it will set NA to all values not being classified as such.
#         d) add.end: Concatenates 'pattern' at the end of the concatenation of the 'merge' columns provided (these ones always joined by '-'). That is, same as 'join', but will add the value given to pattern at the end of the result and always uses '-' to join everything.
#         e) add.beg: Same as add.end, but adding 'pattern' argument at the beggining instead.

# Value ----------------------------->
# Updated seurat object, that is, including the new tags gotten through the merging of the others based on the rules provided.

# Function -------------------------->

add.new.tags <- function(seurat.obj, new.tags.rules, rem.tags=FALSE){

  ### ---------------------- Data preprocessing ----------------------- ###

  # ---> Preflights.
  if(!all(c('tag.name', 'merge', 'rule', 'pattern') %in% colnames(new.tags.rules))) stop(paste0('File describing rules to create new tags was provided (path below), but it does not have the format (necessary columns) required.\n'))
  tmp.tags <- unique(unlist(str_split(string=new.tags.rules$merge, pattern=';', simplify=FALSE)))
  if(!all(tmp.tags %in% colnames(seurat.obj@meta.data))){
    tmp.tags <- tmp.tags[!tmp.tags %in% colnames(seurat.obj@meta.data)]
    stop(paste0('File describing rules to create new tags was provided (path below), but any of the tags\' names provided in the column \'merge\' has not been previously defined in the seurat object as listed below:\n', paste0(tmp.tags, collapse='\n'), '\n'))
  }
  if(any(new.tags.rules$tag.name %in% colnames(seurat.obj@meta.data))){
    if(rem.tags){
      tmp.tags <- new.tags.rules$tag.name[new.tags.rules$tag.name %in% colnames(seurat.obj@meta.data)]
      for(tmp.tag in tmp.tags){
        seurat.obj@meta.data[, tmp.tag] <- NULL
      }
    }else{
      stop('New tag names should not be previosly defined in the suerat object.\n')
    }
  }

  ### ------------------------- Main program -------------------------- ###

  # ---> Get new tags.
  # New tag definition.
  new.tags <- lapply(X=1:nrow(new.tags.rules), FUN=function(idx){
    tag.name <- new.tags.rules[idx, 'tag.name']; merge <- new.tags.rules[idx, 'merge']; rule <- new.tags.rules[idx, 'rule']; pattern <- new.tags.rules[idx, 'pattern']
    new.tag <- get.new.tag(tag.name=tag.name, tags.to.merge=merge, rule=rule, pattern=pattern)
    return(new.tag)
  })
  names(new.tags) <- new.tags.rules$tag.name
  new.tags <- rbindlist(l=new.tags, use.names=TRUE, idcol='tag.name')
  new.tags <- spread(data=new.tags, key=tag.name, val=new.tag)
  new.tags <- as.data.frame(new.tags)
  row.names(new.tags) <- new.tags$cell; new.tags <- new.tags[, new.tags.rules$tag.name]
  new.tags <- new.tags[Cells(seurat.obj), ]
  # Add new tags to previous seurat object.
  meta.data <- merge(seurat.obj@meta.data, new.tags, by='row.names', sort=FALSE)
  row.names(meta.data) <- meta.data$Row.names
  meta.data$Row.names <- NULL
  if(!all(row.names(meta.data)==Cells(seurat.obj))) stop('Unexpected error 1. See program.\n')
  seurat.obj@meta.data <- meta.data
  rm(meta.data)
  rm(new.tags)
  # Return modified seurat object.
  return(seurat.obj)
}


# 2 -------------------------------------------------------------------->
# Name: Get new tag.

# Description ----------------------->
# This function is meant to work for every row for every row in the data.frame 'new.tags.rules'. It does retrieve the information from the original seurat object and based on the rules defined in a single row of the data frame (i.e., the rules meant to create one single new tag), it defines that new tag.

# Arguments ------------------------->
# All taken from the file provided as new.tags.file. Seurat object is taken from the main environment. The rest of the arguments are the ones defined in the data frame 'new.tags.rules' (or its equivalent) per row (see documentation above if needed).
# tag.name
# tags.to.merge - equvalent to 'merge'
# rule
# pattern

# Value ----------------------------->
# Data table listing the new tag's value for each cell.

# Function -------------------------->

get.new.tag <- function(tag.name, tags.to.merge, rule, pattern){

  ### ---------------------- Data preprocessing ----------------------- ###

  if(!any(rule %in% c('join', 'remove', 'keep', 'add.end', 'add.beg'))) stop(paste0('Invalid rule was provided in the file describing rules to create new tags.\nInvalid rule: ', rule, '\n\n'))
  # Get tags' names.
  tags.to.merge <- unlist(str_split(string=tags.to.merge, pattern=';'))
  if(!length(unique(tags.to.merge)) == length(tags.to.merge)) stop(paste0('For new tag named as ', tag.name, ', please provide unique tag names to merge.\n'))


  ### ------------------------- Main program -------------------------- ###

  # Create new tag accordingly.
  if(rule=='join'){
    # ---> Join.
    to.output <- tidyr::unite(data=seurat.obj[[tags.to.merge]], col=new.tag, sep=pattern, remove=TRUE)
  }else{
    if(rule %in% c('remove', 'keep')){
    # ---> Remove or keep.
      # Check tag values provided as pattern are valid in the sense that they are not redundant and are coming from a valid tag (notice that redundant tag values are not supported). Then, order pattern values according to tag origin.
      pattern <- unique(unlist(str_split(string=pattern, pattern=';')))
      tag.vals <- lapply(X=pattern, FUN=function(tmp.val){
        tmp.data <- seurat.obj@meta.data[, tags.to.merge]==tmp.val & !is.na(seurat.obj@meta.data[, tags.to.merge])
        tmp.data <- colSums(tmp.data)
        tmp.data <- names(tmp.data)[tmp.data>1]
        tmp.check <- length(tmp.data) > 1
        if(tmp.check){ tmp.err <- paste0('Tag value \'', tmp.val, '\' was found to be defined for multiple tags, mainly: ', paste0(tmp.data, collapse=', '), '. Such a redundancy cannot be resolved in the current version of the program. Please develope a new version to deal with it. Alternatively, try and rearrange the rules so that tag values are not redundant.\n'); stop(tmp.err) }
        tmp.check <- length(tmp.data) < 1
        if(tmp.check){ tmp.err <- paste0('Tag value \'', tmp.val, '\' is not valid for any of the tags defined for ', tag.name, ', the tag to be created. Please check it is a valid tag value.\n'); stop(tmp.err) }
        return(tmp.data)
      })
      tag.vals <- unlist(tag.vals)
      names(tag.vals) <- pattern
      tag.names <- unique(tag.vals); tag.vals <- lapply(X=tag.names, FUN=function(tmp.name) names(tag.vals)[tag.vals==tmp.name]); names(tag.vals) <- tag.names
      # Pattern values are then replaced by NA for any cell classified (for 'remove') or not classified (for 'keep') with those values.
      cells.to.modify <- lapply(X=names(tag.vals), FUN=function(tmp.tag){
        pattern.vals <- tag.vals[[tmp.tag]]
        to.output <- seurat.obj@meta.data[, tmp.tag] %in% pattern.vals
        if(rule=='keep') to.output <- !to.output
        to.output <- data.table(cell=Cells(seurat.obj), keep=to.output)
        return(to.output)
      }); names(cells.to.modify) <- names(tag.vals)
      cells.to.modify <- rbindlist(l=cells.to.modify, use.names=TRUE, idcol='tag')
      cells.to.modify <- spread(data=cells.to.modify, key=tag, value=keep)
      cells.to.modify <- apply(X=as.matrix(cells.to.modify, rownames='cell'), MARGIN=1, FUN=if(rule=='keep') any else all)
      # Proceed merging columns.
      to.output <- tidyr::unite(data=seurat.obj[[tags.to.merge]], col=new.tag, sep='-', remove=TRUE)
      # Set NAs according to the rules provided by 'remove' or 'keep'
      to.output[cells.to.modify[row.names(to.output)], 'new.tag'] <- NA
    }else{
    # ---> Add.
      to.output <- tidyr::unite(data=seurat.obj[[tags.to.merge]], col=new.tag, sep='-', remove=TRUE)
      if(rule=='add.end'){
        to.output$new.tag <- paste0(to.output$new.tag, pattern)
      }else{
        to.output$new.tag <- paste0(pattern, to.output$new.tag)
      }
    }
  }
  # Define in a established format and output.
  to.output$cell <- row.names(to.output)
  to.output <- as.data.table(to.output[, c('cell', 'new.tag')])
  # to.output[, tag.name:=tag.name]
  return(to.output)
}
