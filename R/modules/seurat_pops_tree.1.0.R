############    ----------   SEURAT ANALYSIS    ---------    ############
############    ---------   Populations tree    ---------    ############

############    -------   Function description   --------    ############
# Name: Populations' tree.
# Dependencies:
#   Seurat, clustree, stringr, R_handy_functions
# Functions dependencies:
#   None
# Description:
# This functions takes a seurat object with clusters calculated for certain resolution values and outputs a clustering tree based on the library clusttree.
# Arguments ------------------------->
# seurat.obj - seurat object.
# reports.path - Absolute path to where result tree should be output.
# gen.norm.method - General normalization method data was normalized with.
# Value ----------------------------->
# Outputs reuslt to reports.path. Mock return.

############    --------   Program definitions   --------    ############
# ---> This program's version-related details:
# Version: 1
# Version updates:
#   First version.
# Subversion: 0
#   First subversion.
# Associations:
# None so far.

out.pops.tree <- function(seurat.obj, reports.path, gen.norm.method='RNA'){
  cat('\n\n')
  cat('############    ---------   Populations tree    ---------    ############\n')

  ### ------------------------- Preprocessing ------------------------- ###
  # Define tags preffix.
  clusts.tag.pffx <- paste0(ifelse(test=gen.norm.method=='SCTransform', yes='SCT', no='RNA'), '_snn_res.')
  # Get resolutions' columns.
  clusts.cols <- grep(x=colnames(seurat.obj@meta.data), pattern=paste0('^', clusts.tag.pffx, '\\d+[\\.\\d+]*$'), value=TRUE, perl=TRUE)
  # Define clusters' info with suitable names.
  clusts.data <- seurat.obj@meta.data[, clusts.cols]
  colnames(clusts.data) <- str_replace(string=colnames(clusts.data), pattern=clusts.tag.pffx, replacement='Resolution')

  ### -------------------------- Tree output -------------------------- ###

  tmp.file.name <- paste0(reports.path, '/PopulationsByResolutionAndTheirRels_ClustTree.pdf')
  pdf(file=tmp.file.name, height=10, width=10)
  print(clustree(x=clusts.data, prefix='Resolution', edge_width=1))
  dev.off()

  # Mock return.
  return(NULL)
}
