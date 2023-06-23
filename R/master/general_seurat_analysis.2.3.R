############    --------   scRNA-seq ANALYSIS    --------    ############
############    ------------   USING SEURAT   -----------    ############

# By Vicente Fajardo
# Some code adapated from Ciro's.

# Whole version: 2.3.8
# --> Version: 2
# Version updates:
#   Check version 2.2.6
# --> Subversion: 3
# Sub-subversion updates:
#   In this new version, we introduce the option to use Harmony to deal with batch effects.
# --> Sub-subversion: 1
#   * When a seurat object is provided as input, we clean the metadata of the previous object effectively. For that, we keep only the  following tags: 'orig.ident', 'nCount_RNA' and 'nFeature_RNA', i.e., the ones that a recently-created seurat object contains.
#   * Now the user can let know the program that it's ok having more batches defined in the aggr table provided than the actual number of batches defined for the seurat object. This implies a new verison (1.2) for the module 'annotate_seurat_obj'.
#   * Minor bug fix: Marker feature plots process is done linearly (new version, 1.2, of module 'seurat_clusters_definition_and_qc').
#   * Minor bug fix: We confirm that we've got valid values for ID tags for which we'll calculate cell counts (e.g., donor ID).
# --> Sub-subversion: 2
#   * We open the possibility to regress a whole bunch of genes when normalizing the data as a single signature instead of individual variables (new version, v2.0, of the module 'seurat_log_normalization').
# --> Sub-subversion: 3
#   * There's now the possibility to exclude B-cell receptor (BCR) genes (similar to the option for removing TCR genes) from the whole analysis (new version, v1.3, of the module 'seurat_qc_analysis').
# --> Sub-subversion: 4
#   * There's now the possibility to provide a set of barcodes to subset an input seurat object.
# --> Sub-subversion: 5
#   * Output violin plots depicting distribution of QC features across clusters with a better resolution (new version, v1.3, of the module 'seurat_clusters_definition_and_qc').
# --> Sub-subversion: 6
#   * There's now the possibility to provide a set of barcodes to subset, not only an input seurat object, but for any kind of input. So just as it worked to presubset a seurat object file, the user may also provide a csv file listing a set of barcodes to be kept for analysis taken from the raw cellranger aggr/count data (see option "--TagsCriteria").
# --> Sub-subversion: 7
#   * Specific adaptations for the script to work on a slurm server as well as in a torque cluster (the latter was the only option in previous versions).
# --> Sub-subversion: 8
#   * New QC analysis script version (1.4). Updates specified within the script.
# --> Sub-subversion: 8
#   * New initial dim reduction script version (1.1). Updates specified within the script.


cat('############    --------   scRNA-seq ANALYSIS    --------    ############\n')
cat('############    ------------   USING SEURAT   -----------    ############\n')
cat('By Vicente Fajardo\n\n')

cat('\n\n')
### --------------------------- Libraries --------------------------- ###
cat('### --------------------------- Libraries --------------------------- ###\n')
cat('Importing libraries...\n\n')
library(dplyr, warn.conflicts=FALSE)
library(Seurat)
library(ggplot2)
library(future)
library(optparse)
library(data.table)
library(ggrepel)
library(parallel)
library(stringr)
library(clustree)
library(tidyr)
library(harmony)
source('/home/vfajardo/scripts/functions/R_handy_functions.0.3.R')
source('/home/vfajardo/scripts/functions/R_visualization_functions.1.5.R')
# Modules and important functions.
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/annotate_seurat_obj.1.2.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/add_new_tags.0.3.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/seurat_qc_analysis.1.4.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/seurat_log_normalization.2.0.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/seurat_initial_dim_reduction.1.1.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/seurat_clustering_and_non-linear_dim_reduction.1.0.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/seurat_populations_definitions.1.2.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/seurat_clusters_definition_and_qc.1.3.R')
source('/home/vfajardo/scripts/seurat_analysis/seurat_analysis_modules/seurat_pops_tree.1.0.R')

cat('Libraries imported!\n\n')

cat('\n\n')
### ------------------------ Parallelisation ------------------------ ###
cat('### ------------------------ Parallelisation ------------------------ ###\n')
# Calculate total number of processors we've got.
total.workers <- as.integer(system(command='echo $PBS_NUM_PPN', intern=TRUE))
if(is.na(total.workers)) total.workers <- as.integer(system(command='echo $SLURM_NTASKS', intern=TRUE))
# Plan regarding number of processors.
plan("multisession")

cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# Declaring arguments to parse from command line ----------------------->
option.list <- list(
  ### ------------------> Generalities.
  make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Absolute path to directory to save analyses reports"),
  make_option(opt_str="--PrjName", type="character", default=NULL, dest="prj.name", help="Project name"),
  make_option(opt_str="--DataFile", type="character", default=NULL, dest="data.file", help="Absolute path to file storing a feature by counts matrix from."),
  make_option(opt_str="--InputType", type="character", default="matrix", dest="input.type", help="Character value indicating what kind of input should be considered, either 'counts' for count matrix or 'seurat' for seurat objetct."),
  make_option(opt_str="--FeatureID", type="character", default="name", dest="feature.id", help="Character defining the feature ID type to take from 10X cellranger output (either from features.tsv file or from the h5 matrix). Default is 'name' (features.tsv column 2), though possible values are 'name' and 'ensembl', the last one corresponding to ENSEMBL ID. This option applies only if a cellranger output path is provided.\n"),
  make_option(opt_str="--Raw10x", type="character", default=NULL, dest="tenx.raw.data", help="Absolute path to 10x raw data path (as a directory instead of an hf5 matrix file). Necessary when \"ensembl\" is picked as the feature ID option.\n"),
  make_option(opt_str="--IsFivePrime", type="logical", dest="is.five.prime", default=FALSE, help="Add flag to indicate that the data comes from a five prime experiment, in which case, antigen receptor (either T-cell receptor or B-cell receptor) genes will be excluded from downstream analyses."),
  make_option(opt_str="--SignificancePlot", type="integer", dest="no.top.sig.feats", default=16, help="Number of top significant features depicted in Ciro's significance plots."),
  make_option(opt_str="--MarkersFile", type="character", dest="markers.file", default=NULL, help="Some markers (mainly related to T cells functions) will be assessed, but you can provide more markers for this purpose as an RData file."),
  make_option(opt_str="--RAM", type="integer", dest="avl.ram", default=NULL, help="Available RAM for the job session."),

  ### ------------------> For Annotations.
  make_option(opt_str="--DoAnnotate", type="logical", default=FALSE, dest="do.annotate", help="Lane tag ID for seurat object."),
  make_option(opt_str="--AggrTable", type="character", default=NULL, dest="aggr.table.file", help="Absolute to fully-annotated aggregation table. Same format as a cellranger aggr table, but with at least one of next tags: hto.tag, lane.tag,  chromium.batch.tag or seq.batch.tag as annotations for each sample. Then, each row describes the info for each 10x lane (lane.tag) per pairs of chromium platform batch and sequencing batch (batch.tag), and hto.tag column lists the absolute paths to the ClassificationPerCell.csv file output by the HTO deconvolution (or similar) pipeline."),
  make_option(opt_str="--DonorsMetaData", type="character", default=NULL, dest="donors.data.file", help="Absolute path to donors metadata file. Set as NULL if there's none."),
  make_option(opt_str="--ExcededBatches", type="logical", default=FALSE, dest="exceded.batches", help="Logical, indicates whether it's ok having more batches defined in the aggr table provided than the actual number of batches defined for the seurat object. If FALSE, the pipeline will get stalled when there are less batches defined (in the seurat object) than the ones expected according to the aggr table."),
  make_option(opt_str="--MergeBy", type="character", default=NULL, dest="donors.merge.by", help="Character, column name in both seurat object meta data and donors meta data in order to merge them and add further annotations. Value considered only if there are annotations and donors metadata provided."),
  make_option(opt_str="--LaneID", type="character", default='lane.id.tag', dest="lane.tag.id", help="Lane tag ID for seurat object. Single tag IDs should not contain semicolons if it is desired to add multiple tags according to lane."),
  make_option(opt_str="--MultAnnPerLane", type="logical", default=FALSE, dest="mult.ann.per.lane", help="Logical, indicates wether the lane tag stands for multiple annotations. If TRUE, they should be separated by semicolon within each aggregation table lane index and there should be the same amount of values per library."),
  make_option(opt_str="--ChromBatchID", type="character", default="chrom.batch.tag", dest="chrom.batch.tag.id", help="10X Chromium Platform Tag ID for seurat object."),
  make_option(opt_str="--SeqBatchID", type="character", default="seq.batch.tag", dest="seq.batch.tag.id", help="Sequencing batch Tag ID for seurat object."),
  make_option(opt_str="--HTOID", type="character", default="hto.tag", dest="hto.tag.id", help="HTO Tag ID for seurat object."),
  make_option(opt_str="--OverlayID", type="character", default="overlay.id.tag", dest="overlay.tag.id", help="Overlay (between chromium batch and HTO tags) Tag ID for seurat object."),
  make_option(opt_str="--NewTags", type="character", default=NULL, dest="new.tags.file", help="Absolute path to file describing the rules to add new tags based on the combination of the ones already defined. See more details in the description above.\n"),

  ### ------------------> Subsetting options.
  make_option(opt_str="--DoPreSubset", type="logical", default=FALSE, dest="do.presubset", help="Logical, only valid when a seurat object is input, either presubsetting should be applied or not."),
  make_option(opt_str="--PreSubsetCriteria", type="character", default=NULL, dest="presubset.criteria.file", help="Character, absolute path to file indicating the criteria used to apply presubsetting. Only valid when a seurat object is input. Two different subset strategies are valid depending on the columns defined in the file. If the file has a single column (named \"barcode\"), then the barcodes listed there will be taken as the subset to be kept for downstream analyses. Otherwise, the file has to have the same format as for Tag and Features criteria and the same subsetting strategy is applied\n."),
  make_option(opt_str="--DoSubset", type="logical", default=FALSE, dest="do.subset", help="Logical, flag indicating wether there should be any kind of filtering regarding the seurat object annotations."),
  make_option(opt_str="--TagsCriteria", type="character", default=NULL, dest="tags.criteria.file", help="Character, absolute path to file indicating the criteria used to apply presubsetting. Two different subset strategies are valid depending on the columns defined in the file. If the file has a single column (named \"barcode\"), then the barcodes listed there will be taken as the subset to be kept for downstream analyses. Otherwise, the file has to have the same format as for Tag and Features criteria and the same subsetting strategy is applied\n."),
  make_option(opt_str="--FeatsCriteria", type="character", default=NULL, dest="feats.criteria.file", help="Character, absolute path to table storing the criteria to filter cells based on for tags required (format: criteria by tag, column name should match that of tag defined in seurat object; criteria should be separated in two different rows, lower and upper threshold respectively)."),

  ### ------------------> For QCs.
  make_option(opt_str="--minCounts", type="integer", default=300, dest="min.nCounts", help="Integer indicating the lower threshold for counts per cell"),
  make_option(opt_str="--maxCounts", type="integer", dest="max.nCounts", help="Integer indicating the upper threshold for counts per cell"),
  make_option(opt_str="--minFeatures", type="integer", default=200, dest="min.nFeatures", help="Integer indicating the lower threshold for number of genes per cell"),
  make_option(opt_str="--maxFeatures", type="integer", default=2500, dest="max.nFeatures", help="Integer indicating the upper threshold for number of genes per cell"),
  make_option(opt_str="--maxMP", type="numeric", default=0.05, dest="max.MP", help="Number indicating the maximum mitochondrial percentage to allow per cell"),
  make_option(opt_str="--FilterOut", type="logical", dest="filter.out", default=TRUE, help="Mark if, given the quality control measurement, likely low quality cells should be filtered out."),

  ### ------------------> For normalization.
  make_option(opt_str="--FVFsMethod", type="character", default="vst", dest="FVFs.method", help="Character string indicating the method to use to find variable features"),
  make_option(opt_str="--FeatsForDSA", type="integer", default=30, dest="feats.for.dsa", help="Either amount of variable features to consider or the percentage of standardized variance to take into consideration to pick a number of variable features for downstream analysis (mainly, as the number of variable features taken for PCA analysis). If the number is smaller than 100, this will be taken as a percentage and the other way around."),
  make_option(opt_str="--MeanCutoff", type="numeric", dest="mean.cutoff", default=0.1, help="Mean cutoff to exclude vairable features."),
  make_option(opt_str="--PropCutoff", type="numeric", dest="prop.cutoff", default=0.001, help="Proportion of cells cutoff to exclude variable features."),
  make_option(opt_str="--GenNormMethod", type="character", default="LogTransform", dest="gen.norm.method", help="Character string indicating the method to use for general data normalization. Possible values (other than these, program will stop): LogTransform (Log-normalization for cell-wise normalization proceeded by other steps such as gene normalization) and SCTransform (normalization and variance stabilization through Regularized Negative Bonimial Regression per gene and use of pearson residuals, see Hafemeister, C., & Satija, R. (2019). Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biology, 20(1), 1-15.)."),
  make_option(opt_str="--NormMethod", type="character", default="LogNormalize", dest="norm.method", help="Character string indicating the method to use for data normalization per cell. If SCTransform selected as the general normalization method, this will be ognored."),
  make_option(opt_str="--VarsToRegress", type="character", dest="vars.to.regress", default="c('nCount_RNA', 'percent.mt')", help="Character, either 1) an absolute path to a csv file (single column, with header 'features') listing a ser of variables to regress during gene normalization (scaling applied with the function ScaleData from Seurat, passed to the argument named the same way) or a character to be evaluated as a charcater vector using parse/eval.\n"),
  make_option(opt_str="--RegressCC", type="character", dest="regress.cc", default=NULL, help="Character, indicates whether cell cycle status information should be regressed while doing cell-wise normalization of expression values (i.e., further values to be regressed while running ScaleData). Possible values are: 1) 'human' for human signatures; 2) 'mouse' for the homologs from the human cell cycle signature; 3) NULL not to regress cell cycle signature info at all. Cell cycle signatures are taken directly from the Seurat package (see vignette 'Cell-Cycle Scoring and Regression').\n"),

  ### ------------------> For Dimensionality reduction and clustering.
  make_option(opt_str="--PCs", type="character", default="30", dest="nos.PCs", help="Number of PCs to take into consideration for downstream analyses after normalization. Many values may be provided as far as they are provided as an R vector."),
  make_option(opt_str="--PCsToComp", type="integer", default=40, dest="pcs.to.comp", help="No matter how many PCs will be taken into account, PCs to compute."),
  make_option(opt_str="--ForHarmony", type="character", default=NULL, dest="do.harmony", help="Character. If different than NULL, it indicates that Harmony batch-effect removal should be applied. Then, this may be evaluated as a character vector of the tags listing the batches to be removed through Harmony.\n"),
  make_option(opt_str="--Resolution", type="character", default="c(0.2, 0.4, 0.6)", dest="resolutions", help="Resolution argument to find clusters. Many values may be provided as far as they are provided as an R vector."),
  make_option(opt_str="--DimReduction", type="character", dest="dim.reduction.methods", default="c('umap', 'tsne')", help="String to be evaluated as a character vector indicating the dimensionality reduction methods that should be applied over the data other than PCA (which will run tried all the same)."),
  make_option(opt_str="--MaxCellsNo", type="integer", dest="max.cells.no", default=5000, help="Maximum number of cells to depict in plots.."),

  ### ------------------> For Differential expression analysis.
  make_option(opt_str="--DEA", type="logical", dest="do.dea", default=TRUE, help="Logical indicating if pairwise comparisons (DEA) should be applied and reported or not."),
  make_option(opt_str="--MaxCellsNo4DEA", type="integer", dest="dea.max.cells", default=200000, help="Maximum number of cells to use for differential gene expressio analysis.")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);

# Usual seurat analysis ------------------------------------------------>
# Moving options to their own variables

# ---> Generalities
reports.path <- opt$reports.path
prj.name <- opt$prj.name
data.file <- opt$data.file
# Input type. Make sure it equals seurat or counts.
input.type <- opt$input.type
if(!(input.type=='seurat' || input.type=='matrix')) stop('Not valid input type provided.')
feature.id <- opt$feature.id
# Check feature ID definition is correct.
if(feature.id!='name' & feature.id!='ensembl') stop('Feature ID not properly defined.\n')
# Check a 10x raw data directory is provided whenever "ensembl" is picked as the feature.id option.
tenx.raw.data <- opt$tenx.raw.data
if(is.null(tenx.raw.data)){
  if(feature.id=='ensembl') stop('Feature ID was defined as ensembl with no 10x raw data directory provided. Notice that even if it\'s the same value as the data.file argument, this should be specified anyway.\n') else warning('10x raw data directory was provided even when not necessary. Did you mean to use "ensembl" as feature.id argument?\n')
}
is.five.prime <- opt$is.five.prime
no.top.sig.feats <- opt$no.top.sig.feats
markers.file <- opt$markers.file
avl.ram <- opt$avl.ram
avl.ram <- avl.ram-(avl.ram*0.1)
# Plan regarding RAM.
options(future.globals.maxSize=avl.ram*1000*1024^2)

# ---> For annotations
do.annotate <- opt$do.annotate
aggr.table.file <- opt$aggr.table.file
exceded.batches <- opt$exceded.batches
donors.data.file <- opt$donors.data.file
donors.merge.by <- opt$donors.merge.by
lane.tag.id <- opt$lane.tag.id
mult.ann.per.lane <- opt$mult.ann.per.lane
chrom.batch.tag.id <- opt$chrom.batch.tag.id
seq.batch.tag.id <- opt$seq.batch.tag.id
hto.tag.id <- opt$hto.tag.id
overlay.tag.id <- opt$overlay.tag.id
new.tags.file <- opt$new.tags.file

# ---> Subsetting options
# Subsetting criteria in case a seurat object is input.
do.presubset <- opt$do.presubset
presubset.criteria.file <- opt$presubset.criteria.file
if(input.type!='seurat' & do.presubset){tmp.warning <- paste0('Pre-subsetting was set to be TRUE though there is no seurat object as input for the program.\n'); warning(tmp.warning)}
if(input.type=='seurat' & do.presubset) if(is.null(presubset.criteria.file) | !file.exists(presubset.criteria.file)) stop(paste0('No appropriate subset criteria file input.\n'))
# Main subsetting criteria.
do.subset <- opt$do.subset
feats.criteria.file <- opt$feats.criteria.file
tags.criteria.file <- opt$tags.criteria.file

# ---> For QCs
min.nCounts <- opt$min.nCounts
max.nCounts <- opt$max.nCounts
min.nFeatures <- opt$min.nFeatures
max.nFeatures <- opt$max.nFeatures
max.MP <- opt$max.MP
filter.out <- opt$filter.out

# ---> For normalization
FVFs.method <- opt$FVFs.method
feats.for.dsa <- opt$feats.for.dsa
# This number, as is suppoused to reduce the number of dimensions for PCA analysis, can't be greater than 2000.
if(feats.for.dsa>2000) stop('The number of variable features to capture cannot be greater than 2000.\n')
mean.cutoff <- opt$mean.cutoff
prop.cutoff <- opt$prop.cutoff
gen.norm.method <- opt$gen.norm.method
norm.method <- opt$norm.method
vars.to.regress <- opt$vars.to.regress
# Either a file or a charcater vector of length 1 to be evaluated as a character vector.
if(!file.exists(vars.to.regress)) vars.to.regress <- eval(expr=parse(text=vars.to.regress))
regress.cc <- opt$regress.cc

# ---> For dimensionality reduction and clustering
nos.PCs <- opt$nos.PCs; nos.PCs <- eval(expr=parse(text=nos.PCs))
pcs.to.comp <- opt$pcs.to.comp
do.harmony <- opt$do.harmony
if(!is.null(do.harmony)) do.harmony <- eval(expr=parse(text=do.harmony))
resolutions <- opt$resolutions; resolutions <- eval(expr=parse(text=resolutions))
# Dimensionality reduction methods to be applied.
dim.reduction.methods <- opt$dim.reduction.methods
dim.reduction.methods <- eval(expr=parse(text=dim.reduction.methods))
dim.reduction.methods <- unique(dim.reduction.methods[dim.reduction.methods=='umap'|dim.reduction.methods=='tsne'])
max.cells.no <- opt$max.cells.no

# ---> For Differential expression analysis.
do.dea <- opt$do.dea
dea.max.cells <- opt$dea.max.cells

# Common variables ----------------------------------------------------->
# --- QC features' thresholds
feats.to.assess <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
feats.to.assess.names <- c('Genes per cell', 'Seq. depth per cell', 'Mitochondrial genes pct.')
QC.tholds <- list(c(min.nFeatures, max.nFeatures), c(min.nCounts, max.nCounts), c(0, max.MP))
# Briefly, set thresholds as lower or upper boundaries.
QC.tholds <- lapply(X=QC.tholds, FUN=function(tholds){names(tholds) <- c('lower.thold', 'upper.thold'); return(tholds)})
names(QC.tholds) <- feats.to.assess.names
names(feats.to.assess) <- feats.to.assess.names
# --- Related to cell count summary.
filt.stages <- c(pre.subset='Pre-subset', pre.filter='Pre-filtering', post.filter='Post-filtering')

# Presenting parameters ------------------------------------------------>
cat('General reports path:', reports.path, '\n')
cat('Abolute path to count matrix file:', data.file, '\n')
cat('Project name:', prj.name, '\n')
cat('Filtering arguments:\n')
cat('Minimum and maximum thresholds for counts number per cell:', min.nCounts, 'and', max.nCounts, '\n')
cat('Minimum and maximum thresholds for genes number per cell:', min.nFeatures, 'and', max.nFeatures, '\n')
cat('Maximum mitochondrial genes percentage per cell:', max.MP, '\n')
cat('Normalization method:', norm.method, '\n')
cat('Variable features selection method:', FVFs.method, '\n')
cat('Number/Percentage of most variable genes to be considered', feats.for.dsa, '\n')
cat('Resolution value(s) for clustering analysis:', paste(resolutions, collapse=', '), '\n')
cat('Values of PCs to be taken into account for dimensionality reduction:', paste(nos.PCs, collapse=', '), '\n')
cat('Should pairwise comparisons among clusters (DEA) be applied?', do.dea, '\n\n')

### --------------------------- Functions --------------------------- ###
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

get.file.name <- function(gen.path, file.desc, extension, added.args=NULL, resolution=NULL, no.PCs=NULL){
  # Temporal before adding extra arguments.
  tmp.file.name <- paste0(gen.path, '/', file.desc)
  # Specify PCs or resolution if required.
  if(!is.null(no.PCs)) tmp.file.name <- paste0(tmp.file.name, '_NoPCs_', no.PCs)
  if(!is.null(resolution)) tmp.file.name <- paste0(tmp.file.name, '_Resolution_', resolution)
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


cat('\n\n')
### ------------------------- Data Loading ------------------------- ###
cat('### ------------------------- Data Loading ------------------------- ###\n')
# ---> Seurat object.
# Read count matrix or seurat object according to input.type varible value.
if(input.type=='matrix'){ # When counts matrix provided.
  data <- read.10X.data(file=data.file, feature.id=feature.id)# Create seurat object.
  cat('Data file imported. If standard errors, you should check the file before continuing with this program.\n\n')
  seurat.obj <- CreateSeuratObject(counts=data, project=prj.name, min.cells=0)
  rm(data)
  cat('Seurat object has been created!\n\n')
}else{ # When a seurat object is provided.
  # --> Load object.
  seurat.obj <- readRDS(file=data.file)
  # --> Subset if indicated.
  if(do.presubset){
    cat('Subsetting seurat object...\n')
    presubset.criteria <- read.csv(file=presubset.criteria.file, stringsAsFactors=FALSE)
    # Subset based on either of two strategies:
    tmp.check <- 'barcode' %in% colnames(presubset.criteria) & ncol(presubset.criteria) == 1
    # @ 1. Set of barcodes.
    if(tmp.check){
      tmp.cells <- presubset.criteria[, 'barcode']
      tmp.check <- all(tmp.cells %in% Cells(seurat.obj))
      if(!tmp.check) stop('A set of barcodes was provided to subset the seurat object but not all of these are defined in the object.\n')
      seurat.obj <- subset(x=seurat.obj, cells=tmp.cells)
    }else{
      presubset.criteria <- read.csv(file=presubset.criteria.file, stringsAsFactors=FALSE, row.names=1) # Re-load as usually.
    # @ 2. Specific criteria based on a predefined tag.
      if(!(colnames(presubset.criteria)[1] %in% colnames(seurat.obj@meta.data))) stop('Tag given as criteria to subset seurat object not appropriately defined in the seurat object metadata.\n')
      if(ncol(presubset.criteria)>1){tmp.warning <- paste0('More than a single column for subsetting in presubsetting criteria file. The program will attempt to subset using the first column.\n'); warning(tmp.warning)}
      seurat.obj <- subset.seurat.obj(seurat.obj=seurat.obj, tags.criteria=presubset.criteria, feats.criteria=NULL)
    }
    cat('Seurat object succesfully subsetted!\n')
  }
  # We'll clean the values, removing the ones that'll be updated somehow during this process.
  tags.to.keep <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
  seurat.obj@meta.data <- seurat.obj@meta.data[, tags.to.keep]
  cat('Seurat object has been loaded and cleaned!\n\n')
}
# Get feature info (name-ENSEMBL ID rels.) whenever it's needed.
# When we need to load the names of the genes.
if(feature.id=='ensembl'){
  feature.info.file <- unique(list.files(path=tenx.raw.data, pattern='features.tsv', full.names=TRUE))
  if(length(feature.info.file)!=1) stop('File with gene names -features.tsv- not properly defined.\n')
  is.zipped <- grepl(x=feature.info.file, pattern='\\.gz')
  if(is.zipped){
    # Then, read relatiosnships between feature names and ENSEMBL IDs.
    tmp.dir <- paste0(tenx.raw.data, '/tmp_with_id_', paste0(sample(x=LETTERS, size=5, replace=TRUE), collapse=''))
    tmp.cmmnd <- paste0('mkdir ',  tmp.dir, ' && cp ', feature.info.file, ' ', tmp.dir)
    system(command=tmp.cmmnd)
    tmp.file.name <- paste0(tmp.dir, '/features.tsv')
    tmp.cmmnd <- paste0('gunzip ', tmp.file.name, '.gz')
    system(tmp.cmmnd)
    feature.info <- read.delim(file=tmp.file.name, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
    tmp.cmmnd <- paste0('rm -r ', tmp.dir)
    system(tmp.cmmnd)
  }else{
    feature.info <- read.delim(file=feature.info.file, header=FALSE, col.names=c('ensembl', 'name', 'type'), stringsAsFactors=FALSE)
  }
}
# ---> Annotations table.
if(do.annotate){
  if(!file.exists(aggr.table.file)) stop('Aggregation table file does not exist.\n')
  aggr.table <- read.csv(file=aggr.table.file, stringsAsFactors=FALSE)
  cat('Annotations for seurat object cells have been loaded!\n\n')
  if(!is.null(donors.data.file)){
    donors.data <- read.csv(file=donors.data.file, stringsAsFactors=FALSE)
  }
}
# ---> Rules to create new tags.
new.tags.rules <- if(!is.null(new.tags.file)) read.csv(file=new.tags.file, stringsAsFactors=FALSE) else NULL
# ---> Subset criteria
if(do.subset){
  # Load criteria according to non-NULL values.
  tags.criteria <- if(!is.null(tags.criteria.file)) read.csv(file=tags.criteria.file, stringsAsFactors=FALSE) else NULL
  feats.criteria <- if(!is.null(feats.criteria.file)) read.csv(file=feats.criteria.file, stringsAsFactors=FALSE, row.names=1) else NULL
  # Make sure we've got any kind of criteria. Otherwise, finish program.
  if(is.null(feats.criteria)&is.null(tags.criteria)) stop('No criteria defined.\n')
  cat('Subset criteria loaded...\n\n')
}
# ---> Variables to regress.
# Load only if it's a valid path to a file.
if(file.exists(vars.to.regress)) vars.to.regress <- read.csv(file=vars.to.regress, stringsAsFactors=FALSE)[, 1]

if(do.annotate){
  cat('\n\n')
  ### ------------------------- Annotations -------------------------- ###
  cat('### ------------------------- Annotations -------------------------- ###\n')
  seurat.obj <- annotate.seurat.obj(seurat.obj=seurat.obj, aggr.table=aggr.table, lane.tag.id=lane.tag.id, chrom.batch.tag.id=chrom.batch.tag.id, seq.batch.tag.id=seq.batch.tag.id, hto.tag.id=hto.tag.id, overlay.tag.id=overlay.tag.id, mult.ann.per.lane=mult.ann.per.lane, exceded.batches=exceded.batches)
  cat('Seurat object has got annotations.\n')
  # ---> Annotations based on donors' data.
  if(!is.null(donors.data.file)){
    if(!(donors.merge.by %in% colnames(donors.data) & donors.merge.by %in% colnames(seurat.obj@meta.data))) stop('Donors metadata was provided but not an appropriate value to add further annotations (i.e., value was not part of seurat annotations.)')
    # ---> Find per-cell further annotations.
    tmp.df <- data.frame(col.1=seurat.obj@meta.data[, donors.merge.by], barcode=Cells(seurat.obj), stringsAsFactors=FALSE)
    colnames(tmp.df)[1] <- donors.merge.by
    tmp.merge <- merge(x=tmp.df, y=donors.data, by=donors.merge.by, all.x=TRUE)
    rownames(tmp.merge) <- tmp.merge$barcode
    tmp.merge$barcode <- NULL
    tmp.merge <- tmp.merge[Cells(seurat.obj), ]
    if(!all(rownames(tmp.merge) == Cells(seurat.obj))) stop('Something went wrong while adding patients\' annotations.')
    if(is.null(dim(tmp.merge[, colnames(tmp.merge)!=donors.merge.by]))){
      tmp.col.name <- colnames(tmp.merge)[colnames(tmp.merge)!=donors.merge.by]
      tmp.merge <- data.frame(tmp.merge[, colnames(tmp.merge)!=donors.merge.by])
      colnames(tmp.merge) <- tmp.col.name
    }else{
      tmp.merge <- tmp.merge[, colnames(tmp.merge)!=donors.merge.by]
    }
    tmp <- cbind(seurat.obj@meta.data, tmp.merge)
    # ---> Add all annotations.
    seurat.obj@meta.data <- tmp
    rm(tmp)
  }
  # ---> New tags.
  # Proceed only if required.
  if(!is.null(new.tags.rules)){
    seurat.obj <- add.new.tags(seurat.obj=seurat.obj, new.tags.rules=new.tags.rules, rem.tags=TRUE)
  }
}

### ---------------------- Preprocessing info ----------------------- ###

# ---> For cell counts.
# Retrieve all ID-related tags.
id.tags <- grep(x=colnames(seurat.obj@meta.data), pattern='\\.id\\.', value=TRUE)
# Confirm we've got valid data for each ID.
tmp.check <- unlist(lapply(X=id.tags, FUN=function(id.tag){
  !all(is.na(seurat.obj@meta.data[, id.tag]))
}))
id.tags <- id.tags[tmp.check]
# Proceed with this section only if any ID-related tags were identified.
if(length(id.tags)>0){
  # Create directory for general reports of the project.
  cell.counts.path <- paste0(reports.path, '/cell_count_summs')
  dir.create(cell.counts.path)
  # ---> Cell count per ID value at the pre-subsetting stage.
  cell.counts <- lapply(X=id.tags, FUN=function(id.tag){
    tmp.data <- as.data.frame(table(seurat.obj@meta.data[, id.tag]), stringsAsFactors=FALSE)
    colnames(tmp.data) <- c(id.tag, 'pre.subset.count')
    return(tmp.data)
  })
  names(cell.counts) <- id.tags
}

if(do.subset){
  cat('\n\n')
  ### ------------------------ Tags Filtering ------------------------ ###
  cat('### ------------------------ Tags Filtering ------------------------ ###\n')
  cat('Subsetting seurat object according to the criteria supplied...\n')
  # Subset based on either of two strategies:
  tmp.check <- 'barcode' %in% colnames(tags.criteria) & ncol(tags.criteria) == 1
  # @ 1. Set of barcodes.
  if(tmp.check){
    tmp.cells <- tags.criteria[, 'barcode']
    tmp.check <- all(tmp.cells %in% Cells(seurat.obj))
    if(!tmp.check) stop('A set of barcodes was provided to subset the seurat object but not all of these are defined in the object.\n')
    seurat.obj <- subset(x=seurat.obj, cells=tmp.cells)
  }else{
  # @ 2. Specific criteria based on a predefined tag.
    tags.criteria <- read.csv(file=tags.criteria.file, stringsAsFactors=FALSE, row.names=1) # Re-load as usually.
    # Check tags and features definition, as well as their possible values.
    if(!all(c(colnames(feats.criteria), colnames(tags.criteria)) %in% colnames(seurat.obj@meta.data))) stop('Some criteria term may not be properly defined in seurat object.')
    # Check rownames for both criteria terms.
    if(!all(rownames(feats.criteria)=='lower'|rownames(feats.criteria)=='upper')) stop('In features criteria, terms not properly defined as lower or upper thresholds.')
    if(!all(rownames(tags.criteria)=='discard'|rownames(tags.criteria)=='keep')) stop('In tags criteria, terms not properly defined as being criteria to discard or keep cells.')
    # Subset.
    seurat.obj <- subset.seurat.obj(seurat.obj=seurat.obj, tags.criteria=tags.criteria, feats.criteria=feats.criteria)
  }
  cat('Seurat object subsetted!\n')
  # ---> Cell count per ID value at the pre-filtering stage.
  # Proceed with this section only if any ID-related tags were identified.
  if(length(id.tags)>0){
    cell.counts <- lapply(X=id.tags, FUN=function(id.tag){
      tmp.data.1 <- cell.counts[[id.tag]]
      tmp.data.2 <- as.data.frame(table(seurat.obj@meta.data[, id.tag]), stringsAsFactors=FALSE)
      colnames(tmp.data.2) <- c(id.tag, 'pre.filter.count')
      tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=id.tag)
      return(tmp.data)
    })
    names(cell.counts) <- id.tags
  }
}


cat('\n\n')
### --------------------- Quality control (QC) --------------------- ###
cat('### --------------------- Quality control (QC) --------------------- ###\n')
# Create a directory to save QC results.
QC.path <- paste0(reports.path, '/QC')
create.dir(QC.path, 'Quality Control')

# General QC analysis.
seurat.obj <- seurat.qc.analysis(this.seurat.obj=seurat.obj)

# ---> Cell count per ID value at the post-filtering stage.
# Proceed with this section only if any ID-related tags were identified.
if(length(id.tags)>0){
  cell.counts <- lapply(X=id.tags, FUN=function(id.tag){
    tmp.data.1 <- cell.counts[[id.tag]]
    tmp.data.2 <- as.data.frame(table(seurat.obj@meta.data[, id.tag]), stringsAsFactors=FALSE)
    colnames(tmp.data.2) <- c(id.tag, 'post.filter.count')
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=id.tag)
    return(tmp.data)
  })
  names(cell.counts) <- id.tags
}

# ---> Comprehensive cell count summary per ID-related tag per pipeline filtering-related stage.
# Proceed with this section only if any ID-related tags were identified.
if(length(id.tags)>0){
  lapply(X=id.tags, FUN=function(id.tag){
    # Retrieve summary.
    tmp.data <- cell.counts[[id.tag]]
    # Output raw counts.
    tag.lab <- str_replace(string=id.tag, pattern='\\.tag$', replacement='')
    tag.lab <- str_replace_all(string=tag.lab, pattern='\\.', replacement='-')
    tmp.file.name <- paste0(cell.counts.path, '/CellCountSummaryForTag-', tag.lab, '.csv')
    write.csv(file=tmp.file.name, x=tmp.data, row.names=FALSE, quote=FALSE)
    # Output visual summary.
    colnames(tmp.data)[1] <- 'tmp.tag'
    tmp.data <- tidyr::gather(data=tmp.data, -`tmp.tag`, key='stage', value='cell.count')
    tmp.data$stage <- str_replace(string=tmp.data$stage, pattern='\\.count$', replacement='')
    filt.stages <- filt.stages[names(filt.stages) %chin% tmp.data$stage]
    tmp.data$stage <- factor(x=filt.stages[tmp.data$stage], levels=filt.stages)
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=stage, y=cell.count, col=tmp.tag, group=tmp.tag)) + geom_point() + geom_line() + labs(x='Filtering stage', y='Cell count', col='ID value')
    tmp.file.name <- paste0(cell.counts.path, '/CellCountSummaryForTag-', tag.lab, '.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot + theme_bw())
    dev.off()
    # Mock return.
    return(NULL)
  })
}

### -------------------- Feature selection and ------------------- ###
### ----------------------- Normalization ------------------------ ###
cat('### -------------------- Feature selection and ------------------- ###\n')
cat('### ----------------------- Normalization ------------------------ ###\n')
# Create a directory to save feature selection reports.
feat.select.path <- paste0(reports.path, '/feature_selection')
create.dir(feat.select.path, 'Feature selection')

# DEPENDING ON NORMALIZATION METHOD SELECTED.
switch(EXPR=gen.norm.method,
### ------------------------ SCTransform ------------------------- ###
'SCTransform'={
  cat('\n### ------------------------ SCTransform ------------------------- ###\n')
  seurat.obj <- SCTransform(seurat.obj, vars.to.regress="percent.mt", verbose=TRUE)
  # Set SCT assat as the working assay.
  wking.assay <- 'SCT'
  # Set final variable features as the default ones identified during the process.
  final.variable.feats <- VariableFeatures(seurat.obj)
},
### ------------------------ LogTransform ------------------------ ###
'LogTransform'={
  cat('\n### ------------------------ LogTransform ------------------------ ###\n')
  tmp.output <- do.normalization() # Here we get a list with 1) the processed seurat object and 2) the vector of final variable features.
  # Set SCT assat as the working assay.
  wking.assay <- 'RNA'
  # Set final variable features as the default ones identified during the process and the seurat object as processed along the process.
  final.variable.feats <- tmp.output[['var.feats']]
  seurat.obj <- tmp.output[['seurat.obj']]
  rm(tmp.output)
},
stop('No appropriate normalization method input.\n')
)
cat('\n\nNormalization completed!\n\n')

### ------------------ Dimensionality reduction ------------------ ###
cat('### ------------------ Dimensionality reduction ------------------ ###\n')
# Create a directory to save dimensionality reduction reports.
dim.reduction.path <- paste0(reports.path, '/dim_reduction')
create.dir(dim.reduction.path, 'Dimensionality reduction')

# Run initial dimensionality reduction (PCA and, if requested, Harmony batch-effect removal)
seurat.obj <- init.dim.reduction()

# Downstream analyses will be performed separately for each value of number of PCs passed to the general program. Therefore, there will also be separate folders for downstream results.
### ----------------------- Save processed object ------------------------ ###
cat('<<< ----------------------- Save processed object ------------------------ >>>\n')
cat('Seurat objects already processed should and will be saved.\n')
seurat.objects.path <- paste0(reports.path, '/seurat_objects')
create.dir(seurat.objects.path, 'Objects')

for(no.PCs in nos.PCs){
cat('### -------------------------- PC value:', no.PCs, '--------------------------- ###\n')
tmp.dim.reduction.path <- paste0(dim.reduction.path, '/dim_reduction_for_', no.PCs, 'PCs')
create.dir(tmp.dim.reduction.path, paste0('Dimensionality reduction for ', no.PCs, ' PCs'))

cat('\n\n')
### -------------------------- Cells clustering -------------------------- ###
### -------------------------------- and --------------------------------- ###
### ------------------ Non-linear dimensional reduction ------------------ ###
cat('### -------------------------- Cells clustering -------------------------- ###\n')
cat('### -------------------------------- and --------------------------------- ###\n')
cat('### ------------------ Non-linear dimensional reduction ------------------ ###\n')

seurat.obj <- dim.reduction.clustering()


cat('\n\n')
### ------------------------ Clusters definition ------------------------- ###
### -------------------------------- and --------------------------------- ###
### ------------------------ Cluster-specific QC ------------------------- ###
cat('### ------------------------ Clusters definition ------------------------- ###\n')
cat('### -------------------------------- and --------------------------------- ###\n')
cat("### ------------------------ Cluster-specific QC ------------------------- ###\n")

clusters.def.and.qc(this.seurat.obj=seurat.obj)

### ----------------------- Save processed object ------------------------ ###
cat('<<< ----------------------- Save processed object ------------------------ >>>\n')
# Saving object.
tmp.file.name <- get.file.name(gen.path=seurat.objects.path, file.desc=paste0('SeuratObjectForPrj', prj.name, '_WithArgs'), extension='RDS', no.PCs=no.PCs)
saveRDS(object=seurat.obj, file=tmp.file.name)
cat('Object saved as an RDS file!\n\n')
} # Here finishes the separate analyses for each value of number of PCs.

### ------------------------ Defining populations ------------------------ ###
cat('### ------------------------ Defining populations ------------------------ ###\n')
cat('These analyses will be saved to cell assignment directory nested to the project directory.\n')
# New directory for DGEA.
sc.assign.path <- paste0(reports.path, '/population_definition')
create.dir(dir.path=sc.assign.path, path.desc='Population definition')

cat('The analysis for population definition for clusters available in this version are next:\n\tDifferential gene expression analyses.\n\tCiro\' original significance plots.\nThis will be applied over the clusters defined by the different resolutions values.\n')

cat('\n')
# 1. DGEA and significance plots ---------------------------------------------
cat('# 1. DGEA and significance plots ---------------------------------------------\n')
# New directory to save the significance plots.
significance.path <- paste0(sc.assign.path, '/significance_plots')
create.dir(dir.path=significance.path, path.desc='Significance plots')
# New directory for DGEA.
DGEA.path <- paste0(sc.assign.path, '/DGEA')
create.dir(dir.path=DGEA.path, path.desc='Differential Gene Expression Analysis')

# This process will be applied for the different sets of clusters according the resolution values applied.
cat('The script will attempt to run this part of the job in a parallel way...\n\n')
# If necessary, downsample cells from the seurat objects.
if(length(Cells(seurat.obj)) > dea.max.cells){
  tmp.sample <- sample(x=Cells(seurat.obj), size=dea.max.cells, replace=FALSE)
  seurat.obj <- subset(x=seurat.obj, cells=tmp.sample)
}
# Proceed with DGEA
for(resolution in resolutions) do.pop.definition(resolution=resolution)
cat('Population definition analyses processed!\n\n')

cat('\n')
# 2. Populations tree --------------------------------------------------------
cat('# 2. Populations tree --------------------------------------------------------\n')
clusts.tree.path <- paste0(sc.assign.path, '/aid_to_pick_res')
create.dir(dir.path=clusts.tree.path, path.desc='Populations tree')

out.pops.tree(seurat.obj=seurat.obj, reports.path=clusts.tree.path, gen.norm.method=gen.norm.method)

cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
