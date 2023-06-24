Seurat-based scRNA-seq Analysis for large-scale datasets (v2.3)
===========

Developed for the analysis of the [DICE database](https://dice-database.org) datasets. Check out our manuscripts:
- [*Single-cell eQTL analysis of activated T cell subsets reveals activation and cell type–dependent effects of disease-risk variants*](https://www.science.org/doi/10.1126/sciimmunol.abm2508#)
- DICE Tissue (TBD)

---
# System requirements 

Please consult the <code>Dependencies</code> file for more details on R dependencies. Originally developed on a Linux platform.

---
# Run example (TO DO)

Please see the <code>RunExample</code> file.

---
# Arguments per step

## General arguments (input definition)
The main input for the program is the UMI counts matrix, an output from 10x's cellranger (see, e.g., [here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count)). The rows in the matrix represent genes and the columns represent cell barcodes. The IDs of the genes might be taken as their gene symbols (referred to here as "names") or their Ensembl IDs. See the options below.
- <code>ReportsPath</code> Char vector, absolute path to directory to save analyses reports.
- <code>PrjName</code> Char vector, name to be given to the project.
- <code>DataFile</code> Char vector, absolute path to file storing a feature by counts matrix (format provided by 10x's cellranger)
- <code>InputType</code> Char vector, indicates what kind of input should be considered, either 'counts' for a count matrix or 'seurat' for a seurat object.
- <code>FeatureID</code> Char vector, defines the feature ID type to take from the counts matrix. Default is 'name' (features.tsv column 2). Possible values are 'name' and 'ensembl', the last one corresponding to the Ensembl ID. This option applies only if a count matrix (and not a seurat object) is provided.
- <code>Raw10x</code> Char vector, absolute path to 10x raw data path (as a directory instead of an hf5 matrix file). Necessary when \"ensembl\" is picked as the feature ID option.
- <code>IsFivePrime</code> Logical, indicates whether the data comes from a five-prime experiment, in which case, antigen receptor (either T-cell receptor or B-cell receptor) genes will be excluded from downstream analyses.
- <code>MarkersFile</code> Char vector, absolute path to RDS file listing gene markers. These markers will be assessed. Format is described elsewhere and an example provided along.
- <code>RAM</code> Int, indicates available RAM for the job session.

## For annotations
If the data from multiple libraries were aggregated for joined analysis, certain knowledge about the barcodes' origin is known apriori based on library of origin. These knowledge can be provided as metadata annotations during the process. To do so, the program depends on an aggregation table, which must have the same format as the one required by [cellranger aggr](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/aggregate) plus a few columns as shown in the example provided here.
- <code>DoAnnotate</code> Logical, indicates whether to perform metadata annotation based on aggr table (format detailed elsewhere).
- <code>AggrTable</code> Char vector, absolute path to fully-annotated aggregation table (format detailed elsewhere and an example provided along).
- <code>DonorsMetaData</code> Char vector, absolute path to donors metadata file. Set as NULL if there's none (format detailed elsewhere and an example provided along).
- <code>ExcededBatches</code> Logical, indicates whether it's ok having more batches defined in the aggr table provided than the actual number of batches defined for the seurat object. If FALSE, the pipeline will get stalled when there are less batches defined (in the seurat object) than the ones expected according to the aggr table.
- <code>MergeBy</code> Char vector, column name in both seurat object meta data and donors meta data in order to merge them and add further annotations. Value considered only if there are annotations and donors metadata provided
- <code>LaneID</code> Char vector, represents teh names to be given to the columns listing the values defined in the column 'lane.tag.id' of the aggregation table. Must be a single string where different column names are separated by semicolons.
- <code>MultAnnPerLane</code> Logical, indicates wether the lane tag column in the aggregation table lists multiple values. If TRUE, they should be separated by semicolon within each aggregation table lane index and there should be the same amount of values per library (i.e., per row in the aggregation table).
- There are more detailed options detailed under this section within the script.

## Subsetting options
- <code>DoPreSubset</code> Logical, indicates whether presubsetting should be applied or not. Only relevant when a seurat object is input. Otherwise, the value is completely ignored.
- <code>PreSubsetCriteria</code> Char vector, absolute path to file indicating the criteria used to apply presubsetting. Only relevant when a seurat object is input. Two different subset strategies are valid depending on the columns defined in the file. If the file has a single column (named \"barcode\"), then the barcodes listed there will be taken as the subset to be kept for downstream analyses. Otherwise, the file has to have the same format as for Tag and Features criteria (see options 'DoSubset', 'TagsCriteria' and 'FeatsCriteria' below) and the same subsetting strategy is applied.
- <code>DoSubset</code> Logical, flag indicates whether the seurat object must be filtered by excluding a subset of cells or a subset of genes/features. Rules to remove or keep cells are specified through options 'TagsCriteria' and 'FeatsCriteria' (see below).
- <code>TagsCriteria</code> Char vector, absolute path to file indicating the criteria used to apply presubsetting according to barcodes. Two different subset strategies are valid depending on the columns defined in the file. If the file has a single column (named \"barcode\"), then the barcodes listed there will be taken as the subset to be kept for downstream analyses. Otherwise, the file has to have the same format as for Tag and Features criteria and the same subsetting strategy is applied.
- <code>FeatsCriteria</code> Char vector, absolute path to table storing the criteria to filter cells based on for tags required (format: criteria by tag, column name should match that of tag defined in seurat object; criteria should be separated in two different rows, lower and upper threshold respectively).

## For quality control (QC)
For QC, three barcode-specific measures are used as features to define low-quality barcodes which are likely to represent one of the following: ambient RNA, dying (perforated) cells or multiplets. These features are: (i) number of UMIs per cell (sometimes also referred to as barcode-specific sequencing depth); (ii) number of genes/features per cell; and (iii) percentage of UMIs that mapped to mitochondrial genes. For further details see [ref](https://www.embopress.org/doi/abs/10.15252/msb.20188746). During the analysis, user may define thresholds for these features to identify and remove low-quality cells. The following options control this.
- <code>minCounts</code> Integer, indicates the lower threshold for number of UMIs per cell.
- <code>maxCounts</code> Integer, indicates the upper threshold for number of UMIs per cell.
- <code>minFeatures</code> Integer, indicates the lower threshold for number of genes per cell.
- <code>maxFeatures</code> Integer, indicates the upper threshold for number of UMIs per cell.
- <code>maxMP</code> Integer, indicates the upper threshold for percentage of mitochondrial UMIs per cell.
- <code>FilterOut</code> Logical, indicates whether the barcodes identified to have low quality (based on the criteria defined with the options above) must be filtered out for downstream analysis.\nFor a pilot run, it is highly recommend to keep even the low-quality cells.

## For normalization
Before applying the dimensionality reduction and clustering steps, UMI counts must be normalized. Several normalization strategies have been proposed, but two in particular have become popular thanks to being available through the Seurat toolkit, namely LogNormalize and SCTransform [(ref)](https://link.springer.com/article/10.1186/s13059-019-1874-1). Both options are available through this pipeline and the following options provide control over this.
- <code>FeatsForDSA</code> Int, indicates one of two options depending on whether the value is over 100. If larger than 100, it indicates the amount of variable features to take into consideration for downstream analysis. If smaller or equal to 100, it indicates the percentage of standardized variance to be considered such that the most variable features that account for that percentage of variance will be taken into consideration.
- <code>GenNormMethod</code> Char vector, indicates the method to use for general data normalization. Possible values (other than these, program will stop): LogTransform (Log-normalization for cell-wise normalization proceeded by other steps such as gene normalization) and SCTransform [(ref)](https://link.springer.com/article/10.1186/s13059-019-1874-1).
- <code>VarsToRegress</code> Character, either 1) an absolute path to a csv file (single column, with header 'features') listing a ser of variables to regress during gene normalization (scaling applied with the function ScaleData from Seurat, passed to the argument named the same way) or a single string that lists the values to regress out during normalization (separated by semicolons).

## For Dimensionality reduction and clustering.
Two kinds of dimensionality reduction are applied during the process, linear dimensionality reduction (essentially, PCA) and non-linear (graph-based) dimensionality reduction. For the latter, a speciifc set of principal components (PCs) might be provided. Previous to clustering, other strategies to correct for potential batch effects have been suggested, [Harmony](https://github.com/immunogenomics/harmony/tree/master) as one of the most popular ones. You can apply harmony for batch effect removal in here. For the non-linear dimensionality reduction process, two methods are available: (i) UMAP embedding and (ii) tSNE embedding, where both options are avaible here (either of them or both). Finally, during clustering, the resolution parameter controls the number of clusters inveiled by the method, such that the higher the resolution is, the larger the number of clusters that will be recovered. The options below control for all of these stages.
- <code>PCs</code> Int, indicates the number of PCs to take into consideration for non-linear dimensionality reduction and clustering.
- <code>ForHarmony</code> Char vector, if different than NULL, it indicates that Harmony batch-effect removal should be applied. Then, this may be evaluated as a character vector of the metadata column listing the batches to be regressed out through Harmony.
- <code>Resolution</code> Char vector. A single string to be interpreted as a numeric vector where each value must be a nonnegative value between 0.1 and 1. For example, the string 'c(0.2, 0.4)' will be interpreted as a numeric vector with entries 0.2 and 0.4. These indicate the resolutions to consider during clustering. There is not a limit in the number of resolutions to try, but be aware that the more you provide, the longer it will take for the script to run.
- <code>DimReduction</code> String to be evaluated as a character vector indicating the dimensionality reduction methods that should be applied over the data other than PCA (which will run tried all the same). Possible values are 'umap' and 'tsne', where either or both can be provided. Default is both.
- <code>MaxCellsNo</code> Int, indicates the maximum number of cells to depict in plots.

---
# Troubleshooting
For further details on the script options or troubleshooting, please feel free to open an issue or to email Vicente directly (vfajrdo@lji.org)
