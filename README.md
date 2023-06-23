Seurat-based scRNA-seq Analysis for large-scale datasets (v2.3)
===========

Developed for the analysis of the [DICE database](https://dice-database.org) datasets. Check out our manuscripts:
- [*Single-cell eQTL analysis of activated T cell subsets reveals activation and cell type–dependent effects of disease-risk variants*](https://www.nature.com/articles/s41592-019-0619-0)
- DICE Tissue (TBD)

# System requirements 

Please consult the <code>Dependencies</code> file for more details on R dependencies. Originally developed on a Linux platform.

# Arguments

### General arguments
- <code>ReportsPath</code> Char vector, absolute path to directory to save analyses reports.
- <code>PrjName</code> Char vector, name to be given to the project.
- <code>DataFile</code> Char vector, absolute path to file storing a feature by counts matrix (format provided by 10x's cellranger)
- <code>InputType</code> Char vector, indicates what kind of input should be considered, either 'counts' for a count matrix or 'seurat' for a seurat object.
- <code>FeatureID</code> Char vector, defines the feature ID type to take from the counts matrix. Default is 'name' (features.tsv column 2). Possible values are 'name' and 'ensembl', the last one corresponding to the Ensembl ID. This option applies only if a count matrix (and not a seurat object) is provided.
- <code>Raw10x</code> Char vector, absolute path to 10x raw data path (as a directory instead of an hf5 matrix file). Necessary when \"ensembl\" is picked as the feature ID option.
- <code>IsFivePrime</code> Logical, indicates whether the data comes from a five-prime experiment, in which case, antigen receptor (either T-cell receptor or B-cell receptor) genes will be excluded from downstream analyses.
- <code>MarkersFile</code> Char vector, absolute path to RDS file listing gene markers. These markers will be assessed. Format is described elsewhere and an example provided along.
- <code>RAM</code> Int, indicates available RAM for the job session.

### For annotations
- <code>DoAnnotate</code> Logical, indicates whether to perform metadata annotation based on aggr table (format detailed elsewhere).
- <code>AggrTable</code> Char vector, absolute path to fully-annotated aggregation table (format detailed elsewhere and an example provided along).

