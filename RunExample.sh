# -------   Examples of how to run the script using the input files listed in this repo     ------- #

# ---> Run standard analysis on the output from the cellranger aggr. QC filtering is applied (because of --FilterOut TRUE) based on the cutoffs provided for the QC features.
Rscript /path/to/repo/Seurat-based_scRNA-seq_Analysis_v2.3/R/master/general_seurat_analysis.2.3.R --ReportsPath /path/to/directory/to/save/output --PrjName AllCellTypes \
    --DataFile /path/to/output/from/cellranger/aggr/outs/count/filtered_feature_bc_matrix --InputType matrix --RAM 350 \
    --DoAnnotate TRUE --AggrTable /path/to/aggr_table_annotated_example.csv --LaneID 'species.tag;cell.type.tag;activation.status.tag' --MultAnnPerLane TRUE \
    --DonorsMetaData /path/to/SamplesMetadata.csv --MergeBy donor.tag \
    --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag \ # Usually there's no need to modify the values in this row.
    --DoSubset TRUE --TagsCriteria /path/to/TagsSubsetCriteria.csv \
    --FilterOut TRUE --minCounts 1500 --maxCounts 20000 --minFeatures 800 --maxFeatures 4400 --maxMP 1500 \
    --FeatsForDSA 25 --PCs 30 --Resolution 'c(0.2, 0.3)' --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")' --DimReduction 'c("umap")' \
    --MarkersFile /path/to/T-cell_markers.1.0.RData

# ---> Run standard analysis when QC filtering is not applied (because of --FilterOut FALSE) based on the cutoffs provided for the QC features.
Rscript /path/to/repo/Seurat-based_scRNA-seq_Analysis_v2.3/R/master/general_seurat_analysis.2.3.R --ReportsPath /path/to/directory/to/save/output --PrjName AllCellTypes \
    --DataFile /path/to/output/from/cellranger/aggr/outs/count/filtered_feature_bc_matrix --InputType matrix --RAM 350 \
    --DoAnnotate TRUE --AggrTable /path/to/aggr_table_annotated_example.csv --LaneID 'species.tag;cell.type.tag;activation.status.tag' --MultAnnPerLane TRUE \
    --DonorsMetaData /path/to/SamplesMetadata.csv --MergeBy donor.tag \
    --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag \ # Usually there's no need to modify the values in this row.
    --DoSubset TRUE --TagsCriteria /path/to/TagsSubsetCriteria.csv \
    --FilterOut FALSE --minCounts 1500 --maxCounts 20000 --minFeatures 800 --maxFeatures 4400 --maxMP 1500 \
    --FeatsForDSA 25 --PCs 30 --Resolution 'c(0.2, 0.3)' --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")' --DimReduction 'c("umap")' \
    --MarkersFile /path/to/T-cell_markers.1.0.RData

# ---> Run analysis when harmony is applied for batch effect correction based on the sequencing batch info.
Rscript /path/to/repo/Seurat-based_scRNA-seq_Analysis_v2.3/R/master/general_seurat_analysis.2.3.R --ReportsPath /path/to/directory/to/save/output --PrjName AllCellTypes \
    --DataFile /path/to/output/from/cellranger/aggr/outs/count/filtered_feature_bc_matrix --InputType matrix --RAM 350 \
    --DoAnnotate TRUE --AggrTable /path/to/aggr_table_annotated_example.csv --LaneID 'species.tag;cell.type.tag;activation.status.tag' --MultAnnPerLane TRUE \
    --DonorsMetaData /path/to/SamplesMetadata.csv --MergeBy donor.tag \
    --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag \ # Usually there's no need to modify the values in this row.
    --DoSubset TRUE --TagsCriteria /path/to/TagsSubsetCriteria.csv \
    --FilterOut TRUE --minCounts 1500 --maxCounts 20000 --minFeatures 800 --maxFeatures 4400 --maxMP 1500 \
    --FeatsForDSA 25 --PCs 30 --Resolution 'c(0.2, 0.3)' --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")' --DimReduction 'c("umap")' \
    --ForHarmony 'c("seq.batch.tag")'
    --MarkersFile /path/to/T-cell_markers.1.0.RData

# ---> Run standard analysis but gettng the Ensembl IDs instead of gene names as gene IDs (--FeatureID ensembl).
# NOTE: Here you have to modify the markers file accordingly (--MarkersFile).
Rscript /path/to/repo/Seurat-based_scRNA-seq_Analysis_v2.3/R/master/general_seurat_analysis.2.3.R --ReportsPath /path/to/directory/to/save/output --PrjName AllCellTypes \
    --DataFile /path/to/output/from/cellranger/aggr/outs/count/filtered_feature_bc_matrix --InputType matrix --FeatureID ensembl --RAM 350 \
    --DoAnnotate TRUE --AggrTable /path/to/aggr_table_annotated_example.csv --LaneID 'species.tag;cell.type.tag;activation.status.tag' --MultAnnPerLane TRUE \
    --DonorsMetaData /path/to/SamplesMetadata.csv --MergeBy donor.tag \
    --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag \ # Usually there's no need to modify the values in this row.
    --DoSubset TRUE --TagsCriteria /path/to/TagsSubsetCriteria.csv \
    --FilterOut TRUE --minCounts 1500 --maxCounts 20000 --minFeatures 800 --maxFeatures 4400 --maxMP 1500 \
    --FeatsForDSA 25 --PCs 30 --Resolution 'c(0.2, 0.3)' --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")' --DimReduction 'c("umap")' \
    --MarkersFile /path/to/T-cell_markers_ensembl.1.0.RData

# ----> Run analysis for a specific cluster identified in a previous analysis by taking the seurat object output from that analysis as input in the new one (--DataFile). 
Rscript /path/to/repo/Seurat-based_scRNA-seq_Analysis_v2.3/R/master/general_seurat_analysis.2.3.R --ReportsPath /path/to/directory/to/save/output --PrjName AllCellTypes \
    --DataFile /path/to/output/from/previous/analysis/SeuratObj.RDS --InputType seurat --RAM 350 \
    --DoAnnotate TRUE --AggrTable /path/to/aggr_table_annotated_example.csv --LaneID 'species.tag;cell.type.tag;activation.status.tag' --MultAnnPerLane TRUE \
    --DonorsMetaData /path/to/SamplesMetadata.csv --MergeBy donor.tag \
    --ChromBatchID chrom.batch.tag --SeqBatchID seq.batch.tag --HTOID hashtag.tag --OverlayID donor.tag \ # Usually there's no need to modify the values in this row.
    --DoSubset TRUE --TagsCriteria /path/to/TagsSubsetCriteria.csv \
    --FilterOut TRUE --minCounts 1500 --maxCounts 20000 --minFeatures 800 --maxFeatures 4400 --maxMP 1500 \
    --FeatsForDSA 25 --PCs 30 --Resolution 'c(0.2, 0.3)' --MeanCutoff 0.01 --PropCutoff 0.001 --VarsToRegress 'c("nCount_RNA", "percent.mt")' --DimReduction 'c("umap")' \
    --MarkersFile /path/to/T-cell_markers.1.0.RData