# SingleCell_IL13
Code for "Single cell and population transcriptomics reveal pan-epithelial remodeling in type 2-high asthma"

Scripts are in the scripts directory and data referenced in the scripts are in the data directory. However note that original count matrices will first need to be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145013).

## Scripts
01_SeuratAnalysis_48hr.R - data QC, count normalization, and cluster inference for the acute (48 hour) IL-13 single cell sequencing experiment.

02_SeuratAnalysis_11d.R - data QC, count normalization, dataset alignment, and cluster inference for the chronic (11 day) IL-13 single cell sequencing experiment.
