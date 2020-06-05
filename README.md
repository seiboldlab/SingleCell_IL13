# SingleCell_IL13
Herein lies code from the paper, *Single cell and population transcriptomics reveal pan-epithelial remodeling in type 2-high asthma*.

Scripts are in the 'scripts' directory and data referenced in the scripts are in the 'data' directory. However, when starting from the beginning, the original count matrices referenced will first need to be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145013). For the single cell datasets, subsequent analyses can be run by first loading saved R files containing the relevant datasets (Seurat objects), which can also be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145013).

## Scripts
**01_SeuratAnalysis_48hr.R** - Data QC, count normalization, and cluster inference for the acute (48 hour) IL-13 single cell sequencing experiment.

**02_SeuratAnalysis_11d.R** - Data QC, count normalization, dataset alignment, and cluster inference for the chronic (11 day) IL-13 single cell sequencing experiment.

**03_MonocleAnalysis_48hr.R** - Pseudotime trajectory analysis for defense and mucus secretory populations from the acute (48 hour) IL-13 single cell sequencing experiment.
