In this folder we use a deconvolution method to estimate the abundance of tissue-infiltrating immune cell populations, than we perform the DEA taking into account these estimates.

Before running the DEA_deconvolution.R script, the user needs to create the B-lineage, Lymphocytes, Monocytic, Myeloid_cells, Neutrophils and Tcells subfolders.

We provided the output after running the script. Each subfolder includes the following files:

- edgeR_LUAD_’name CellType’.csv
- edgeR_LUSC_’name CellType’.csv
- down_edgeR_LUAD_’name CellType’.txt
- down_edgeR_LUSC_’name CellType’.txt
- up_edgeR_LUAD_’name CellType’.txt
- up_edgeR_LUSC_’name CellType’.txt

The user needs to have in this folder the files:
(a) MCPcounter_functions.R
(b) TCGAbiolinks_functions.R
(c) LUAD_PreprocessedData_all_tumorPurity.rda
(d) LUSC_PreprocessedData_all_tumorPurity.rda  



