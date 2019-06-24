In this folder we use a deconvolution method to estimate the abundance of tissue-infiltrating immune cell populations, than we perform the DEA taking into account these estimates. We also perform the pathway and GO enrichment analysis.


1) 1-DEA
We provide the output after running the 1-DEA_deconvolution.R script. Each subfolder includes the following files:

- edgeR_LUAD_’name CellType’.csv
- edgeR_LUSC_’name CellType’.csv
- down_edgeR_LUAD_’name CellType’.txt
- down_edgeR_LUSC_’name CellType’.txt
- up_edgeR_LUAD_’name CellType’.txt
- up_edgeR_LUSC_’name CellType’.txt

The 2-difference_LUAD_LUSCgenes.R script selects the non-overlap DEA genes between LUAD and LUSC, separating the up- and down-regulated genes. Each subfolder includes also the following files :

- up_LUAD_only.txt
- down_LUAD_only.txt
- up_LUSC_only.txt
- down_LUSC_only.txt

2) 2-Reactome
In this folder we perform the pathway enrichment analysis for each immune population, separating the up- and down-regulated genes and the cancer type (LUAD and LUSC).
Each subfolder includes the following files which contain the significant pathways (FDR < 0.05):

- pathway_table_down_LUAD.txt
- pathway_table_up_LUAD.txt
- pathway_table_down_LUSC.txt
- pathway_table_up_LUSC.txt

If one or more of these files are not present, it means no significant pathways have been detected for the specific immune population.

3) 3-GO
In this folder we perform the GO enrichment analysis. It includes the GO results stored in the RData and cycle plots (png).
The limma_LUAD_all_tss_tumorPurity.csv and limma_LUSC_all_tss_tumorPurity.csv files are used to get the logFC used in the circle plots


Before running the 1-DEA_deconvolution.R and pathwayEnrichhment.R scripts, the user needs to create the B-lineage, Lymphocytes, Monocytic, Myeloid_cells, Neutrophils and Tcells subfolders.

Therefore the user needs to have in this folder the files:
(a) MCPcounter_functions.R
(b) TCGAbiolinks_functions.R
(c) LUAD_PreprocessedData_all_tumorPurity.rda
(d) LUSC_PreprocessedData_all_tumorPurity.rda
(e) LUAD_Illumina_HiSeq_all_tumorPurity.rda
(f) LUSC_Illumina_HiSeq_all_tumorPurity.rda



