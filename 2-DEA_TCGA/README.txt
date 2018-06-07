In this folder, the scripts to carry DEA with the three pipelines used in the study are reported.
Moreover, we provide the outputs from the analyses.
The function used by TCGAbiolinks at the time of these analyses has been currently replaced so if the
results for the EdgeR_TCGAbiolink pipeline need to be reproduced we provide also the original function in the folder. It is called:

TCGAanalyze_DEA_old.R

To carry out the different DEA pipelines on the TCGA data, it is possible to use the following scripts:
(1)edgeR_TCGAb.R
(2)edgeR.R
(3)limma.R

For analyses of the data and production of figures such as the ones used in Figure 1 and 2 of the manuscript:
(1)  UpSetR_script_dataset.R 
(2)  UpSetR_script_methods.R 
(3)  comparison_FC_FDR.R 
(4)  scatterplot.R

NB In the folder LUAD_LUSC_DEA_revised the alternative scripts for the LIMMA-voom step described in the Methods of the article is also provided.
