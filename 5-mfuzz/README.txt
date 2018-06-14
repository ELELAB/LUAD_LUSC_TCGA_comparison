In this folder, the scripts to carry out mfuzz clustering are reported
Please first unzip the files LUAD/data/LUAD_all_consensus.csv.gz and LUSC/data/LUSC_all_consensus.csv.gz

To carry out the different DEA pipelines on the TCGA data, it is possible to use the following scripts:
(1)edgeR_TCGAb.R
(2)edgeR.R
(3)limma.R

For analyses of the data and production of figures such as the ones used in Figure 1 and 2 of the manuscript:
(1)  UpSetR_script_dataset.R 
(2)  UpSetR_script_methods.R 
(3)  VennDiagram_dataset.R
(4)  VennDiagram_three_methods.R
(5)  comparison_FC_FDR.R 
(6)  scatterplot.R
(7) final_DEgene_sets.R

NB In the folder LUAD_LUSC_DEA_revised the alternative scripts for the LIMMA-voom step described in the Methods of the article is also provided.
