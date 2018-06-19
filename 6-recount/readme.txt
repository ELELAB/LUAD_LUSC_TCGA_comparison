README:

* The script files unifiedLUAD.R and unifiedLUSC.R should be run to download Recount2 project
data (GTEx+ TCGA). Meta data and gene expression data are saved into .RData files after running the script

* Gene names in the data produced have been strapped off their version number ( content after "." in ENSG ids), please
edit script if you want gene IDs to remain untouched.

* After preparation, data is discordance-free and contain only samples with 60% tumor purity or more

(1) run first the scripts unifiedLUAD.R and unifiedLUSC.R in LUAD AND LUSC directory, respectively. We also provided the option to apply the Rail transformation suggested by recount developers.
(2) go to LUAD or LUSC sub-folders to continue with the DE analyses 
(3) then go the the VennDiagram_unified_LUAD-LUSC  subfolder to define the genes unique for LUAD and LUSC, respectively
(4) finally use the DEGcomparison_unified-tcga subfolder to compare these results with the one from the DEA on TCGA datasets

NOTES
(1) in the LUSC folder, the user will need to have the file with information on the LUSC discordant samples.

