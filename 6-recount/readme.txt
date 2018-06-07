READ ME:

* The script files unifiedLUAD.R and unifiefLUSC.R should be run to download Recount2 project
data (GTEx+ TCGA). Meta data and gene expression data are saved into .RData files after running the script

* Gene names in the data produced have been strapped off their version number ( content after "." in ENSG ids), please
edit script if you want gene IDs to remain untouched.
---> comment line 90 in both files

* After preparation, data is discordance-free and contain only samples with 60% tumor purity or more