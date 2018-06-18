In this folder the script and steps to run the coexpression analyses are included.

We also provided the output after runnning the scripts as a reference
Before running the scripts the user need to have the proper structure of sub-directories in place - comments are provided about this within the scripts.

The scripts to be run are, in order:
(1) download_script_tp_LUAD.R
(2) download_script_tp_LUSC.R
(2) Processing_script_LUAD_TP_CEMiTool.R
(3) Processing_script_LUSC_TP_CEMiTool.R 


You need to have in this folder also the files:
(a)TCGAbiolinks_functions.R
(b) i2d_geneName.txt, which contains the I2D human protein-protein interaction database with UNIPROT ID converted to gene names 

The following files contain the pre-processed summarized experiment datasets for each case of study:
LUAD_Illumina_HiSeq_TP.rda 
LUSC_Illumina_HiSeq_TP.rda


