In this folder the script and steps to assembly the TCGA datasets used in the study 
are illustrated.
We also provided the output after runnning the scripts as a reference
Before running the script 1-download.R, the user needs to create the LUAD and LUSC subfolders, each of them containing additional subfolders, i.e. all, paired and unpaired.
The scripts to be run are, in order:
(1) 1-download_script.R
(2) 2-preprocessing_LUAD.R
(2) 2-preprocessing_LUSC.R


You need to have in this folder also the files:
(a)TCGAbiolinks_functions.R
(b) discordant_LUSC.txt

The following files contain the pre-processed summarized experiment datasets for each case of study:
./LUAD/all/LUAD_Illumina_HiSeq_all.rda
./LUAD/paired/LUAD_Illumina_HiSeq_paired.rda
./LUSC/all/LUSC_Illumina_HiSeq_all.rda
./LUSC/paired/LUSC_Illumina_HiSeq_paired.rda

The following files contain the summarized experiment datasets upon filtering of tumor purity:
./LUAD/all/LUAD_Illumina_HiSeq_all_tumorPurity.rda
./LUAD/paired/LUAD_Illumina_HiSeq_paired_tumorPurity.rda
./LUAD/unpaired/LUAD_Illumina_HiSeq_unpaired_tumorPurity.rda
./LUSC/all/LUSC_Illumina_HiSeq_all_tumorPurity.rda
./LUSC/paired/LUAD_Illumina_HiSeq_paired_tumorPurity.rda
/LUSC/unpaired/LUSC_Illumina_HiSeq_unpaired_tumorPurity.rda

The following files contain the summarized experiment datasets upon the normalization steps:
./LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda
./LUAD/paired/LUAD_PreprocessedData_paired_tumorPurity.rda
./LUAD/unpaired/LUAD_PreprocessedData_unpaired_tumorPurity.rda
./LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda
./LUSC/paired/LUSC_PreprocessedData_paired_tumorPurity.rda
./LUSC/unpaired/LUSC_PreprocessedData_unpaired_tumorPurity.rda

