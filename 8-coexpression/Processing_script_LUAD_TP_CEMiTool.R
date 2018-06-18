#before running this script, the user needs to run the download_script_tp_LUAD.R script
#remember to setwd("~/8-coexpression")
source("TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(limma)
library(CEMiTool)
library(ggplot2)

#---------------------------------------------------------------------------------
#                           LUAD-all
#------------------------------------------------------------------------------

#if you want to reproduce exactly the analyses of this publication you can just start from here and load the SE object
#already provided with the folder (i.e. the download step can be skipped)
SE_LUAD <- get(load("LUAD_Illumina_HiSeq_TP.rda"))

#check the samples number - NT should be 0 - we only need TP for this analysis
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]
#save(SE_LUAD, file="LUAD_Illumina_HiSeq_TP_tumorPurity.rda")

#Preprocessing - should be coherent with what was done in the pre-processing for DEA
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

#save(dataFilt, file = "LUAD_PreprocessedData_TP_tumorPurity.rda")

#data normalised with the voom methodology
dataFilt_voom = voom(dataFilt)$E
dataFilt_ok <- as.data.frame(dataFilt_voom)

#here the coexpression analysis starts
cem <- cemitool(dataFilt_ok)
cem #to have basics information on the modules and selected genes
filter_expr(cem)

# inspect modules
nmodules(cem)
head(module_genes(cem))
modules<- module_genes(cem)
#verify that in the folder you have the CEMiTool_summary directory already created
write.table(modules, file="./CEMiTool_summary/modules_luad.txt", sep=" ", quote = FALSE)

#first 10 hubs with highest connnectivity
hubs <- get_hubs(cem,10)
hubs_M1 <- as.data.frame(hubs[["M1"]])
hubs_M2 <- as.data.frame(hubs[["M2"]])
hubs_M3 <- as.data.frame(hubs[["M3"]])
hubs_M4 <- as.data.frame(hubs[["M4"]])
write.table(hubs_M1, file="./CEMiTool_summary/hubs_M1_luad.txt",sep=" ", quote = FALSE)
write.table(hubs_M2, file="./CEMiTool_summary/hubs_M2_luad.txt",sep=" ", quote = FALSE)
write.table(hubs_M3, file="./CEMiTool_summary/hubs_M3_luad.txt",sep=" ", quote = FALSE)
write.table(hubs_M4, file="./CEMiTool_summary/hubs_M4_luad.txt",sep=" ", quote = FALSE)
summary <- mod_summary(cem)
generate_report(cem, output_format=c("pdf_document", "html_document"), force=TRUE)
write_files(cem, directory="./CEMiTool_summary/Tables", force=TRUE)

# plot gene expression within each module - in the example reported for module 1
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[1]

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)

# perform over representation analysis
cem <- mod_ora(cem, gmt_in)
# plot ora results - in the example for module 1
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[1]

# read interactions with the CEMiTool default 
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# plot interactions - in the example for module 1
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]

# read interactions with i2d 
i2d <- read.table("i2d_geneName.txt")
i2d_df <- as.data.frame(i2d)


# plot interactions with i2d - in the example for module 1
library(ggplot2)
interactions_data(cem) <- i2d_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[1]
