source("TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)


#---------------------------------------------------------------------------------
#                           LUAD-all
#------------------------------------------------------------------------------

SE_LUAD <- get(load("elena/LUAD_Illumina_HiSeq_TP.rda"))

#check the samples number
length(which(colData(SE_LUAD)$shortLetterCode =="TP"))
length(which(colData(SE_LUAD)$shortLetterCode =="NT"))

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUAD),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUAD <- SE_LUAD[,colnames(SE_LUAD) %in% union(list$pure_barcodes,list$filtered)]
#save(SE_LUAD, file="LUAD_Illumina_HiSeq_TP_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUAD, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

#save(dataFilt, file = "LUAD_PreprocessedData_TP_tumorPurity.rda")
#expression estimates with counts in less than 20% of cases
#dataFilt = dataFilt[apply(dataFilt,1,function(x) sum(x==0))<ncol(dataFilt)*0.8,]
#data normalised with the voom methodology
library(limma)
dataFilt_voom = voom(dataFilt)$E
#source("https://bioconductor.org/biocLite.R")
#biocLite("CEMiTool")
library("CEMiTool")
dataFilt_ok <- as.data.frame(dataFilt_voom)
cem <- cemitool(dataFilt_ok)
cem
filter_expr(cem)
# inspect modules
nmodules(cem)
head(module_genes(cem))
modules<- module_genes(cem)
write.table(modules, file="elena/CEMiTool_summary/modules_luad.txt", sep=" ", quote = FALSE)
hubs <- get_hubs(cem,10)
hubs_M1 <- as.data.frame(hubs[["M1"]])
hubs_M2 <- as.data.frame(hubs[["M2"]])
hubs_M3 <- as.data.frame(hubs[["M3"]])
hubs_M4 <- as.data.frame(hubs[["M4"]])
write.table(hubs_M1, file="elena/CEMiTool_summary/hubs_M1_luad.txt",sep=" ", quote = FALSE)
write.table(hubs_M2, file="elena/CEMiTool_summary/hubs_M2_luad.txt",sep=" ", quote = FALSE)
write.table(hubs_M3, file="elena/CEMiTool_summary/hubs_M3_luad.txt",sep=" ", quote = FALSE)
write.table(hubs_M4, file="elena/CEMiTool_summary/hubs_M4_luad.txt",sep=" ", quote = FALSE)
summary <- mod_summary(cem)
generate_report(cem, output_format=c("pdf_document", "html_document"), force=TRUE)
write_files(cem, directory="./Tables", force=TRUE)
# plot gene expression within each module
cem <- plot_profile(cem)
plots <- show_plot(cem, "profile")
plots[2]

# read GMT file
gmt_fname <- system.file("extdata", "pathways.gmt", package = "CEMiTool")
gmt_in <- read_gmt(gmt_fname)
# perform over representation analysis
cem <- mod_ora(cem, gmt_in)
# plot ora results
cem <- plot_ora(cem)
plots <- show_plot(cem, "ora")
plots[5]

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[4]

# read interactions with i2d
i2d <- read.table("i2d_geneName_OK.txt")
i2d_df <- as.data.frame(i2d)


# plot interactions with i2d
library(ggplot2)
interactions_data(cem) <- i2d_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[4]
