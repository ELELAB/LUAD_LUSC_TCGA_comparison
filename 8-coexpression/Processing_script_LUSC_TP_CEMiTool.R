source("elena/TCGAbiolinks_functions.R")
library(TCGAbiolinks)
library(SummarizedExperiment)
#source("http://bioconductor.org/biocLite.R")
#biocLite("pamr")

#---------------------------------------------------------------------------------
#                           LUSC-all
#------------------------------------------------------------------------------

SE_LUSC <- get(load("elena/LUSC_Illumina_HiSeq_TP.rda"))

#check the samples number
length(which(colData(SE_LUSC)$shortLetterCode =="TP"))
length(which(colData(SE_LUSC)$shortLetterCode =="NT"))

#remove discordant-LUSC samples
discordant.samples <- read.table("elena/discordant_samples.txt",col.names = "barcode")
SE_LUSC <- discordant_removing(discordant.samples,SE_LUSC)

#remove samples with low tumor purity
list <- TCGAtumor_purity(colnames(SE_LUSC),0,0,0,0,0.6)
length(list$pure_barcodes)
length(list$filtered)
SE_LUSC <- SE_LUSC[,colnames(SE_LUSC) %in% union(list$pure_barcodes,list$filtered)]
#save(SE_LUAD, file="LUAD_Illumina_HiSeq_TP_tumorPurity.rda")

#Preprocessing
dataPrep <- TCGAanalyze_Preprocessing(object = SE_LUSC, cor.cut = 0.6)


dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                      geneInfo = geneInfo,
                                      method = "gcContent")                

dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile", 
                                  qnt.cut =  0.25)

#save(dataFilt, file = "LUSC_PreprocessedData_TP_tumorPurity.rda")
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
write.table(modules, file="elena/CEMiTool_summary/modules_lusc.txt", sep=" ", quote = FALSE)
#first 10 hubs with highest connnectivity
hubs <- get_hubs(cem,10)
hubs_M1 <- as.data.frame(hubs[["M1"]])
hubs_M2 <- as.data.frame(hubs[["M2"]])
hubs_M3 <- as.data.frame(hubs[["M3"]])
hubs_M4 <- as.data.frame(hubs[["M4"]])
hubs_M5 <- as.data.frame(hubs[["M5"]])
hubs_M6 <- as.data.frame(hubs[["M6"]])
write.table(hubs_M1, file="elena/CEMiTool_summary/hubs_M1_lusc.txt",sep=" ", quote = FALSE)
write.table(hubs_M2, file="elena/CEMiTool_summary/hubs_M2_lusc.txt",sep=" ", quote = FALSE)
write.table(hubs_M3, file="elena/CEMiTool_summary/hubs_M3_lusc.txt",sep=" ", quote = FALSE)
write.table(hubs_M4, file="elena/CEMiTool_summary/hubs_M4_lusc.txt",sep=" ", quote = FALSE)
write.table(hubs_M5, file="elena/CEMiTool_summary/hubs_M5_lusc.txt",sep=" ", quote = FALSE)
write.table(hubs_M6, file="elena/CEMiTool_summary/hubs_M6_lusc.txt",sep=" ", quote = FALSE)
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
plots[7]

# read interactions
int_fname <- system.file("extdata", "interactions.tsv", package = "CEMiTool")
int_df <- read.delim(int_fname)
head(int_df)

# plot interactions
library(ggplot2)
interactions_data(cem) <- int_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[5]

# read interactions with i2d
i2d <- read.table("i2d_geneName_OK.txt")
i2d_df <- as.data.frame(i2d)


# plot interactions with i2d
library(ggplot2)
interactions_data(cem) <- i2d_df # add interactions
cem <- plot_interactions(cem) # generate plot
plots <- show_plot(cem, "interaction") # view the plot for the first module
plots[6]
