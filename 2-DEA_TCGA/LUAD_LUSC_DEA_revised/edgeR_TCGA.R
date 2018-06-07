library(TCGAbiolinks)

#-------------------------------------------------------------------------------
#                     edgeR-TCGAbiolinks
#--------------------------------------------------------------------------------


edgeR_TCGA <- function(samplesNT,samplesTP,dataFilt,edgeR_TCGA_name,up_name,down_name){

dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,samplesNT],
                            mat2 = dataFilt[,samplesTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")  

# add information about the genes modulation
my.tags <- dataDEGs
index.up <- which(my.tags$logFC >= 1 & my.tags$FDR < 0.01)
index.down <- which(my.tags$logFC <= -1 & my.tags$FDR < 0.01)
direction <- c()
direction[index.up] <- "up"
direction[index.down] <- "down"
direction[!(1:nrow(my.tags) %in% union(index.up,index.down))] <- "no DE"
my.tags <- cbind(my.tags,direction)
write.csv(my.tags,edgeR_TCGA_name, quote = FALSE)

#number of up-regulated genes cancer vs normal
print(paste0("the analysis has detected ",length(which(my.tags$direction == "up"))," up-regulated genes"))
#number of down-regulated genes cancer vs normal
print(paste0("the analysis has detected ",length(which(my.tags$direction== "down")), " down-regulated genes"))

# up and down regulated genes
up <- data.frame(rownames(my.tags[my.tags$direction == "up", ]))
down <- data.frame(rownames(my.tags[my.tags$direction == "down", ]))
write.table(up, up_name, sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down,down_name, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
}


#----------------------------------------------------------------------------------
# LUSC paired
#--------------------------------------------------------------
dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/LUSC_PreprocessedData_woBatch-voom_paired.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/edgeR_TCGA/edgeR_TCGA_LUSC_paired.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/edgeR_TCGA/up_edgeR_TCGA_LUSC_paired.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/paired/edgeR_TCGA/down_edgeR_TCGA_LUSC_paired.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)

#----------------------------------------------------------------------------------
# LUSC paired after removing low tumor purity samples
#--------------------------------------------------------------
dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/LUSC_PreprocessedData_paired_TumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/edgeR_TCGA/edgeR_TCGA_LUSC_paired_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/edgeR_TCGA/up_edgeR_TCGA_LUSC_paired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/paired/edgeR_TCGA/down_edgeR_TCGA_LUSC_paired_tumorPurity.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)


#----------------------------------------------------------------------------
# LUAD paired
#--------------------------------------------------------------------------

dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/paired/LUAD_PreprocessedData_woBatch-voom_paired.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/paired/edgeR_TCGA/edgeR_TCGA_LUAD_paired.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/paired/edgeR_TCGA/up_edgeR_TCGA_LUAD_paired.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/paired/edgeR_TCGA/down_edgeR_TCGA_LUAD_paired.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)

#----------------------------------------------------------------------------
# LUAD paired after removing low tumor purity samples
#--------------------------------------------------------------------------

dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/LUAD_PreprocessedData_paired_TumorPurity.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/edgeR_TCGA/edgeR_TCGA_LUAD_paired_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/edgeR_TCGA/up_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/paired/edgeR_TCGA/down_edgeR_TCGA_LUAD_paired_tumorPurity.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUAD all
#--------------------------------------------------------------------------------

dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/LUAD_PreprocessedData_all.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/edgeR_TCGA/edgeR_TCGA_LUAD_all.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/edgeR_TCGA/up_edgeR_TCGA_LUAD_all.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/all/edgeR_TCGA/down_edgeR_TCGA_LUAD_all.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)

#--------------------------------------------------------------------------------
# LUAD all after removing low tumor purity samples
#--------------------------------------------------------------------------------
dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/LUAD_PreprocessedData_all_TumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/edgeR_TCGA/edgeR_TCGA_LUAD_all_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/edgeR_TCGA/up_edgeR_TCGA_LUAD_all_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/all/edgeR_TCGA/down_edgeR_TCGA_LUAD_all_tumorPurity.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)


#---------------------------------------------------------------------------------
# LUAD unpaired
#---------------------------------------------------------------------------------

dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/LUAD_PreprocessedData_unpaired.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/edgeR_TCGA/edgeR_TCGA_LUAD_unpaired.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/edgeR_TCGA/up_edgeR_TCGA_LUAD_unpaired.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD/unpaired/edgeR_TCGA/down_edgeR_TCGA_LUAD_unpaired.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)

#---------------------------------------------------------------------------------
# LUAD unpaired after removing low tumor purity samples
#---------------------------------------------------------------------------------

dataframe_LUAD <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/LUAD_PreprocessedData_unpaired_TumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUAD),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/edgeR_TCGA/edgeR_TCGA_LUAD_unpaired_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/edgeR_TCGA/up_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUAD_tss_tumorPurity/unpaired/edgeR_TCGA/down_edgeR_TCGA_LUAD_unpaired_tumorPurity.txt"
edgeR_TCGA(samplesNT,samplesTP,dataframe_LUAD,edgeR_TCGA_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUSC all
#------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/LUSC_PreprocessedData_all.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/edgeR_TCGA/edgeR_TCGA_LUSC_all.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/edgeR_TCGA/up_edgeR_TCGA_LUSC_all.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/all/edgeR_TCGA/down_edgeR_TCGA_LUSC_all.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)

#-------------------------------------------------------------------------------
# LUSC all after removing low tumor purity
#------------------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/LUSC_PreprocessedData_all_TumorPurity.rda"))
# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/edgeR_TCGA/edgeR_TCGA_LUSC_all_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/edgeR_TCGA/up_edgeR_TCGA_LUSC_all_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/all/edgeR_TCGA/down_edgeR_TCGA_LUSC_all_tumorPurity.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)

#----------------------------------------------------------------
# LUSC unpaired
#----------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/LUSC_PreprocessedData_unpaired.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/edgeR_TCGA/edgeR_TCGA_LUSC_unpaired.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/edgeR_TCGA/up_edgeR_TCGA_LUSC_unpaired.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC/unpaired/edgeR_TCGA/down_edgeR_TCGA_LUSC_unpaired.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)

#----------------------------------------------------------------
# LUSC unpaired after removing low tumor purity samples
#----------------------------------------------------------------

dataframe_LUSC <- get(load("/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/LUSC_PreprocessedData_unpaired_TumorPurity.rda"))

# selection of normal samples "NT"
samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("NT"))
# selection of tumor samples "TP"
samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataframe_LUSC),typesample = c("TP"))
edgeR_TCGA_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/edgeR_TCGA/edgeR_TCGA_LUSC_unpaired_tumorPurity.csv"
up_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/edgeR_TCGA/up_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"
down_name <- "/data/user/marta/pipeline/DE/new_LUAD-LUSC/LUSC_tss_tumorPurity/unpaired/edgeR_TCGA/down_edgeR_TCGA_LUSC_unpaired_tumorPurity.txt"

edgeR_TCGA(samplesNT,samplesTP,dataframe_LUSC,edgeR_TCGA_name,up_name,down_name)



