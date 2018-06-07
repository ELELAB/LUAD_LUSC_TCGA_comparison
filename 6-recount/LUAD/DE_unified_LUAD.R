# created by Marta 10/01/2018

library(SummarizedExperiment)
library(TCGAbiolinks)
library(biomaRt)
source("/data/user/shared_projects/luad_lusc_2018/marta/TCGAbiolinks_functions.R")

load("gtex.RData")
dim(eset.gtex)
load("tcga_cancer_luad.RData")
dim(eset.tcga.cancer)

eset.all<-cbind(eset.tcga.cancer, eset.gtex)
dim(eset.all)
#remove ".number" in rownames(eset.all) otherwise the normalization step and bioMart
# do not work
rownames(eset.all)<-gsub("\\..*", "",rownames(eset.all))

# get list of all protein-coding genes
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
listOfGenes <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol", "gene_biotype"),filters = c("biotype"),values = list(biotype="protein_coding"), mart = mart)
head(listOfGenes)
dim(listOfGenes)

eset.all <- subset(eset.all,rownames(eset.all) %in% listOfGenes$ensembl_gene_id)
dim(eset.all)

#normalization and filtering
eset.all.N<-TCGAanalyze_Normalization(tabDF = eset.all,
                                      geneInfo = geneInfoHT,
                                      method = "gcContent")

dim(eset.all.N)
save(eset.all.N,file = "./unifiedData_normalized_luad.RData")


eset.all.N <- get(load("unifiedData_normalized_luad.RData"))

#convert into gene symbol
rse <- get(load("TCGA-LUADTranscriptome_ProfilingTue_Mar__6_17:45:53_2018.RData"))
rownames(eset.all.N) <- rowData(rse)[match(rownames(eset.all.N), rowData(rse)[,"ensembl_gene_id"]),"external_gene_name"]

eset.all.NF<- TCGAanalyze_Filtering(tabDF = eset.all.N,
                                    method = "quantile", 
                                    qnt.cut =  0.25)
dim(eset.all.NF)


save(eset.all.NF,file = "unifiedData_normalized_filtered_luad.RData")

eset.all.NF <- get(load("unifiedData_normalized_filtered_luad.RData"))
dim(eset.all.NF)


#----------------------------- DEA limma---------------------------------------

# preliminary steps

load("tcga_meta_luad.RData")
tcga.metadata.cancer <- subset(tcga.metadata,tcga.metadata$uid %in% colnames(eset.all.NF),select=c("uid","tss"))
tcga.metadata.cancer$condition <- as.factor("cancer")
load("gtex_meta_luad.RData")
gtex.metadata$condition <- as.factor("normal")
colnames(gtex.metadata)[2] <- "uid"
gtex.metadata <- gtex.metadata[,-1]

joint_IDs <- rbind(tcga.metadata.cancer,gtex.metadata)
rownames(joint_IDs) <- seq(1:nrow(joint_IDs))
condition <- joint_IDs$condition[match(colnames(eset.all.NF),as.character(joint_IDs$uid))]
tss <- as.factor(as.character(joint_IDs$tss[match(colnames(eset.all.NF),as.character(joint_IDs$uid))]))


# limma

limma_name <- "./limma/limma_unified_LUAD_tss.csv"
up_name <- "./limma/up_limma_unified_LUAD_tss.txt"
down_name <- "./limma/down_limma_unified_LUAD_tss.txt"
limma_unified(eset.all.NF,condition,tss,limma_name,up_name,down_name)
