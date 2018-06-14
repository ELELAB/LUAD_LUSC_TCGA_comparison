#remember to setwd("~/5-mfuzz")
#unzip the .gz files in the LUAD/data and LUSC/data folders that are used as inputs
# Build consensus matrix (all genes expression and log2fc values in LUAD)
# load raw tables

library(TCGAbiolinks)
library(Mfuzz)

#extract the NT samples
data <- get(load("./LUAD/data/LUAD_Illumina_HiSeq_all.rda"))
NT_barcodes <- TCGAquery_SampleTypes(colnames(data),"NT")


query.luad2 <- GDCquery(project = "TCGA-LUAD",
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts",
                        barcode = NT_barcodes)


GDCdownload(query=query.luad2)

dataPrep1.luad <- GDCprepare(query = query.luad2, 
                             save = TRUE )

dataPrep<-TCGAanalyze_Preprocessing(object = dataPrep1.luad, 
                                    cor.cut = 0.6)

# write the preprocessed NT expression matrix
write.table(file="LUAD_NT_Preprocessed.tsv",dataPrep,col.names=T,row.names=T,sep='\t')

M1 <- read.table("./LUAD/data/data_LUAD_all_noDE.csv",header=T,sep=',',row.names=1,check.names=F)
M2 <- read.table("./LUAD/data/data_LUAD_all_DE.csv",header=T,sep=',',row.names=1,check.names=F) 
M3 <- rbind(M1,M2)
# reorder concatenated matrix
M3 <- M3[order(rownames(M3)),]
tail(M3[,1:5])

write.table(file='LUAD_all_consensus.csv',M3,sep=',',col.names=NA, row.names=TRUE)

# get the barcodes from the preprocessed NT expression matrix 
LUAD_NT_preprocessed <- read.table(file="LUAD_NT_Preprocessed.tsv",header=T,row.names=1,sep='\t',check.names=F)

NT_barcodes <- colnames(LUAD_NT_preprocessed)

consensus_matrix_LUAD <- read.table('LUAD_all_consensus.csv',sep=',',header=TRUE,row.names=1,check.names=F)

NT_index <- which(colnames(consensus_matrix_LUAD) %in% NT_barcodes)

# load clinical data to map barcode with stages
clinical_LUSC <- read.csv("./LUSC/data/TCGA-LUSC_clinical-Apr6.csv",sep=',',header=TRUE)
clinical_LUAD <- read.csv("./LUAD/data/TCGA-LUAD_clinical-Apr6.csv",sep=',',header=TRUE)
clinical <- rbind(clinical_LUSC, clinical_LUAD)

# extract first 3 elements of sample ID corresponding to barcode
barcode_vect <- c()
for (i in 1:length(colnames(consensus_matrix_LUAD)[-c(NT_index)])){
barcode_vect <- rbind(barcode_vect,c(colnames(consensus_matrix_LUAD)[-c(NT_index)][i],paste(as.matrix(unlist(strsplit(colnames(consensus_matrix_LUAD)[-c(NT_index)][i],"[-]")))[1:3],collapse='-')))
}
barcode_vect <- as.matrix(barcode_vect)
colnames(barcode_vect) <- c('sample','barcode')
sample_matrix_LUAD <- merge(barcode_vect,clinical[c('tumor_stage','bcr_patient_barcode')],by.x='barcode',by.y='bcr_patient_barcode',all.x=TRUE)
colnames(sample_matrix_LUAD) <- c("barcode","sample","stage")


barcode_vect_NT <- c()
for (i in 1:length(colnames(consensus_matrix_LUAD)[c(NT_index)])){
  barcode_vect_NT <- rbind(barcode_vect_NT,c(paste(as.matrix(unlist(strsplit(colnames(consensus_matrix_LUAD)[c(NT_index)][i],"[-]")))[1:3],collapse='-'),colnames(consensus_matrix_LUAD)[c(NT_index)][i]))
}
barcode_vect_NT <- cbind(barcode_vect_NT,rep("NT",dim(barcode_vect_NT)[1]))
colnames(barcode_vect_NT) <- c("barcode","sample","stage")

sample_matrix_LUAD <- rbind(sample_matrix_LUAD,barcode_vect_NT)

# collapse cancer substages (a + b) into 1 stage
sample_matrix_2 <- as.matrix(sample_matrix_LUAD)

sample_matrix_2[which(sample_matrix_2[,3] == 'stage ia'),3] <- 'stage i'
sample_matrix_2[which(sample_matrix_2[,3] == 'stage ib'),3] <- 'stage i'

sample_matrix_2[which(sample_matrix_2[,3] == 'stage iia'),3] <- 'stage ii'
sample_matrix_2[which(sample_matrix_2[,3] == 'stage iib'),3] <- 'stage ii'

sample_matrix_2[which(sample_matrix_2[,3] == 'stage iiia'),3] <- 'stage iii'
sample_matrix_2[which(sample_matrix_2[,3] == 'stage iiib'),3] <- 'stage iii'

table(sample_matrix_2[,3])

sample_matrix_LUAD <- sample_matrix_2
sample_matrix_LUAD <- sample_matrix_LUAD[order(sample_matrix_LUAD[,2]),]


# sample columns reordered in alphabetic order to match the order in sample_matrix
consensus_matrix_LUAD_ordered <- consensus_matrix_LUAD[,order(colnames(consensus_matrix_LUAD))]
sample_matrix_LUAD <- as.matrix(sample_matrix_LUAD)

#samples with 'not reported' stage
which(as.matrix(sample_matrix_LUAD[,3]) == 'not reported')


# taking out samples with not reported stage
consensus_matrix_LUAD_ordered <- consensus_matrix_LUAD_ordered[,-which(as.matrix(sample_matrix_LUAD[,3]) == 'not reported')]
sample_matrix_LUAD <- sample_matrix_LUAD[-which(as.matrix(sample_matrix_LUAD[,3]) == 'not reported'),]


dge_LUAD <- consensus_matrix_LUAD_ordered
getAverageExpSampleSet<- function(dataSE, set){
  
  set_samples <- which(sample_matrix_LUAD[,3] == set)
  
  #calculate average expression
  set_mean <- matrix(nrow=dim(dataSE)[1], ncol=1)
  if (length(set_samples) > 1){
  set_mean <- as.matrix(rowMeans(dataSE[,set_samples]))
  }else{
  if (length(set_samples) == 1){
  set_mean <- as.matrix(dataSE[,set_samples])}
  }
  return(set_mean)
}

stage_averages_LUAD <- matrix(nrow=dim(dge_LUAD)[1], ncol=5)                          
colnames(stage_averages_LUAD)<-c('NT','stage i','stage ii','stage iii','stage iv')

rownames(stage_averages_LUAD)<- rownames(dge_LUAD)

# computing mean gene expr value per tumor stage
for (i in colnames(stage_averages_LUAD)){
  this_average<-getAverageExpSampleSet(dge_LUAD, i)
  stage_averages_LUAD[,i]<- this_average}
head(stage_averages_LUAD)

write.table(file='stage_averages_LUAD.csv',stage_averages_LUAD,sep=',',col.names=NA, row.names=TRUE)
stage_averages_LUAD <- read.table('stage_averages_LUAD.csv',sep=',',header=TRUE,row.names=1,check.names=F)
stage_averages_LUAD <- as.matrix(stage_averages_LUAD)




exprSet_LUAD=ExpressionSet(assayData=stage_averages_LUAD)

exprSet_LUAD.s <- standardise(exprSet_LUAD)  #NB this step makes expression between genes COMPARABLE
dim(exprSet_LUAD.s)

m1_LUAD=mestimate(exprSet_LUAD.s) 
m1_LUAD

cl <- mfuzz(exprSet_LUAD.s,c=6,m=m1_LUAD)


mfuzz.plot(exprSet_LUAD.s,cl=cl,mfrow=c(3,2))


current_test <- "./LUAD/results/Mfuzz_plots/"
pdf(file=paste0(current_test,"LUAD_NT_4stages_6clusters",".pdf"))
lab <-c( 'NT','I', 'II', 'III', 'IV') 

mfuzz.plot(exprSet_LUAD.s,cl=cl,mfrow=c(3,2),min.mem=0.56, time.labels=lab, new.window = F)
dev.off()

cores<-acore(exprSet_LUAD.s,cl=cl,min.acore=0.56)

write.csv(file = "./LUAD/results/Mfuzz_clusters/LUAD_NT_4stages_cluster1.csv", x = cores[[1]], quote = FALSE, row.names=FALSE)
write.csv(file = "./LUAD/results/Mfuzz_clusters/LUAD_NT_4stages_cluster2.csv", x = cores[[2]], quote = FALSE, row.names=FALSE)
write.csv(file = "./LUAD/results/Mfuzz_clusters/LUAD_NT_4stages_cluster3.csv", x = cores[[3]], quote = FALSE, row.names=FALSE)
write.csv(file = "./LUAD/results/Mfuzz_clusters/LUAD_NT_4stages_cluster4.csv", x = cores[[4]], quote = FALSE, row.names=FALSE)
write.csv(file = "./LUAD/results/Mfuzz_clusters/LUAD_NT_4stages_cluster5.csv", x = cores[[5]], quote = FALSE, row.names=FALSE)
write.csv(file = "./LUAD/results/Mfuzz_clusters/LUAD_NT_4stages_cluster6.csv", x = cores[[6]], quote = FALSE, row.names=FALSE)


#make a list object with gene names for each cluster
cluster_data<-list()
# get total number of genes after filetering : std = 0.3, membership > 0.56
gene_sum=0

for ( i in 1:dim(summary(cores))[1]){
  print (paste( "Cluster", i, ":" , dim(cores[[i]])[1] ))
  this_cluster<-as.character(paste0("cluster",i))
  cluster_data[[this_cluster]]<-rownames(cores[[i]])
  gene_sum=gene_sum+dim(cores[[i]])[1]
}


gene_sum


 



