setwd("~/cross-validation/GSE33532/DEA")
library(WGCNA)
#load and install (if not) R packages 
source("../../install_packages_script.R")

# add normal samples to the data
data <- read.csv("../final_matrix.csv",row.names = 1)
NL1 <- read.table("NL_GSM835272.txt")
NL2 <- read.table('NL_GSM835277.txt')
NL3 <- read.table('NL_GSM835297.txt')
NL4 <- read.table("NL_GSM835302.txt")
NL5 <- read.table('NL_GSM835307.txt')
NL6 <- read.table('NL_GSM835312.txt')
NL7 <- read.table('NL_GSM835322.txt')
NL8 <- read.table('NL_GSM835337.txt')
NL9 <- read.table('NL_GSM835342.txt')
NL10 <- read.table('NL_GSM835352.txt')
NL11 <- read.table('NL_GSM835317.txt')
NL12 <- read.table('NL_GSM835327.txt')
NL13 <- read.table('NL_GSM835357.txt')
NL14 <- read.table('NL_GSM835362.txt')
data <- cbind(data,NL1,NL2,NL3,NL4,NL5,NL6,NL7,NL8,NL9,NL10,NL11,NL12,NL13,NL14)

rowID <- rownames(data)
dim(data)
if(file.exists("../conversion_file.csv")){
  conversion <- read.csv("../conversion_file.csv",row.names = 1)} else{
    conversion <- gconvert(rowID, organism = "hsapiens", target="AFFY_HG_U133_PLUS_2")}

gene_symbol <- conversion$name[match(rowID,tolower(conversion$alias))]
length(gene_symbol)
# remove all the non-converted probes
conversion_data <- na.omit(cbind(rowID,as.character(gene_symbol)))
new_data <- subset(data,rownames(data) %in% conversion_data[,"rowID"])
dim(new_data)

rowGroup <- conversion_data[,2]
rowID <- rownames(new_data)

# collapse
collapse.object=collapseRows(datET=new_data, rowGroup=rowGroup, rowID=rowID)

dat1Collapsed=data.frame(collapse.object$group2row, collapse.object$datETcollapsed)
# check if the gene names are unique
nrow(dat1Collapsed)
length(unique(dat1Collapsed$group))
length(unique(dat1Collapsed$selectedRowID))
# remove group and selectedRowID columns
dat1Collapsed <- dat1Collapsed[,- c(1,2)]
write.csv(dat1Collapsed,file="collapsed_data_LUAD_LUSC_Normal.csv")

#------------------------------------------------------------------------------------------
#                               DEA- limma
#------------------------------------------------------------------------------------------

# make condition factor (levels: NL, LUAD, LUSC)
condition <- vector()
index_nl <- grep("NL",colnames(dat1Collapsed))
index_luad <- grep("LUAD",colnames(dat1Collapsed))
index_lusc <- grep("LUSC",colnames(dat1Collapsed))
condition[index_nl] <- "NL"
condition[index_luad] <- "LUAD"
condition[index_lusc] <- "LUSC"
condition <- as.factor(condition)

#design matrix
design.matrix <- model.matrix(~0+condition)
colnames(design.matrix) <- c("LUAD","LUSC","NL")

# make contrasts: adenocarcinoma vs normal and squamous vs normal
contrast <- makeContrasts("LUAD-NL","LUSC-NL", levels= design.matrix)

fit <- lmFit(dat1Collapsed, design.matrix)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

LUADvsNL <- topTable(fit2, coef=1, adjust="fdr", number = nrow(dat1Collapsed))
LUSCvsNL <- topTable(fit2, coef=2, adjust="fdr", number = nrow(dat1Collapsed))

#add column which specifies if genes are up- down-regulated
make_DEA_table <- function(tt,FDR,LFC){
  index.up <- which(tt$logFC >= LFC & tt$adj.P.Val < FDR)
  index.down <- which(tt$logFC <= -LFC & tt$adj.P.Val < FDR)
  direction <- c()
  direction[index.up] <- "up"
  direction[index.down] <- "down"
  direction[!(1:nrow(tt) %in% union(index.up,index.down))] <- "no DE"
  tt <- cbind(tt,direction)
  return(tt)
}

LUADvsNL <- make_DEA_table(LUADvsNL,0.01,1)
# differentially expressed genes file
write.csv(LUADvsNL, "limma_DEA_LUADvsNL.csv", quote = FALSE)

#number of up-regulated genes cancer vs normal
print(paste0("the analysis of adenocarcinoma vs normal lung detected ",length(which(LUADvsNL$direction == "up"))," up-regulated genes"))
#number of down-regulated genes cancer vs normal
print(paste0("the analysis of adenocarcinoma vs normal lung detected ",length(which(LUADvsNL$direction== "down")), " down-regulated genes"))

up <- data.frame(rownames(LUADvsNL[LUADvsNL$direction == "up", ]))
down <- data.frame(rownames(LUADvsNL[LUADvsNL$direction == "down", ]))
write.table(up, "up_LUADvsNL.txt", sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down, "down_LUADvsNL.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



LUSCvsNL <- make_DEA_table(LUSCvsNL,0.01,1)
# differentially expressed genes file
write.csv(LUSCvsNL, "limma_DEA_LUSCvsNL.csv", quote = FALSE)

#number of up-regulated genes cancer vs normal
print(paste0("the analysis of squamous vs normal lung detected ",length(which(LUSCvsNL$direction == "up"))," up-regulated genes"))
#number of down-regulated genes cancer vs normal
print(paste0("the analysis of squamous vs normal lung detected ",length(which(LUSCvsNL$direction== "down")), " down-regulated genes"))

# up and down regulated genes
up <- data.frame(rownames(LUSCvsNL[LUSCvsNL$direction == "up", ]))
down <- data.frame(rownames(LUSCvsNL[LUSCvsNL$direction == "down", ]))
write.table(up, "up_LUSCvsNL.txt", sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down, "down_LUSCvsNL.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#-----------------------------------------------------------------------------------
#             compare the DEA results (LUADvsNT and LUSCvsNT)
#-----------------------------------------------------------------------------------

# up-regulated genes
up_luad <- read.table("up_LUADvsNL.txt")
up_lusc <- read.table("up_LUSCvsNL.txt")
overlap <- intersect(up_luad$V1,up_lusc$V1)

up_luad_only <- setdiff(up_luad$V1,overlap)
write.table(up_luad_only,"up_LUAD_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

up_lusc_only <- setdiff(up_lusc$V1,overlap)
write.table(up_lusc_only,"up_LUSC_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


# down-regulated genes
down_luad <- read.table("down_LUADvsNL.txt")
down_lusc <- read.table("down_LUSCvsNL.txt")
overlap <- intersect(down_luad$V1,down_lusc$V1)

down_luad_only <- setdiff(down_luad$V1,overlap)
write.table(down_luad_only,"down_LUAD_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

down_lusc_only <- setdiff(down_lusc$V1,overlap)
write.table(down_lusc_only,"down_LUSC_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

