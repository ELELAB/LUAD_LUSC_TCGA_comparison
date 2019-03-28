#set the working directory

setwd("~/cross-validation/Bhattacharjee/DEA")

# load and install (if not) R packages 
source("../../install_packages_script.R")

#######
# gconvert

# convert probe set to gene names and collapse them 
# (multiple probe set match same gene)


data <- read.csv("../DatasetA_12600gene.csv")
#remove "gene description" column
data <- data[,-2]
rownames(data) <- data$probe.set
#remove probe set column
data <- data[,-1]
rowID <- rownames(data)
dim(data)

#check if conversion file exists (already created in cross-validation/Bhattacharjee/script.R)
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

# retain adenocarcinoma, squamous and normal samples
dat1Collapsed <- dat1Collapsed[,union(grep("NL",colnames(dat1Collapsed)),union(grep("AD",colnames(dat1Collapsed)),grep("SQ",colnames(dat1Collapsed))))]
write.csv(dat1Collapsed,"collapsed_data_AD_SQ_NL.csv")

#check number of samples
length(grep("NL",colnames(dat1Collapsed)))
length(grep("AD",colnames(dat1Collapsed)))
length(grep("SQ",colnames(dat1Collapsed)))

#------------------------------------------------------------------------------------------
#                               DEA- limma
#------------------------------------------------------------------------------------------

# make condition factor (levels: NL, AD, SQ)
condition <- vector()
index_nl <- grep("NL",colnames(dat1Collapsed))
index_ad <- grep("AD",colnames(dat1Collapsed))
index_sq <- grep("SQ",colnames(dat1Collapsed))
condition[index_nl] <- "NL"
condition[index_ad] <- "AD"
condition[index_sq] <- "SQ"
condition <- as.factor(condition)

#design matrix
design.matrix <- model.matrix(~0+condition)
colnames(design.matrix) <- c("AD","NL","SQ")

# make contrasts: adenocarcinoma vs normal and squamous vs normal
contrast <- makeContrasts("AD-NL","SQ-NL", levels= design.matrix)

fit <- lmFit(dat1Collapsed, design.matrix)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

ADvsNL <- topTable(fit2, coef=1, adjust="fdr", number = nrow(dat1Collapsed))
SQvsNL <- topTable(fit2, coef=2, adjust="fdr", number = nrow(dat1Collapsed))

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

ADvsNL <- make_DEA_table(ADvsNL,0.01,1)
# differentially expressed genes file
write.csv(ADvsNL, "limma_DEA_ADvsNL.csv", quote = FALSE)

#number of up-regulated genes cancer vs normal
print(paste0("the analysis of adenocarcinoma vs normal lung detected ",length(which(ADvsNL$direction == "up"))," up-regulated genes"))
#number of down-regulated genes cancer vs normal
print(paste0("the analysis of adenocarcinoma vs normal lung detected ",length(which(ADvsNL$direction== "down")), " down-regulated genes"))

up <- data.frame(rownames(ADvsNL[ADvsNL$direction == "up", ]))
down <- data.frame(rownames(ADvsNL[ADvsNL$direction == "down", ]))
write.table(up, "up_ADvsNL.txt", sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down, "down_ADvsNL.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)



SQvsNL <- make_DEA_table(SQvsNL,0.01,1)
# differentially expressed genes file
write.csv(SQvsNL, "limma_DEA_SQvsNL.csv", quote = FALSE)

#number of up-regulated genes cancer vs normal
print(paste0("the analysis of squamous vs normal lung detected ",length(which(SQvsNL$direction == "up"))," up-regulated genes"))
#number of down-regulated genes cancer vs normal
print(paste0("the analysis of squamous vs normal lung detected ",length(which(SQvsNL$direction== "down")), " down-regulated genes"))

# up and down regulated genes
up <- data.frame(rownames(SQvsNL[SQvsNL$direction == "up", ]))
down <- data.frame(rownames(SQvsNL[SQvsNL$direction == "down", ]))
write.table(up, "up_SQvsNL.txt", sep = "\t",col.names = FALSE,row.names = FALSE, quote = FALSE)
write.table(down, "down_SQvsNL.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#-----------------------------------------------------------------------------------
#             compare the DEA results (ADvsNT and SQvsNT)
#-----------------------------------------------------------------------------------

# up-regulated genes
up_luad <- read.table("up_ADvsNL.txt")
up_lusc <- read.table("up_SQvsNL.txt")
overlap <- intersect(up_luad$V1,up_lusc$V1)

up_luad_only <- setdiff(up_luad$V1,overlap)
write.table(up_luad_only,"up_LUAD_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

up_lusc_only <- setdiff(up_lusc$V1,overlap)
write.table(up_lusc_only,"up_LUSC_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


# down-regulated genes
down_luad <- read.table("down_ADvsNL.txt")
down_lusc <- read.table("down_SQvsNL.txt")
overlap <- intersect(down_luad$V1,down_lusc$V1)

down_luad_only <- setdiff(down_luad$V1,overlap)
write.table(down_luad_only,"down_LUAD_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

down_lusc_only <- setdiff(down_lusc$V1,overlap)
write.table(down_lusc_only,"down_LUSC_only.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
