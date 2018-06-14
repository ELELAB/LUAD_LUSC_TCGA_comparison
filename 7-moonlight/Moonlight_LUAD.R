#remember to setwd("~/7-moonlight")
library("MoonlightR")
#####prediction TSG and OGs######

dataFilt <- getDataTCGA(cancerType = "LUAD", 
                        dataType = "Gene expression",
                        directory = "data",
                        nSample= 59)


dataDEGs <- DPA(dataFilt = dataFilt,
                dataType = "Gene expression", fdr.cut = 0.01,
                logFC.cut = 2, diffmean.cut = 0.25)

dataFEA <- FEA(DEGsmatrix = dataDEGs)

dataGRN <- GRN(TFs = rownames(dataDEGs)[1:1485], 
               DEGsmatrix = dataDEGs,
               DiffGenes = TRUE,
               normCounts = dataFilt)

dataURA <- URA(dataGRN = dataGRN, 
               DEGsmatrix = dataDEGs, 
               BPname = c("apoptosis",
                          "proliferation of cells"))

dataDual <- PRA(dataURA = dataURA, 
                BPname = c("apoptosis",
                           "proliferation of cells"),
                thres.role = 0)
CancerGenes <- list("TSG"=names(dataDual$TSG), "OCG"=names(dataDual$OCG))
#plotURA(dataURA = dataURA[c(names(dataDual$TSG), names(dataDual$OCG)),, drop = FALSE], additionalFilename = "_LUAD")
TSG <- as.data.frame(CancerGenes[["TSG"]])
OCG <- as.data.frame(CancerGenes[["OCG"]])
write.table(TSG, file ="Moonlight/TSG_LUAD.txt", sep=" ", quote = FALSE)
write.table(OCG, file ="Moonlight/OCG_LUAD.txt", sep=" ", quote = FALSE)
###dual roles across stages ####
listMoonlight <- NULL
for (i in 1:4){
  dataDual <- moonlight(cancerType = "LUAD", 
                        dataType = "Gene expression",
                        directory = "data",
                        nSample = 26,
                        nTF = 1485,
                        DiffGenes = TRUE,
                        BPname = c("apoptosis","proliferation of cells"),
                        stage = i)
  listMoonlight <- c(listMoonlight, list(dataDual))
  save(dataDual, file = paste0("dataDual_stage",as.roman(i), ".Rdata"))
}
names(listMoonlight) <- c("stage1", "stage2", "stage3", "stage4")
save(listMoonlight, file = paste0("listMoonlight_LUAD_stages.Rdata"))
plotCircos(listMoonlight = listMoonlight, additionalFilename = "_LUAD_stages")
#none MG found

