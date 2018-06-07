

#' @title Differentially expression analysis (DEA) using edgeR package.
#' @description
#'    TCGAanalyze_DEA allows user to perform Differentially expression analysis (DEA),
#'    using edgeR package to identify differentially expressed genes (DEGs).
#'     It is possible to do a two-class analysis.
#'
#'     TCGAanalyze_DEA performs DEA using following functions from edgeR:
#'     \enumerate{
#'     \item edgeR::DGEList converts the count matrix into an edgeR object.
#'     \item edgeR::estimateCommonDisp each gene gets assigned the same dispersion estimate.
#'     \item edgeR::exactTest performs pair-wise tests for differential expression between two groups.
#'     \item edgeR::topTags takes the output from exactTest(), adjusts the raw p-values using the
#'     False Discovery Rate (FDR) correction, and returns the top differentially expressed genes.
#'     }
#' @param mat1 numeric matrix, each row represents a gene,
#' each column represents a sample with Cond1type
#' @param mat2 numeric matrix, each row represents a gene,
#' each column represents a sample with Cond2type
#' @param Cond1type a string containing the class label of the samples in mat1
#'  (e.g., control group)
#' @param Cond2type a string containing the class label of the samples in mat2
#' (e.g., case group)
#' @param method is 'glmLRT' (1) or 'exactTest' (2).
#' (1) Fit a negative binomial generalized log-linear model to
#' the read counts for each gene
#' (2) Compute genewise exact tests for differences in the means between
#' two groups of negative-binomially distributed counts.
#' @param  fdr.cut is a threshold to filter DEGs according their p-value corrected
#' @param logFC.cut is a threshold to filter DEGs according their logFC
#' @param elementsRatio is number of elements processed for second for time consumation estimation
#' @importFrom edgeR DGEList estimateCommonDisp exactTest topTags estimateGLMCommonDisp
#' estimateGLMTagwiseDisp glmFit glmLRT
#' @export
#' @examples
#' dataNorm <- TCGAbiolinks::TCGAanalyze_Normalization(dataBRCA, geneInfo)
#' dataFilt <- TCGAanalyze_Filtering(tabDF = dataBRCA, method = "quantile", qnt.cut =  0.25)
#' samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("NT"))
#' samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt), typesample = c("TP"))
#' dataDEGs <- TCGAanalyze_DEA(dataFilt[,samplesNT],
#'                       dataFilt[,samplesTP],"Normal", "Tumor")
#'
#' @return table with DEGs containing for each gene logFC, logCPM, pValue,and FDR
TCGAanalyze_DEA <- function(mat1,
                            mat2,
                            Cond1type,
                            Cond2type,
                            method = "exactTest",
                            fdr.cut = 1,
                            logFC.cut = 0,
                            elementsRatio = 30000) {
  
  TOC <- cbind(mat1,mat2)
  Cond1num <- ncol(mat1)
  Cond2num <- ncol(mat2)
  
  message("----------------------- DEA -------------------------------")
  message(message1 <- paste( "there are Cond1 type", Cond1type ,"in ",
                             Cond1num, "samples"))
  message(message2 <- paste( "there are Cond2 type", Cond2type ,"in ",
                             Cond2num, "samples"))
  message(message3 <- paste( "there are ", nrow(TOC) ,
                             "features as miRNA or genes "))
  
  timeEstimated <- format(ncol(TOC)*nrow(TOC)/elementsRatio,digits = 2)
  message(messageEstimation <- paste("I Need about ", timeEstimated,
                                     "seconds for this DEA. [Processing 30k elements /s]  "))
  
  # Reading in the data and creating a DGEList object
  colnames(TOC) <- paste0('s',1:ncol(TOC))
  #DGE <- DGEList(TOC,group=rep(c("Normal","Tumor"),c(NormalSample,
  #TumorSample)))
  
  if (method == "exactTest"){
    DGE <- edgeR::DGEList(TOC,group = rep(c(Cond1type,Cond2type),
                                          c(Cond1num,Cond2num)))
    # Analysis using common dispersion
    disp <- edgeR::estimateCommonDisp(DGE) # Estimating the common dispersion
    #tested <- exactTest(disp,pair=c("Normal","Tumor")) # Testing
    tested <- edgeR::exactTest(disp,pair = c(Cond1type,Cond2type)) # Testing
    # Results visualization
    logFC_table <- tested$table
    tableDEA <- edgeR::topTags(tested,n = nrow(tested$table))$table
    tableDEA <- tableDEA[tableDEA$FDR <= fdr.cut,]
    tableDEA <- tableDEA[abs(tableDEA$logFC) >= logFC.cut,]
  }
  
  if (method == "glmLRT"){
    tumorType <- factor(x =  rep(c(Cond1type,Cond2type),
                                 c(Cond1num,Cond2num)),
                        levels = c(Cond1type,Cond2type))
    design <- model.matrix(~tumorType)
    aDGEList <- edgeR::DGEList(counts = TOC, group = tumorType)
    aDGEList <- edgeR::estimateGLMCommonDisp(aDGEList, design)
    aDGEList <- edgeR::estimateGLMTagwiseDisp(aDGEList, design)
    aGlmFit <- edgeR::glmFit(aDGEList, design, dispersion = aDGEList$tagwise.dispersion,
                             prior.count.total=0)
    aGlmLRT <- edgeR::glmLRT(aGlmFit, coef = 2)
    
    tableDEA <- cbind(aGlmLRT$table, FDR = p.adjust(aGlmLRT$table$PValue, "fdr"))
    tableDEA <- tableDEA[tableDEA$FDR < fdr.cut,]
    tableDEA <- tableDEA[abs(tableDEA$logFC) > logFC.cut,]
  }
  if(all(grepl("ENSG",rownames(tableDEA)))) tableDEA <- cbind(tableDEA,map.ensg(genes = rownames(tableDEA))[,2:3])
  message("----------------------- END DEA -------------------------------")
  
  return(tableDEA)
  
}
