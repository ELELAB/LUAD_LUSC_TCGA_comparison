source('MCPcounter_functions.R')
source('TCGAbiolinks_functions.R')


############################################## Monocytic lineage ###################################

# luad

luad <- get(load('LUAD_PreprocessedData_all_tumorPurity.rda'))
estimation_luad <- MCPcounter.estimate(luad,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_luad, 'Monocytic lineage'))
edgeR_tss_immune(luad,'./Monocytic/edgeR_LUAD_Monocytic.csv','./Monocytic/up_edgeR_LUAD_Monocytic.txt',
                 './Monocytic/down_edgeR_LUAD_Monocytic.txt',Tcells)

### lusc
lusc <- get(load('LUSC_PreprocessedData_all_tumorPurity.rda'))
estimation_lusc <- MCPcounter.estimate(lusc,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_lusc, 'Monocytic lineage'))
edgeR_tss_immune(lusc,'./Monocytic/edgeR_LUSC_Monocytic.csv','./Monocytic/up_edgeR_LUSC_Monocytic.txt',
                 './Monocytic/down_edgeR_LUSC_Monocytic.txt',Tcells)


############################################## T cells ###################################

# luad

luad <- get(load('LUAD_PreprocessedData_all_tumorPurity.rda'))
estimation_luad <- MCPcounter.estimate(luad,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_luad, 'T cells'))
edgeR_tss_immune(luad,'./Tcells/edgeR_LUAD_Tcells.csv','./Tcells/up_edgeR_LUAD_Tcells.txt',
                 './Tcells/down_edgeR_LUAD_Tcells.txt',Tcells)

### lusc
lusc <- get(load('LUSC_PreprocessedData_all_tumorPurity.rda'))
estimation_lusc <- MCPcounter.estimate(lusc,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_lusc, 'T cells'))
edgeR_tss_immune(lusc,'./Tcells/edgeR_LUSC_Tcells.csv','./Tcells/up_edgeR_LUSC_Tcells.txt',
                 './Tcells/down_edgeR_LUSC_Tcells.txt',Tcells)


########################################## Cytotoxic lymphocytes ########################
# luad

luad <- get(load('../1-download_preprocessing/LUAD/all/LUAD_PreprocessedData_all_tumorPurity.rda'))
estimation_luad <- MCPcounter.estimate(luad,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_luad, 'Cytotoxic lymphocytes'))
edgeR_tss_immune(luad,'./Lymphocytes/edgeR_LUAD_lymphocytes.csv','./Lymphocytes/up_edgeR_LUAD_lymphocytes.txt',
                 './Lymphocytes/down_edgeR_LUAD_lymphocytes.txt',Tcells)

#lusc
lusc <- get(load('../1-download_preprocessing/LUSC/all/LUSC_PreprocessedData_all_tumorPurity.rda'))
estimation_lusc <- MCPcounter.estimate(lusc,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_lusc, 'Cytotoxic lymphocytes'))
edgeR_tss_immune(lusc,'./Lymphocytes/edgeR_LUSC_lymphocytes.csv','./Lymphocytes/up_edgeR_LUSC_lymphocytes.txt',
                 './Lymphocytes/down_edgeR_LUSC_lymphocytes.txt',Tcells)

########################################## B lineage ########################
# luad

luad <- get(load('LUAD_PreprocessedData_all_tumorPurity.rda'))
estimation_luad <- MCPcounter.estimate(luad,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_luad, 'B lineage'))
edgeR_tss_immune(luad,'./B_lineage/edgeR_LUAD_B_lineage.csv','./B_lineage/up_edgeR_LUAD_B_lineage.txt',
                 './B_lineage/down_edgeR_LUAD_B_lineage.txt',Tcells)

#lusc
lusc <- get(load('LUSC_PreprocessedData_all_tumorPurity.rda'))
estimation_lusc <- MCPcounter.estimate(lusc,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_lusc, 'B lineage'))
edgeR_tss_immune(lusc,'./B_lineage/edgeR_LUSC_B_lineage.csv','./B_lineage/up_edgeR_LUSC_B_lineage.txt',
                 './B_lineage/down_edgeR_LUSC_B_lineage.txt',Tcells)

##########################################   Neutrophils ########################
# luad

luad <- get(load('LUAD_PreprocessedData_all_tumorPurity.rda'))
estimation_luad <- MCPcounter.estimate(luad,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_luad, 'Neutrophils'))
edgeR_tss_immune(luad,'./Neutrophils/edgeR_LUAD_Neutrophils.csv','./Neutrophils/up_edgeR_LUAD_Neutrophils.txt',
                 './Neutrophils/down_edgeR_LUAD_Neutrophils.txt',Tcells)

#lusc
lusc <- get(load('LUSC_PreprocessedData_all_tumorPurity.rda'))
estimation_lusc <- MCPcounter.estimate(lusc,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_lusc, 'Neutrophils'))
edgeR_tss_immune(lusc,'./Neutrophils/edgeR_LUSC_Neutrophils.csv','./Neutrophils/up_edgeR_LUSC_Neutrophils.txt',
                 './Neutrophils/down_edgeR_LUSC_Neutrophils.txt',Tcells)

##########################################  Myeloid dendritic cells ########################
# luad

luad <- get(load('LUAD_PreprocessedData_all_tumorPurity.rda'))
estimation_luad <- MCPcounter.estimate(luad,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_luad, 'Myeloid dendritic cells'))
edgeR_tss_immune(luad,'./Myeloid_cells/edgeR_LUAD_Myeloid_cells.csv','./Myeloid_cells/up_edgeR_LUAD_Myeloid_cells.txt',
                 './Myeloid_cells/down_edgeR_LUAD_Myeloid_cells.txt',Tcells)

#lusc
lusc <- get(load('LUSC_PreprocessedData_all_tumorPurity.rda'))
estimation_lusc <- MCPcounter.estimate(lusc,featuresType="HUGO_symbols")
Tcells <- as.factor(abundanceRange(estimation_lusc, 'Myeloid dendritic cells'))
edgeR_tss_immune(lusc,'./Myeloid_cells/edgeR_LUSC_Myeloid_cells.csv','./Myeloid_cells/up_edgeR_LUSC_Myeloid_cells.txt',
                 './Myeloid_cells/down_edgeR_LUSC_Myeloid_cells.txt',Tcells)
