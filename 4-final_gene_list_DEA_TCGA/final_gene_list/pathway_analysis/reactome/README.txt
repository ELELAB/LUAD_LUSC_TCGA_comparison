
pathway_analysis.R script performs the pathway enrichment analysis using the ReactomePA R package. We converted the gene symbols into entrez IDs using the SummarizedExperiment object. We performed the analysis for up- and down-regulated genes, separately.

As output, we have a table (csv format) and a plot (png format) for each analysis including the significant pathways and the corrisponding genes.