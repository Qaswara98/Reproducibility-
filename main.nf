#!/usr/bin/env nextflow

params.inputFile = "$baseDir/data/input.txt"

process DataPreprocessing  {
    publishDir "$baseDir/results", mode: 'copy'

    input:
file inputFile


    output:
    file 'normCounts_res.CSV'

    script:
    """
    Rscript - <<END
    library(dplyr)
    library(tidyr)
    library(GEOquery) 
    library(DESeq2) 
    library(edgeR)
    library(openxlsx) 
    library(biomaRt)
    options(timeout=8000) 

    counts <- read.table("\${inputFile}", header=TRUE, row.names=1)
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    attributes <- c("ensembl_gene_id", "gene_biotype")
    filters <- "ensembl_gene_id"
    gene_ids <- rownames(counts)
    gene_data <- getBM(attributes=attributes, filters=filters, values=gene_ids, mart=mart)
    
    END
    """
}
workflow {
    DataPreprocessing(file(params.inputFile))
}
