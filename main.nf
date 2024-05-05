#!/usr/bin/env nextflow

params.inputFile = "$baseDir/data/input.txt"

process DataPreprocessing  {
    publishDir "$baseDir/results", mode: 'copy'

    input:
    file inputFile from file(params.inputFile)

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
    gene_counts <- table(gene_data$gene_biotype)
    filter_biotypes <- c("artifact","TEC","pseudogene", "processed_pseudogene", "unprocessed_pseudogene", "TR_C_gene","sRNA","ribozyme","TR_D_gene",
                     "transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene", "vault_RNA","miRNA","scaRNA","IG_D_gene",
                     "unitary_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", "lncRNA","IG_V_gene","vault_RNA","TR_V_gene",
                     "IG_V_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "translated_processed_pseudogene","scRNA","IG_C_gene","IG_J_gene","TR_J_gene",
                     "rRNA_pseudogene", "rRNA", "Mt_rRNA", "Mt_tRNA", "snRNA", "snoRNA", "misc_RNA","transcribed_unitary_pseudogene","IG_pseudogene")
    gene_biotypes <- gene_data$gene_biotype
    filtered_gene_data <- gene_data[!gene_biotypes %in% filter_biotypes, ]
    counts_filtered <- counts[rownames(counts) %in% filtered_gene_data$ensembl_gene_id, ]
    remaining_biotypes <- unique(filtered_gene_data$gene_biotype)
    gse <- getGEO("GSE216738", GSEMatrix = TRUE, getGPL = FALSE)
    metadata <- pData(gse[[1]])
    metadata.subset <- metadata[, c(1, 46, 47, 48)]
    colnames(metadata.subset) <- c("samples", "eln_group", "Lncrna_score", "tissue")
    counts_filtered <- counts_filtered[, !(names(counts_filtered) %in% c("gene_name", "gene_biotype"))]
    rownames(metadata.subset) <- metadata.subset$samples
    counts_filtered_sample_names <- colnames(counts_filtered)
    metadata.subset_aligned <- metadata.subset[counts_filtered_sample_names, ]
    most_frequent_eln_group <- names(sort(table(metadata.subset_aligned$eln_group), decreasing = TRUE))[1]
    metadata.subset_aligned$eln_group[is.na(metadata.subset_aligned$eln_group)] <- most_frequent_eln_group
    most_frequent_Lncrna_score <- names(sort(table(metadata.subset_aligned$Lncrna_score), decreasing = TRUE))[1]
    metadata.subset_aligned$Lncrna_score[is.na(metadata.subset_aligned$Lncrna_score)] <- most_frequent_Lncrna_score
    meanLog2CPM <- rowMeans(log2(cpm(counts_filtered) + 1))
    counts_filtered <- counts_filtered[meanLog2CPM > 1, ]
    dds <- DESeqDataSetFromMatrix(countData = counts_filtered, colData = metadata.subset_aligned, design = ~ eln_group)  
    normCounts <- vst(dds)
    countData <- assay(normCounts)
    write.csv(countData, file = "normCounts_res.CSV")
    END
    """
}
workflow {
    DataPreprocessing (params.inputFile)
}
