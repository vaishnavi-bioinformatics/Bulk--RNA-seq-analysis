library(DESeq2)
library(ggplot2)
library(pheatmap)


 # ================================
 # Load counts matrix
 # ================================
 counts <- read.csv("gene_count_matrix.csv", row.names = 1, check.names = FALSE)
 counts <- as.matrix(counts)
 # ================================
 # Load metadata
 # ================================
 sample_info <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)
 colnames(sample_info) <- trimws(colnames(sample_info))
 rownames(sample_info) <- sample_info$sample
 sample_info$condition <- as.factor(sample_info$condition)
 sample_info$condition <- relevel(sample_info$condition, ref = "control")
 # ================================
 # Ensure sample matching
 # ================================
 if (!all(colnames(counts) %in% rownames(sample_info))) {
   stop("Sample names in counts and metadata do not match!")
 }
 
 # Reorder metadata
 sample_info <- sample_info[colnames(counts), ]
 
 # ================================
 # Create DESeq2 dataset
 # ================================
 dds <- DESeqDataSetFromMatrix(
   countData = counts,
   colData = sample_info,
   design = ~ condition
 )
 
 # ================================
 # Run DESeq2
 # ================================
 dds <- DESeq(dds)
 # ================================
 # Extract results
 # ================================
 res <- results(dds)
 res_df <- as.data.frame(res)
 res_df$gene <- rownames(res_df)
 
 
 # ==================
 # Normalized counts 
 # ==================
 norm_counts <- counts(dds, normalized = TRUE)
 
 # ====================
 # Convert to dataframe
 # ====================
 norm_df <- as.data.frame(norm_counts)
 norm_df$gene <- rownames(norm_df)
 
 # ================================
 # Calculate mean expression
 # ================================
 control_samples <- rownames(sample_info[sample_info$condition == "control", ])
 case_samples <- rownames(sample_info[sample_info$condition == "diseased", ])
 
 mean_df <- data.frame(
   gene = rownames(norm_counts),
   mean_control = rowMeans(norm_counts[, control_samples, drop = FALSE]),
   mean_case = rowMeans(norm_counts[, case_samples, drop = FALSE])
 )
 
 # ================================
 # Merge all together
 # ================================
 final_df <- merge(res_df, mean_df, by = "gene")
 final_df <- merge(final_df, norm_df, by = "gene")
 
 # ================================
 # Save output
 # ================================
 write.csv(final_df, "DESeq2_results_with_means_and_samples.csv", row.names = FALSE)
 
 # Preview
 head(final_df)
 
 # ============================================================
 # 11. Significant features
 # ============================================================
 
 significant_DEGs <- res_df[
   !is.na(res_df$padj) &
     res_df$padj < 0.05 &
     abs(res_df$log2FoldChange) > 1,
 ]
 
 write.csv(
   significant_DEGs,
   file.path("Significant_DEGs.csv"),
   row.names = FALSE
 )
 
 
 # ============================================================
 # 12. Upregulated features
 # ============================================================
 
 upregulated <- res_df[
   !is.na(res_df$padj) &
     res_df$padj < 0.05 &
     res_df$log2FoldChange > 1,
 ]
 
 write.csv(
   upregulated,
   file.path("Upregulated_features.csv"),
   row.names = FALSE
 )
 
 # ============================================================
 # 13. Downregulated features
 # ============================================================
 
 downregulated <- res_df[
   !is.na(res_df$padj) &
     res_df$padj < 0.05 &
     res_df$log2FoldChange < -1,
 ]
 
 write.csv(
   downregulated,
   file.path("Downregulated_features.csv"),
   row.names = FALSE
 )
 
 print("Number of significant features:")
 print(nrow(significant_DEGs))
 
 print("Number upregulated:")
 print(nrow(upregulated))
 
 print("Number downregulated:")
 print(nrow(downregulated))
 
 up_gene_list_clean <- sub("\\|.*","",upregulated$gene)
 up_gene_list_clean <- unique(up_gene_list_clean)
 
 down_gene_list_clean <- sub("\\|.*","", downregulated$gene)
 down_gene_list_clean <- unique(down_gene_list_clean)
 
 write.table(
   up_gene_list_clean,
   "GO_UPREGULATED_GENES.txt",
   row.names = FALSE,
   col.names = FALSE,
   quote = FALSE 
 )
 
 write.table(
   down_gene_list_clean,
   "GO_DOWNREGULATED_GENES.txt",
   row.names = FALSE,
   col.names = FALSE,
   quote = FALSE 
 )
 
 # ================================
 # rlog transformation
 # ================================
 rld <- rlog(dds, blind = TRUE)
 
 # ================================
 # PCA Plot
 # ================================
 pcaData <- plotPCA(
   rld,
   intgroup = "condition",
   returnData = TRUE
 )
 
 percentVar <- round(
   100 * attr(pcaData, "percentVar")
 )
 
 pdf(
   file.path("PCA_plot.pdf"),
   width = 7,
   height = 6
 )
 
 print(
   ggplot(
     pcaData,
     aes(
       x = PC1,
       y = PC2,
       color = condition,
       label = name
     )
   ) +
     geom_point(size = 4) +
     geom_text(
       vjust = -1,
       size = 3
     ) +
     xlab(
       paste0(
         "PC1: ",
         percentVar[1],
         "% variance"
       )
     ) +
     ylab(
       paste0(
         "PC2: ",
         percentVar[2],
         "% variance"
       )
     ) +
     theme_minimal(base_size = 14) +
     labs(
       color = "Condition"
     )
 )
 
 dev.off()
 # ================================
 # Volcano Plot
 # ================================
 volcano_df <- res_df[
   !is.na(res_df$padj),
 ]
 
 volcano_df$significance <- "Not Significant"
 
 volcano_df$significance[
   volcano_df$padj < 0.05 &
     volcano_df$log2FoldChange > 1
 ] <- "Upregulated"
 
 volcano_df$significance[
   volcano_df$padj < 0.05 &
     volcano_df$log2FoldChange < -1
 ] <- "Downregulated"
 
 pdf(
   file.path("Volcano_plot.pdf"),
   width = 8,
   height = 6
 )
 
 print(
   ggplot(
     volcano_df,
     aes(
       x = log2FoldChange,
       y = -log10(padj),
       color = significance
     )
   ) +
     geom_point(
       alpha = 0.7,
       size = 1.5
     ) +
     scale_color_manual(
       values = c(
         "Upregulated" = "red",
         "Downregulated" = "blue",
         "Not Significant" = "grey"
       )
     ) +
     geom_vline(
       xintercept = c(-1, 1),
       linetype = "dashed"
     ) +
     geom_hline(
       yintercept = -log10(0.05),
       linetype = "dashed"
     ) +
     theme_minimal(base_size = 14) +
     labs(
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value"
     )
 )
 
 dev.off()
 
 
 # ================================
 # Heatmap of Top 50 DE Genes
 # ================================

  res_ordered <- res[
   order(res$padj, na.last = NA),
 ]
 
 topGenes <- head(
   rownames(res_ordered),
   min(50, nrow(res_ordered))
 )
 
 print("Top features:")
 print(topGenes)
 
 if (length(topGenes) > 1) {
   
   mat <- assay(rld)[
     topGenes,
     ,
     drop = FALSE
   ]
   
   # Z-score by feature
   mat <- t(
     scale(
       t(mat)
     )
   )
   
   annotation_col <- data.frame(
     condition = sample_info$condition
   )
   
   rownames(annotation_col) <- rownames(sample_info)
   
   pdf(
     file.path(
       "Heatmap_Top50_DE_features.pdf"
     ),
     width = 8,
     height = 10
   )
   
   pheatmap(
     mat,
     annotation_col = annotation_col,
     show_rownames = TRUE,
     fontsize_row = 6,
    main = "Top 50 Differentially expressed genes"
   )
   
   dev.off()
 }
 
 

 