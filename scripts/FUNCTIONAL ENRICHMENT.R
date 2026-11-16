
library(gprofiler2)
library(ggplot2)

setwd("D:/RNA SEQ")

# Ensure output directories exist
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

run_gprofiler_enrichment <- function(gene_file, prefix, org = "hsapiens") {
  if (!file.exists(gene_file)) {
    stop(paste("File not found:", gene_file))
  }
  
  genes <- readLines(gene_file)
  genes <- trimws(genes[genes != ""])
  
  if (length(genes) == 0) {
    message(paste("No genes found in", gene_file, "- skipping."))
    return(NULL)
  }
  
  message(paste("Running g:Profiler on", length(genes), "genes for:", prefix))
  
  # Run enrichment analysis
  # Change org = "mmusculus" for Mouse, "rnorvegicus" for Rat, etc.
  gost_res <- gost(
    query = genes,
    organism = org,
    ordered_query = FALSE,
    multi_query = FALSE,
    significant = TRUE,
    exclude_iea = FALSE,
    measure_underrepresentation = FALSE,
    evcodes = TRUE, # Includes the exact matching gene symbols in output
    user_threshold = 0.05,
    correction_method = "fdr",
    sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC")
  )
  
  if (!is.null(gost_res$result) && nrow(gost_res$result) > 0) {
    res_df <- gost_res$result
    
    # Flatten list columns so write.csv exports cleanly
    res_df$parents <- sapply(res_df$parents, paste, collapse = ";")
    
    # Save full enrichment results table
    out_csv <- paste0("results/tables/gProfiler_", prefix, "_enrichment.csv")
    write.csv(res_df, out_csv, row.names = FALSE)
    
    # Plot top 15 enriched terms
    top_terms <- head(res_df[order(res_df$p_value), ], 15)
    
    out_pdf <- paste0("results/figures/gProfiler_", prefix, "_dotplot.pdf")
    pdf(out_pdf, width = 9.5, height = 6.5)
    print(
      ggplot(
        top_terms,
        aes(
          x = -log10(p_value),
          y = reorder(term_name, -log10(p_value)),
          size = intersection_size,
          color = source
        )
      ) +
        geom_point(alpha = 0.85) +
        theme_minimal(base_size = 12) +
        labs(
          title = paste("Top Enriched Terms -", prefix),
          subtitle = paste("Organism:", org, "| FDR < 0.05"),
          x = "-Log10 Adjusted P-value (FDR)",
          y = NULL,
          size = "Gene Count",
          color = "Database"
        )
    )
    dev.off()
    
    message(paste("Saved:", out_csv, "and", out_pdf))
  } else {
    message(paste("No significant terms found at FDR < 0.05 for:", prefix))
  }
}

# ================================
# Run for Up & Downregulated Lists
# ================================
run_gprofiler_enrichment("GO_UPREGULATED_GENES.txt", "Upregulated", org = "hsapiens")
run_gprofiler_enrichment("GO_DOWNREGULATED_GENES.txt", "Downregulated", org = "hsapiens")

