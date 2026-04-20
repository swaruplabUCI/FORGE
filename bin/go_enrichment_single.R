#!/usr/bin/env Rscript
# GO Enrichment Analysis for single DE results file

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript go_enrichment_single.R <input_file>")
}

INPUT_FILE <- args[1]

suppressPackageStartupMessages({
  library(enrichR)
  library(tidyverse)
})

# Configuration
P_VALUE_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 0
TOP_N_ENRICHMENT <- 10
ENRICHMENT_PLOT_WIDTH <- 10
ENRICHMENT_PLOT_HEIGHT <- 8

# EnrichR databases (species-dependent - add logic here)
ENRICHR_DBS <- c(
  'GO_Biological_Process_2023',
  'GO_Cellular_Component_2023',
  'GO_Molecular_Function_2023',
  'KEGG_2021_Human',
  'WikiPathway_2023_Human'
)

# Parse filename
filename <- basename(INPUT_FILE)
base_name <- gsub("_DEResults\\.csv$", "", filename)

cat("Processing:", filename, "\n")

# Create output directory
output_dir <- dirname(INPUT_FILE)

# Load DE results
de_data <- read_csv(INPUT_FILE, show_col_types = FALSE)

# Filter significant genes
sig_genes <- de_data %>%
  filter(p_val_adj < P_VALUE_CUTOFF)

up_genes <- sig_genes %>%
  filter(avg_log2FC > LOG2FC_CUTOFF) %>%
  pull(gene) %>%
  na.omit() %>%
  unique()

down_genes <- sig_genes %>%
  filter(avg_log2FC < -LOG2FC_CUTOFF) %>%
  pull(gene) %>%
  na.omit() %>%
  unique()

cat("  Found", length(up_genes), "upregulated and", 
    length(down_genes), "downregulated genes\n")

# Text wrapping function
wrapText <- function(x, len) {
  sapply(x, function(y) paste(strwrap(y, len), collapse = "\n"), USE.NAMES = FALSE)
}

# Function to run enrichR analysis with retry on rate limit
run_enrichr_analysis <- function(gene_list, output_prefix, title_prefix, max_retries = 3) {
  if(length(gene_list) == 0) {
    cat("No genes for", title_prefix, "- skipping\n")
    return(NULL)
  }

  # Run enrichR with retry on rate limit (429)
  enriched <- NULL
  for (attempt in 1:max_retries) {
    enriched <- tryCatch(enrichr(gene_list, ENRICHR_DBS), error = function(e) NULL)
    # Check if result looks valid (has expected columns)
    if (!is.null(enriched) && length(enriched) > 0 &&
        "Adjusted.P.value" %in% colnames(enriched[[1]])) {
      break
    }
    cat(sprintf("  EnrichR attempt %d/%d returned malformed results, retrying in %ds...\n",
                attempt, max_retries, attempt * 30))
    Sys.sleep(attempt * 30)
  }

  # Final validation — if still malformed, write summary and exit gracefully
  if (is.null(enriched) || length(enriched) == 0 ||
      !"Adjusted.P.value" %in% colnames(enriched[[1]])) {
    cat("  WARNING: EnrichR API returned invalid results after", max_retries,
        "attempts — skipping", title_prefix, "\n")
    writeLines(paste("EnrichR failed for", title_prefix, "- API rate limit or unavailable"),
               file.path(output_dir, paste0(output_prefix, "_summary.txt")))
    return(NULL)
  }

  # Save results
  for(i in seq_along(ENRICHR_DBS)) {
    if(nrow(enriched[[i]]) > 0 && "Adjusted.P.value" %in% colnames(enriched[[i]])) {
      significant_results <- enriched[[i]] %>%
        filter(Adjusted.P.value < P_VALUE_CUTOFF)

      if(nrow(significant_results) > 0) {
        db_name <- gsub("[^a-zA-Z0-9]", "_", ENRICHR_DBS[i])
        csv_file <- file.path(output_dir, paste0(output_prefix, "_", db_name, ".csv"))
        write.csv(significant_results, csv_file, row.names = FALSE)
      }
    }
  }

  # Create plots
  pdf_file <- file.path(output_dir, paste0(output_prefix, "_enrichment.pdf"))
  pdf(pdf_file, width = ENRICHMENT_PLOT_WIDTH, height = ENRICHMENT_PLOT_HEIGHT)

  par(mfrow=c(3,1), mar=c(2,22,2,2))

  for(i in seq_along(ENRICHR_DBS)) {
    db_name <- ENRICHR_DBS[i]
    if (!"Adjusted.P.value" %in% colnames(enriched[[i]])) next
    significant_results <- enriched[[i]] %>% filter(Adjusted.P.value < P_VALUE_CUTOFF)
    
    if(nrow(significant_results) > 0) {
      n_terms <- min(TOP_N_ENRICHMENT, nrow(significant_results))
      tmp_results <- significant_results[1:n_terms, c("Term", "Combined.Score")]
      tmp_results <- tmp_results[order(tmp_results$Combined.Score, decreasing = FALSE),]
      
      barplot(tmp_results$Combined.Score, horiz = TRUE, col = "blue",
              names.arg = wrapText(tmp_results$Term, 45), 
              cex.names = 1.0, las = 1,
              main = paste(title_prefix, "-", db_name),
              xlab = "Combined Score")
      abline(v = 2, col = "red")
    } else {
      plot(1, type = "n", xlab = "", ylab = "", axes = FALSE)
      text(1, 1, paste("No significant terms in", db_name), cex = 1.2)
    }
    
    if(i %% 3 == 0 && i < length(ENRICHR_DBS)) {
      par(mfrow=c(3,1), mar=c(2,22,2,2))
    }
  }
  
  dev.off()
  return(enriched)
}

# Run enrichment
if (length(up_genes) > 0) {
  run_enrichr_analysis(up_genes, paste0(base_name, "_UP"), "Upregulated")
}

if (length(down_genes) > 0) {
  run_enrichr_analysis(down_genes, paste0(base_name, "_DOWN"), "Downregulated")
}

cat(" GO enrichment complete\n")
