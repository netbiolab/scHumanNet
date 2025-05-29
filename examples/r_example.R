#!/usr/bin/env Rscript
#' Example usage of scHumanNet R package
#'
#' This script demonstrates how to use the scHumanNet R package for
#' cell-type-specific network analysis.

library(scHumanNet)  # Assuming package is installed
library(igraph)
library(dplyr)

# Function to create example data
create_example_data <- function() {
  set.seed(42)
  
  # Generate example gene names
  genes <- paste0("GENE", sprintf("%04d", 1:1000))
  
  # Create networks for different cell types and conditions
  networks <- list()
  
  for (celltype in c('Neuron', 'Astrocyte', 'Microglia')) {
    for (condition in c('Control', 'Disease')) {
      # Generate random network
      n_edges <- sample(200:500, 1)
      
      network_data <- data.frame(
        gene1 = sample(genes, n_edges, replace = TRUE),
        gene2 = sample(genes, n_edges, replace = TRUE),
        LLS = runif(n_edges, 0.1, 1.0),
        scinet_weight = runif(n_edges, 0.1, 1.0),
        stringsAsFactors = FALSE
      )
      
      # Remove self-loops
      network_data <- network_data[network_data$gene1 != network_data$gene2, ]
      
      network_name <- paste(condition, celltype, sep = "_")
      networks[[network_name]] <- network_data
    }
  }
  
  # Create example metadata
  n_cells <- 3000
  metadata <- data.frame(
    cell_id = paste0("Cell_", sprintf("%04d", 1:n_cells)),
    celltype = sample(c('Neuron', 'Astrocyte', 'Microglia'), n_cells, replace = TRUE),
    condition = sample(c('Control', 'Disease'), n_cells, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  return(list(networks = networks, metadata = metadata))
}

# Main analysis function
main <- function() {
  cat("scHumanNet R Example\n")
  cat(paste(rep("=", 40), collapse = ""), "\n")
  
  # Create example data
  cat("Creating example data...\n")
  example_data <- create_example_data()
  networks <- example_data$networks
  metadata <- example_data$metadata
  
  cat("Created", length(networks), "networks\n")
  cat("Created metadata for", nrow(metadata), "cells\n")
  
  # Sort networks and add LLS weights
  cat("\nSorting networks and adding LLS weights...\n")
  sorted_networks <- SortAddLLS(networks)
  
  # Calculate network statistics
  cat("\nCalculating network statistics...\n")
  for (name in names(sorted_networks)) {
    net_df <- sorted_networks[[name]]
    cat(sprintf("%s: %d edges\n", name, nrow(net_df)))
  }
  
  # Calculate centralities
  cat("\nCalculating degree centrality...\n")
  centralities <- GetCentrality('degree', sorted_networks)
  
  # Combine percentile ranks
  cat("Combining percentile ranks...\n")
  rank_df <- CombinePercRank(centralities)
  cat(sprintf("Combined centrality matrix: %d x %d\n", nrow(rank_df), ncol(rank_df)))
  
  # Find hub genes
  cat("\nFinding significant hub genes...\n")
  tryCatch({
    hub_results <- FindAllHub(
      net.list = sorted_networks,
      centrality = 'degree',
      threshold = 0.05
    )
    
    cat("Found", nrow(hub_results), "significant hub genes\n")
    if (nrow(hub_results) > 0) {
      cat("Top 5 hub genes:\n")
      top_hubs <- head(hub_results[order(hub_results$Centrality_PR, decreasing = TRUE), ], 5)
      for (i in 1:min(5, nrow(top_hubs))) {
        row <- top_hubs[i, ]
        cat(sprintf("  %s (%s): PR=%.3f, q=%.3e\n", 
                   row$gene, row$celltype, row$Centrality_PR, row$qvalue))
      }
    }
  }, error = function(e) {
    cat("Hub analysis failed:", e$message, "\n")
  })
  
  # Differential analysis
  cat("\nRunning differential network analysis...\n")
  tryCatch({
    diff_results <- FindDiffHub(
      rank.df.final = rank_df,
      meta = metadata,
      celltypes = 'celltype',
      condition = 'condition',
      control = 'Control',
      net.list = sorted_networks,
      min.cells = 100  # Lower threshold for example
    )
    
    cat("Found", nrow(diff_results), "differential results\n")
    if (nrow(diff_results) > 0) {
      # Show top differential genes
      diff_results$abs_diffPR <- abs(diff_results$diffPR)
      top_diff <- head(diff_results[order(diff_results$abs_diffPR, decreasing = TRUE), ], 5)
      cat("Top 5 differential genes:\n")
      for (i in 1:min(5, nrow(top_diff))) {
        row <- top_diff[i, ]
        cat(sprintf("  %s (%s): diffPR=%.3f, q=%.3e\n", 
                   row$gene, row$celltype, row$diffPR, row$qvalue))
      }
    }
  }, error = function(e) {
    cat("Differential analysis failed:", e$message, "\n")
  })
  
  cat("\nExample completed successfully!\n")
}

# Run the example
if (!interactive()) {
  main()
}