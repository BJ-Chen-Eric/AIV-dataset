#!/usr/bin/env Rscript

# Load necessary libraries and suppress startup messages
suppressPackageStartupMessages({
  library(castor)
  library(argparse)
  library(stringr)
  library(dplyr)
  library(tibble)
})

# Source external functions (ensure this file contains necessary helper functions)
source('~/R/aiv/function.R')

# Function to parse command-line arguments
parse_args_custom <- function(test_args = NULL) {
  parser <- ArgumentParser(description = "Compute patristic distance matrix and assign groups based on MPD thresholds.")
  
  # Define all your arguments
  parser$add_argument("-seg", "--segment", type = 'character', required = TRUE,
                      help = "Segment to compute the distance matrix for.")
  
  parser$add_argument("-tree", "--tree_dir", type = 'character', required = TRUE,
                      help = "Directory containing tree files.")
  
  parser$add_argument("-p", "--prefix", type = 'character', required = TRUE,
                      help = "Prefix of the tree files (e.g., '_fasttree.nexus').")
  
  parser$add_argument("-pd", "--pd_dir", type = 'character', required = FALSE,
                      help = "Path to precomputed patristic distance RData file. If not provided, it will be computed from the tree.")
  
  parser$add_argument("-of", "--out_fix", type = 'character', required = TRUE,
                      help = "Output file prefix.")
  
  parser$add_argument("-o", "--out_dir", type = 'character', required = TRUE,
                      help = "Directory to store the generated distance matrix and group assignments.")
  
  if (is.null(test_args)) {
    args <- parser$parse_args()
  } else {
    args <- parser$parse_args(test_args)
  }
  
  return(args)
}

# Function to read and clean the tree
read_and_clean_tree <- function(tree_dir, segment, prefix) {
  # Construct the file pattern to match
  tree_path_pattern <- paste0(segment, prefix, "$")
  
  # Read the tree using read_as_list (assumed to handle multiple trees if present)
  tree_list <- read_as_list(path = tree_dir, file_type = 'nexus', prefix = tree_path_pattern)
  
  if (length(tree_list) == 0) {
    stop("No tree files matched the specified pattern.")
  }
  
  # Assuming only one tree is matched; adjust if multiple trees are expected
  tree <- tree_list[[1]]
  
  # Clean tip labels by removing single quotes
  tree$tip.label <- gsub("'", "", tree$tip.label)
  
  return(tree)
}

# Corrected Function to compute patristic distances
compute_patristic_distances <- function(tree) {
  cat("Computing all pairwise patristic distances...\n")
  Ntips <- length(tree$tip.label)
  distances <- get_all_pairwise_distances(tree)[1:Ntips, 1:Ntips]
  colnames(distances) <- rownames(distances) <- tree$tip.label
  
  # Extract distances for tip labels only
  tip_labels <- tree$tip.label
  distances_tips <- distances[tip_labels, tip_labels]
  
  return(distances_tips)
}

# Function to extract leaf nodes for each internal node
extract_leaf_nodes <- function(tree) {
  cat("Extracting leaf nodes for each internal node...\n")
  Ntips <- length(tree$tip.label)
  Nnodes <- tree$Nnode
  
  node_leafe <- vector("list", Nnodes)  # Pre-allocate list for performance
  
  for (i in 1:Nnodes) {
    n_node <- Ntips + i
    n_desc <- getDescendants(tree, n_node)
    
    # Filter to include only tips
    leaf_indices <- n_desc[n_desc <= Ntips]
    
    # Extract and clean tip labels
    node_leafe[[i]] <- tree$tip.label[leaf_indices] %>% 
      str_remove_all(pattern = "'")
  }
  
  names(node_leafe) <- as.character(1:Nnodes)
  
  return(node_leafe)
}

# Function to save the patristic distances and leaf nodes
save_patristic_results <- function(distances, node_leafe, out_dir, segment, out_fix) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    cat(paste("Created output directory:", out_dir, "\n"))
  }
  
  output_path <- file.path(out_dir, paste0(segment, out_fix, "_patristic_dist_matrix.RData"))
  saveRDS(list(d = distances, leafe = node_leafe), file = output_path)
  cat(paste("Saved patristic distances and leaf nodes to:", output_path, "\n"))
}

# Function to extract strain number from header
strain_number_extract <- function(header) {
  str_extract(header, '[0-9]+_H|[0-9]+\\|H') %>% str_remove('_H|\\|H')
}

# Corrected Function to perform MPD threshold filtering and group assignment
perform_threshold_filter <- function(j, mean_all, leafe_df, node_leafe) {
  mean_filtered <- mean_all %>% filter(MPD < j)
  node_list <- c()
  tips <- leafe_df$epi %>% unique() %>% length
  leafe <- leafe_df[leafe_df$node %in% mean_filtered$index, ]
  
  node_list <- c()
  x <- 0
  while(nrow(mean_filtered) >0) {
    x <- x+1
    # mean_filtered[1, 1]
    seed <- leafe[leafe$node %in% mean_filtered[1, 1], 2] %>% unique() #
    seed_node <- leafe[leafe$epi %in% seed, 1] %>% unique()
    seed_node_epi <- leafe[leafe$node %in% seed_node, 2] %>% unique()
    node <- leafe[leafe$epi %in% seed_node_epi, 1] %>% table_DF() %>% arrange(desc(Freq))
    # print(node)
    # node <- leafe_all[leafe_all$epi %in% seed_node_epi, 1] %>% table_DF() %>% arrange(desc(Freq))
    print(node[1, ])
    mean_filtered <- mean_filtered[!mean_filtered$index %in% (node$x %>% unlist), ]
    leafe <- leafe[!leafe$node %in% (node$x %>% unlist), ]
    
    node_list[x] <- node$x[1]
    node_list
  }
  
  
  A <- node_leafe[node_list]
  A <- A[which((lapply(A, length) %>% unlist) > 2)] # remove those unclassified
  for(k in names(A)) {
    A[[k]] <- data.frame(clade=k, leafe=A[[k]])
  }
  A <- A %>% do.call(what=rbind) %>% 
    mutate(Strain_number=str_extract(leafe, pattern = '[0-9]+_H|[0-9]+\\|H') %>% str_remove(pattern = '_H|\\|H')) 
  colnames(A)[1] <- paste('d_', j, sep = '')
  return(A)
  # colnames(out) <- j
  # out[, 2] %>% table()
  # if(((out[, y] %>% unique() %>% length()))!=1) {
  #   if(j == 0.01) {out <- out[, 1] %>% as.data.frame() %>% rename('Strain_number'='.')}
  #  
  # }else(break)
  
}

# Function to save MPD and group assignments
save_mpd_results <- function(node_MPD, group, out_dir, segment, out_fix) {
  dir.create(out_dir)
  output_path <- file.path(out_dir, paste0(segment, out_fix, "_MPD_groups.RData"))
  saveRDS(list(node_MPD = node_MPD, group = group), file = output_path)
  cat(paste("Saved MPD and group assignments to:", output_path, "\n"))
}




# Main execution ----------------------------------------------------------
# Main execution function
main <- function(args = NULL) {
  if (is.null(args)) {
    args <- parse_args_custom()
  }
  
  cat(paste("Starting patristic distance computation for segment:", args$segment, "\n"))
  
  
  # Read and clean the tree
  tree <- read_and_clean_tree(args$tree_dir, args$segment, args$prefix)
  
  # Compute or load patristic distances
  if(!is.null(args$pd_dir) && file.exists(args$pd_dir)) {
    # Load precomputed patristic distances
    cat("Loading precomputed patristic distances...\n")
    pd_list <- readRDS(args$pd_dir)
    distances <- pd_list$d
    node_leafe <- pd_list$leafe
    
  } else {
    # Compute patristic distances using the corrected function
    distances <- compute_patristic_distances(tree)
    
    # Extract leaf nodes for each internal node
    node_leafe <- extract_leaf_nodes(tree)
    
    # Save patristic distances and leaf nodes
    save_patristic_results(distances, node_leafe, args$out_dir, args$segment, args$out_fix)
  }
  
  # Compute node_MPD
  cat("Calculating median pairwise distances (MPD) for each node...\n")
  node_MPD <- lapply(seq_along(node_leafe), function(j) {
    A <- node_leafe[[j]]
    n <- which(colnames(distances) %in% A)
    if(length(n) == 0) {
      warning(paste("No matching columns found for node", j))
      return(NA)
    }
    median(distances[n, n], na.rm = TRUE)
  })
  
  # Clean up memory
  rm(distances)
  gc()
  
  # Create a dataframe for MPD and leafe
  mean_all <- data.frame(index = 1:length(node_MPD), 
                         MPD = unlist(node_MPD)) %>%
    mutate(leafe = sapply(node_leafe, length)) %>%
    arrange(leafe)
  
  # Extract strain number and prepare output dataframe
  # Note: Assuming node_leafe[[1]] contains strain headers
  out <- data.frame(header = node_leafe[[1]]) %>%
    mutate(Strain_number = strain_number_extract(header)) %>%
    select(Strain_number)
  
  # Prepare 'leafe' dataframe for threshold filtering
  leafe_df <- node_leafe %>%
    lapply(FUN = function(x) { data.frame(epi = x) }) %>%
    do.call(what=rbind) %>% as.data.frame() %>% 
    rownames_to_column('node') %>%
    mutate(node = as.numeric(str_remove(node, '\\.[0-9]+')))  # Ensure 'node' is numeric
  
  # Verify 'leafe_df' structure
  if(!all(c("node", "epi") %in% colnames(leafe_df))) {
    stop("leafe_df does not contain the required 'node' and 'epi' columns.")
  }
  
  # Define the sequence of thresholds
  thresholds <- seq(0.01, max(mean_all$MPD), 0.01) %>% round(digits = 2)
  y <- 0
  
  for (j in thresholds) {
    gc()
    cat(paste('Processing threshold', j, '\n'))
    y <- y + 1
    
    A <- perform_threshold_filter(j, mean_all, leafe_df, node_leafe)
    
    if(is.null(A)) {
      cat(paste("No groups found for threshold", j, ". Stopping.\n"))
      break
    }
    out <- merge(out, A[, c(1, 3)], by = 'Strain_number', all = T)
    # out <- merge(out, A, by = 'Strain_number', all = TRUE)
  }
  
  # Ensure all NAs are filled appropriately
  rank <- tree$tip.label %>% strain_number_extract()
  rownames(out) <- out$Strain_number
  out <- out[rank, ]
  rownames(out) <- 1:nrow(out)
  
  
  for(j in 2:ncol(out)) {
    for(k in is.na(out[, j]) %>% which()) {
      if(k==1) {out[k, j] <- out[k+1, j]}
      else(out[k, j] <- out[k-1, j])
    }
  }
  
  # Save MPD and group assignments
  save_mpd_results(node_MPD, out, args$out_dir, args$segment, args$out_fix)
  
  cat("MPD calculation and group assignment completed successfully.\n")
}

# Define a manual args list for interactive testing
# Uncomment the following lines to test interactively
# gisaid_IRD_merged_N2_aligned_iqtree.nexus
# args <- list(
#   segment = "N9",
#   tree_dir = "/home/eric/Analysis/aiv/merge/0307/iq_tree/",
#   # prefix = "_aligned_fasttree.tree",
#   prefix = "_aligned_iqtree.nexus",
#   pd_dir = '~/Analysis/aiv/merge/0307/distance/iq/N9_iq_patristic_dist_matrix.RData',  # Set to NULL to compute from tree
#   out_fix = "_iq",
#   out_dir = "/home/eric/Analysis/aiv/merge/0307/distance/NA_1009/"
# )

# Run the main function with test_args (for interactive testing)
# Uncomment the following line when testing interactively
# main(args = test_args)

# For command-line execution, keep the main call without arguments
main()
