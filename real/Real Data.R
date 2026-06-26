###############################################################################
# Real Data Experiment 
# -----------------------------------------------------------------------------
# This script applies multiple structural FDR-control methods to brain imaging
# data from the ADNI (Alzheimer's Disease Neuroimaging Initiative) dataset.
# The data include subject-level sMRI features (Level 5 parcellation).
# 
# Note: Due to privacy and data-sharing restrictions, the ADNI dataset is not
# publicly distributed. One may follow the steps in README.md to obtain
# and preprocess the data.
#
# Workflow:
#   1) Load sMRI feature matrix (X) and response (y).
#   2) Construct transformational matrices D1, D2, D3:
#        - D1: identity (zero-order)
#        - D2: first-order transformation
#        - D3: second-order transformation (trimmed to remove zero rows)
#   3) Apply multiple methods:
#        - SATLasso + StatKnock
#        - BY for structural sparsity
#        - High-dimensional SplitKnockoff
#   4) Count feature/edge selection frequencies across 100 iterations
#        for each method–structure combination.
#
# Outputs:
#   - frequency_Level5.csv: combined summary of selection frequencies
#     (columns: Method, Group [1–3], Index, Count)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
source("../function/SATLasso+StatKnock.R")
source("../function/other methods/HD_SplitKnockoff.R")
source("../function/other methods/BY_Structural.R")
source("../simulation/helper/Graph Difference Operator.R")
library(GKnockoff)     # package version 0.1.0
library(SplitKnockoff) # package version 1.2
library(genlasso)
library(hdi)
library(Rfast)
library(dplyr)
library(readr)
library(readxl)
################################################################################
################################################################################
################################################################################
set.seed(2025)
nIterations = 100
q = 0.2
################################################################################
################################################################################
################################################################################
## Load Data
data <- read_csv("ADNI_sMRI_Level5.csv") 
y <- as.matrix(data[, 4])
X <- as.matrix(data[, 5:ncol(data)])

n = dim(X)[1]
p = dim(X)[2]
######################################
######################################
######################################
## Create the Graph
data_lookup <- read.csv("lookup/multilevel_lookup_table.csv")
pairs <- list()
# Loop over each possible pair
for (i in 1:(nrow(data_lookup) - 1)) {
  for (j in (i + 1):nrow(data_lookup)) {
    # Check if levels match
    if (data_lookup$Level_3[i] == data_lookup$Level_3[j] &&
        data_lookup$Level_2[i] == data_lookup$Level_2[j] && data_lookup$Level_1[i] == data_lookup$Level_1[j]) {
      # Store indices of matching pairs
      pairs[[length(pairs) + 1]] <- c(i, j)
    }
  }
}
pair_indices <- do.call(rbind, pairs)
pair_indices <- as.matrix(pair_indices)

col1 <- data_lookup[[1]]
col2 <- data_lookup[[2]]
src_endpoint <- paste(col1[pair_indices[,1]], col2[pair_indices[,1]], sep = "_")
dst_endpoint <- paste(col1[pair_indices[,2]], col2[pair_indices[,2]], sep = "_")
edge_df <- data.frame(
  source = src_endpoint,
  target = dst_endpoint,
  stringsAsFactors = FALSE
)
#write.csv(edge_df, "edge_endpoints.csv", row.names = T)
edge_df <- read.csv("lookup/edge_endpoints.csv")
######################################
######################################
######################################
## Different Transformation Matricies

D_1 <- diag(p)

D_2 <- matrix(0, nrow = nrow(pair_indices), ncol = nrow(data_lookup))
for (k in 1:nrow(pair_indices)) {
  D_2[k, pair_indices[k, 1]] <- -1
  D_2[k, pair_indices[k, 2]] <- 1
}

D_3_raw <- t(D_2) %*% D_2

process_matrix <- function(D) {
  keep_rows <- which(!apply(D, 1, function(row) all(row == 0)))
  
  D_trimmed <- D[keep_rows, , drop = FALSE]
  
  list(
    matrix = D_trimmed,
    rows_kept = keep_rows
  )
}

D_3 <- process_matrix(D_3_raw)$matrix
D_3_index = process_matrix(D_3_raw)$rows_kept
######################################
######################################
######################################
m_1 <- nrow(D_1)
print(m_1)
m_2 <- nrow(D_2)
print(m_2)
m_3 <- nrow(D_3)
print(m_3)

methods <- c("SATLasso_StatKnock", "BY", "SplitKnockoff")

for (i in 1:3) {
  m_val <- get(paste0("m_", i))
  for (method in methods) {
    assign(paste0(method, "_frequency_", i), numeric(m_val))
  }
}

for (i in 1:nIterations) {
  print(i)
  ### SATLasso+StatKnock
  SATLasso_StatKnock_Fit_1 <- SATLasso_StatKnock(X = X, y = y, D = D_1, q = q)
  SATLasso_StatKnock_Sel_1 <- SATLasso_StatKnock_Fit_1$S
  SATLasso_StatKnock_Fit_2 <- SATLasso_StatKnock(X = X, y = y, D = D_2, q = q)
  SATLasso_StatKnock_Sel_2 <- SATLasso_StatKnock_Fit_2$S
  SATLasso_StatKnock_Fit_3 <- SATLasso_StatKnock(X = X, y = y, D = D_3, q = q)
  SATLasso_StatKnock_Sel_3 <- D_3_index[SATLasso_StatKnock_Fit_3$S]
  SATLasso_StatKnock_Sel_1
  SATLasso_StatKnock_Sel_2
  SATLasso_StatKnock_Sel_3
  
  ### B-Y
  BY_Sel_1 = BY_GLasso(X = X, y = y, D = D_1, q = q)
  BY_Sel_2 = BY_GLasso(X = X, y = y, D = D_2, q = q)
  BY_Sel_3 = D_3_index[BY_GLasso(X = X, y = y, D = D_3, q = q)]
  BY_Sel_1
  BY_Sel_2
  BY_Sel_3
  
  ### Split Knockoff
  SplitKnockoff_Fit_1 = HD_SplitKnockoff(X, y, D = D_1, q = q)
  SplitKnockoff_Sel_1 <- SplitKnockoff_Fit_1$Sel
  SplitKnockoff_Fit_2 = HD_SplitKnockoff(X, y, D = D_2, q = q)
  SplitKnockoff_Sel_2 <- SplitKnockoff_Fit_2$Sel
  SplitKnockoff_Fit_3 = HD_SplitKnockoff(X, y, D = D_3, q = q)
  SplitKnockoff_Sel_3 <- D_3_index[SplitKnockoff_Fit_3$Sel]
  SplitKnockoff_Sel_1
  SplitKnockoff_Sel_2
  SplitKnockoff_Sel_3
  
  for (method in methods) {
    for (g in 1:3) {
      freq_name <- paste0(method, "_frequency_", g)
      sel_name <- paste0(method, "_Sel_", g)
      
      if (exists(freq_name) && exists(sel_name)) {
        freq <- get(freq_name)
        sel <- get(sel_name)
        
        freq[sel] <- freq[sel] + 1
        
        assign(freq_name, freq)
      }
    }
  }
}

summary_list <- list()

for (method in methods) {
  for (i in 1:3) {
    freq_name <- paste0(method, "_frequency_", i)
    if (exists(freq_name)) {
      freq <- get(freq_name)
      nonzero_index <- which(freq > 0)
      nonzero_freq <- freq[nonzero_index]
      
      if (length(nonzero_index) > 0) {
        sorted_order <- order(nonzero_freq, decreasing = TRUE)
        sorted_index <- nonzero_index[sorted_order]
        sorted_freq <- nonzero_freq[sorted_order]
        
        df <- data.frame(
          Method = method,
          Group = i,
          Index = sorted_index,
          Count = sorted_freq
        )
        
        summary_list[[paste0(method, "_", i)]] <- df
      }
    }
  }
}

final_summary <- do.call(rbind, summary_list)

write.csv(final_summary, file = "frequency_Level5.csv", row.names = FALSE)