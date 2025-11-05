###############################################################################
# Graph Difference Operator
# -----------------------------------------------------------------------------
# This script constructs graph difference operators D for common graphs:
#   (1) create.chain(p):       chain graph with p nodes and p-1 edges.
#   (2) create.lattice(s,t):   s×t lattice graph.
#   (3) create.complete(p):    complete graph with p nodes and p(p−1)/2 edges.
#   (4) create.general(edges): general undirected graph from edge list.
#
# Notes:
#   - Rows correspond to edges; columns correspond to nodes.
#   - Each row has exactly two nonzeros: (+1 at tail, −1 at head).
#   - For chain/lattice, the orientation follows increasing index order.
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
### (1) Chain Graph
### Note: same as getDtf(p,0)/get1D(p) function in genlasso package
create.chain <- function(p){
  assoc_matrix <- matrix(0, nrow = p-1, ncol = p)
  for(i in 1:(p-1)){
    assoc_matrix[i,i] <- -1
    assoc_matrix[i, i+1] <- 1
  }
  return(assoc_matrix)
}
###############################################################################
###############################################################################
###############################################################################
### (2) Lattice Graph (2D Grid)
### Note: same as getD2d(t,s) function in genlasso package (order changed)
create.lattice <- function(s, t) {
  num_edges <- (s - 1) * t + (t - 1) * s
  num_nodes <- s * t
  assoc_matrix <- matrix(0, nrow = num_edges, ncol = num_nodes)
  edge_counter <- 1
  # Horizontal edges
  for (row in 1:s) {
    for (col in 1:(t-1)) {
      node <- (row - 1) * t + col
      assoc_matrix[edge_counter, node] <- -1
      assoc_matrix[edge_counter, node + 1] <- 1
      edge_counter <- edge_counter + 1
    }
  }
  # Vertical edges
  for (col in 1:t) {
    for (row in 1:(s-1)) {
      node <- (row - 1) * t + col
      assoc_matrix[edge_counter, node] <- -1
      assoc_matrix[edge_counter, node + t] <- 1
      edge_counter <- edge_counter + 1
    }
  }
  return(assoc_matrix)
}
################################################################################
################################################################################
################################################################################
### (3) Complete Graph
create.complete <- function(p) {
  num_edges <- p * (p - 1) /2
  assoc_matrix <- matrix(0, nrow = num_edges, ncol = p)
  edge_counter <- 1
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      assoc_matrix[edge_counter, i] <- 1
      assoc_matrix[edge_counter, j] <- -1
      edge_counter <- edge_counter + 1
    }
  }
  return(assoc_matrix)
}
################################################################################
################################################################################
################################################################################
### P.S. General Graph
# example: edges <- list(c(1, 2), c(2, 3), c(1, 3))
create.general <- function(edges) {
  num_nodes <- max(unlist(edges))
  num_edges <- length(edges)
  assoc_matrix <- matrix(0, nrow = num_edges, ncol = num_nodes)
  for (edge_index in 1:num_edges) {
    u <- edges[[edge_index]][1]
    v <- edges[[edge_index]][2]
    assoc_matrix[edge_index, u] <- 1  
    assoc_matrix[edge_index, v] <- -1   
  }
  return(assoc_matrix)
}