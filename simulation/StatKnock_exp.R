###############################################################################
# Simulation Experiments for StatKnock
# -----------------------------------------------------------------------------
# This script reproduces the results in Figure S.4 and Figure S.5 in the Supplement
################################################################################
################################################################################
################################################################################
#install.packages("GKnockoff_0.1.0.tar.gz", repos = NULL, type = "source")
#install.packages("SplitKnockoff_1.2.tar.gz", repos = NULL, type = "source")
source("helper/Graph Difference Operator.R")
source("helper/Beta Structure Generation.R")
source("../function/StatKnock.R")
source("../function/other methods/BY_Structural.R")
source("../function/other methods/DS_Structural.R")
library(GKnockoff)     # package version 0.1.0
library(SplitKnockoff) # package version 1.2
library(genlasso)
library(hdi)
library(Rfast)
################################################################################
################################################################################
################################################################################
### Helper function
fdp_power <- function(selected_index, signal_index){
  fdp=0
  power=0
  ssr=0
  if (length(selected_index)!=0){
    num_selected <- length(selected_index)
    tp <- length(intersect(selected_index, signal_index))
    fp <- num_selected - tp
    fdp <- fp / max(num_selected, 1)
    power <- tp / length(signal_index)
    if (power==1){
      ssr=1
    }else{
      ssr=0
    }
  }  
  return(list(fdp = fdp, power = power, ssr=ssr))
}

mean.se = function(x, ns = 1) { 
  paste0(format(round(mean(x)*100, ns), nsmall=ns), ' (', 
         format(round(100*sd(x)/sqrt(length(x)), ns), nsmall=ns), ')')}
################################################################################
################################################################################
################################################################################
set.seed(2025)
nIterations = 100
q = 0.2
################################################################################
################################################################################
################################################################################
### Figure S.4 (a): piecewise constant over chain (n=500, p=120)
set.seed(2025)
p = 120; n = 500
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- c(0.2, 0.5, 0.8)
ChainPiecewiseConstant_result_all_n500 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      y <- X %*% Beta0 + rnorm(n, sd = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-10)
      
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ### DS & MDS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel  <- DS_Fit$DS_feature
      MDS_Sel <- DS_Fit$MDS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      fdp_power5 = fdp_power(selected_index = MDS_Sel, signal_index = signal_ind)
      
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ### Split Knockoff
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = 10.^seq(0, 2, length.out = 10),
        sign = 'disabled'
      )
      SplitKnockoff_plus_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_plus_Sel <- SplitKnockoff_plus_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_plus_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power5$fdp, fdp_power5$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_result_all_n500[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_result_all_n500[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_result_all_n500, file = "StatKnock_ChainPiecewiseConstant_n500signal12.RData")
################################################################################
################################################################################
################################################################################
### Figure S.4 (b): piecewise constant over chain (n=250, p=120)
set.seed(2025)
p = 120; n = 250
p1 = 120
m1 = 12
D <- create.chain(p)
m <- nrow(D)

A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
rho_values <- c(0.2, 0.5, 0.8)
ChainPiecewiseConstant_result_all_n250 <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))
  
  ChainPiecewiseConstant_list <- list()
  
  for (A in A_values) {
    cat("  Running for A =", A, "\n")
    
    fdp_power_list_ChainPiecewiseConstant = sapply(1:nIterations, function(it) {
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      y <- X %*% Beta0 + rnorm(n, sd = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-10)
      
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      fdp_power1
      fdp_power2
      
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      fdp_power3
      
      ### DS & MDS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel  <- DS_Fit$DS_feature
      MDS_Sel <- DS_Fit$MDS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      fdp_power5 = fdp_power(selected_index = MDS_Sel, signal_index = signal_ind)
      
      ### GKnockoff
      GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      GKnockoff_Sel <- GKnockoff_Fit$selected
      fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ### Split Knockoff
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = 10.^seq(0, 2, length.out = 10),
        sign = 'disabled'
      )
      SplitKnockoff_plus_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_plus_Sel <- SplitKnockoff_plus_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_plus_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power5$fdp, fdp_power5$power,
            fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    ChainPiecewiseConstant_list[[paste0("A=", A)]] <- fdp_power_list_ChainPiecewiseConstant
  }
  
  ChainPiecewiseConstant_result_all_n250[[paste0("rho=", rho)]] <- ChainPiecewiseConstant_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- ChainPiecewiseConstant_result_all_n250[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(ChainPiecewiseConstant_result_all_n250, file = "StatKnock_ChainPiecewiseConstant_n250signal12.RData")
################################################################################
################################################################################
################################################################################
### Figure S.5: 1D Sparse fused Lasso

set.seed(2025)
p = 150; n = 600
p1 = 40
m1 = 12
D <- rbind(diag(p), create.chain(p))
m <- nrow(D)

rho_values <- c(0.5, 0.8)
A_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)

SPL_1D_result_all <- list()

for (rho in rho_values) {
  cat("Running for rho =", rho, "\n")
  
  Sigma = toeplitz(rho^(1:p - 1))

  SPL_1D_list <- list()
  
  for (A in A_values) {
    cat("Running for A =", A, "\n")
    
    fdp_power_list_SPL_1D = sapply(1:nIterations, function(it) {
      ### Generate the Data
      Beta0 = Chain_PiecewiseConstant(p, p1, m1, A, min_interval = 3, edge_buffer = 3)
      X <- matrix(rnorm(n * p), n) %*% chol(Sigma)
      y <- X %*% Beta0 + rnorm(n, sd = 1)
      signal_ind = which(abs(D %*% Beta0) > 1e-10)
      
      ### StatKnock
      StatKnock_Fit <- StatKnock(X = X, y = y, D = D, q = q)
      StatKnock_Sel <- StatKnock_Fit$S
      StatKnock_plus_Sel <- StatKnock_Fit$S_plus
      fdp_power1 = fdp_power(selected_index = StatKnock_Sel, signal_index = signal_ind)
      fdp_power2 = fdp_power(selected_index = StatKnock_plus_Sel, signal_index = signal_ind)
      
      ### B-Y
      BY_Sel = BY_GLasso(X = X, y = y, D = D, q = q)
      fdp_power3 = fdp_power(selected_index = BY_Sel, signal_index = signal_ind)
      
      ### DS & MDS
      DS_Fit = DS_D(X, y, D = D, num_split = 50, q = q)
      DS_Sel  <- DS_Fit$DS_feature
      MDS_Sel <- DS_Fit$MDS_feature
      fdp_power4 = fdp_power(selected_index = DS_Sel, signal_index = signal_ind)
      fdp_power5 = fdp_power(selected_index = MDS_Sel, signal_index = signal_ind)
      
      # ### GKnockoff
      # GKnockoff_Fit = GKnockoff(X, y, D = D, fdr = q)
      # GKnockoff_Sel <- GKnockoff_Fit$selected
      # fdp_power6 = fdp_power(selected_index = GKnockoff_Sel, signal_index = signal_ind)
      
      ### Split Knockoff
      option <- list(
        q = q,
        method = 'knockoff+',
        eta = 0.1,
        normalize = 'true',
        copy = 'true',
        lambda = 10.^seq(0, -6, by = -0.01),
        nu = 10.^seq(0, 2, length.out = 10),
        sign = 'disabled'
      )
      SplitKnockoff_plus_Fit = sk.filter(X, D, y, option)
      SplitKnockoff_plus_Sel <- SplitKnockoff_plus_Fit$results$data
      fdp_power7 = fdp_power(selected_index = SplitKnockoff_plus_Sel, signal_index = signal_ind)
      
      r = c(fdp_power1$fdp, fdp_power1$power,
            fdp_power2$fdp, fdp_power2$power,
            fdp_power3$fdp, fdp_power3$power,
            fdp_power4$fdp, fdp_power4$power,
            fdp_power5$fdp, fdp_power5$power,
            #fdp_power6$fdp, fdp_power6$power,
            fdp_power7$fdp, fdp_power7$power
      )
      return(r)
    })
    
    SPL_1D_list[[paste0("A=", A)]] <- fdp_power_list_SPL_1D
  }
  
  SPL_1D_result_all[[paste0("rho=", rho)]] <- SPL_1D_list
}

for (rho in rho_values) {
  cat("\nResults for rho =", rho, "\n")
  for (A in A_values) {
    result_matrix <- SPL_1D_result_all[[paste0("rho=", rho)]][[paste0("A=", A)]]
    cat("A =", A, "\n")
    print(apply(result_matrix, 1, mean.se))
  }
}

save(SPL_1D_result_all, file = "StatKnock_SPL_1D.RData")