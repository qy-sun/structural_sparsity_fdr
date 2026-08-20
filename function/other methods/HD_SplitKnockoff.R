###############################################################################
# High-dimensional SplitKnockoff and SplitLasso
# -----------------------------------------------------------------------------
# This script implements the high-dimensional SplitKnockoff and SplitLasso.
#
# Dependencies: SplitKnockoff (version 1.2)
#
# Reference:
#   Cao, Y., Sun, X., & Yao, Y. (2024).
#   Controlling the false discovery rate in transformational sparsity: Split knockoffs.
#   Journal of the Royal Statistical Society: Series B (Statistical Methodology),
#   86(2), 386–410.
###############################################################################
library(glmnet)
#install.packages("SplitKnockoff_1.2.tar.gz", repos = NULL, type = "source")
library(SplitKnockoff)
###############################################################################
###############################################################################
###############################################################################
## (1) SplitLasso
SplitLasso <- function(X, y, D, nu_values = 10^seq(-1, 3, length.out = 10)) {
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- nrow(D)
  
  objective_function <- function(nu) {
    creat.result  <- sk.create(X, y, D, nu)
    A_beta  <- creat.result$A_beta
    A_gamma <- creat.result$A_gamma
    tilde_y <- creat.result$tilde_y
    
    penalty <- matrix(1, p + m, 1)
    penalty[1:p, 1] <- 0
    
    cv <- cv.glmnet(cbind(A_beta, A_gamma), tilde_y, penalty.factor = penalty, alpha = 1, nfolds = 5)
    return(min(cv$cvm))
  }
  
  tryCatch({
    cv_errors <- sapply(nu_values, objective_function)
    
    best_nu <- nu_values[which.min(cv_errors)]
    print(paste("best nu =", best_nu))
    
    creat.result  <- sk.create(X, y, D, best_nu)
    A_beta  <- creat.result$A_beta
    A_gamma <- creat.result$A_gamma
    tilde_y <- creat.result$tilde_y
    
    penalty <- matrix(1, p + m, 1)
    penalty[1:p, 1] <- 0
    
    cv <- cv.glmnet(cbind(A_beta, A_gamma), tilde_y, penalty.factor = penalty, alpha = 1)
    best_lambda <- cv$lambda.min
    
    screen_step <- glmnet(cbind(A_beta, A_gamma), tilde_y, alpha = 1, lambda = best_lambda, penalty.factor = penalty)
    coefs <- screen_step$beta
    betas <- coefs[1:p, ]
    gammas <- coefs[(p+1):(p+m), ]
    
    S <- which(gammas != 0)
    
    list(S = S, gammas = gammas, betas = betas, best_nu = best_nu, best_lam = best_lambda)
    
  }, error = function(e) {
    print("error")
    list(S = NULL, gammas = NULL, betas = NULL, best_nu = NULL, best_lam = NULL)
  })
}

sk.create <- function(X, y, D, nu) {
  n <- nrow(X)
  m <- nrow(D)
  
  A_beta <- rbind(X/sqrt(2*n),D/sqrt(2*nu))
  A_gamma <- rbind(matrix(0,n,m),-diag(m)/sqrt(2*nu))
  
  tilde_y <- c(y/sqrt(2*n),rep(0, m))
  
  structure(list(call = match.call(),
                 A_beta = A_beta,
                 A_gamma = A_gamma,
                 tilde_y = tilde_y))
}
###############################################################################
###############################################################################
###############################################################################
## (2) HD_SplitKnockoff
HD_SplitKnockoff <- function(
    X, y, D, q = 0.1,
    seed        = NULL,
    fold        = 10,
    eta         = 0.1,
    method      = "knockoff+",
    normalize   = TRUE,
    copy        = TRUE,
    lambda_seq  = 10^seq(0, -6, by = -0.01),
    nu          = 10.^seq(0, 2, length.out = 10),
    sign        = "disabled"
) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  
  #—— sample splitting ——#
  sample_index1 <- sample(seq_len(n), size = floor(n/2), replace = FALSE)
  sample_index2 <- setdiff(seq_len(n), sample_index1)
  
  #—— screening step 1: Lasso ——#
  cvfit    <- cv.glmnet(
    x           = X[sample_index1, , drop = FALSE],
    y           = y[sample_index1],
    type.measure= "mse",
    nfolds      = fold
  )
  lambda_min <- cvfit$lambda.min
  beta1      <- as.vector(glmnet(
    x      = X[sample_index1, , drop = FALSE],
    y      = y[sample_index1],
    family = "gaussian",
    alpha  = 1,
    lambda = lambda_min
  )$beta)
  
  S_beta <- sort(which(beta1 != 0))
  k_beta <- length(sample_index2) - 1
  if (length(S_beta) > k_beta) {
    idx_order <- order(-abs(beta1))
    S_beta    <- sort(idx_order[1:k_beta])
  }
  
  #—— screening step 2: SplitLasso ——#
  split_lasso_res <- SplitLasso(
    X[sample_index1, S_beta, drop = FALSE],
    y[sample_index1],
    D[, S_beta, drop = FALSE],
    nu_values = nu
  )
  gamma1 <- split_lasso_res$gammas
  S_gamma <- sort(which(gamma1 != 0))
  if (length(S_gamma) == 0) {
    SplitKnockoff_Sel <- integer(0)
    message("No gamma selected in the screening step! Return Nothing.")
    return(list(Sel = SplitKnockoff_Sel))
  }
  
  k_gamma <- length(sample_index2) - length(S_beta)
  if (length(S_gamma) > k_gamma) {
    idx_order <- order(abs(gamma1), decreasing = TRUE)
    S_gamma   <- idx_order[1:k_gamma]
  }
  
  D_scr <- D[S_gamma, S_beta, drop = FALSE]
  
  # print(dim(D_scr))
  # print(all(which(abs(D %*% Beta0)> 1e-4) %in% S_gamma)) ## Sure screening of gamma
  # print(all(which(abs(Beta0)> 1e-4) %in% S_beta))        ## Sure screening of beta
  # print(which(abs(D %*% Beta0)> 1e-4))                             ## Original active set S
  # print(S_gamma[which(abs(D_scr %*% Beta0[S_beta])> 1e-4)])    ## S after screening
  
  #—— selection step: SplitKnockoff ——#
  option <- list(data = NA, dim = length(data), dimnames = NULL)
  option$q <- q
  # choice on threshold
  option$method <- method
  # degree of separation between original design and its split knockoff copy
  # in the range of [0, 2], the less the more separated
  option$eta <- eta
  # whether to normalize the dataset
  option$normalize <- ifelse(normalize, "true", "false")
  # whether to create a knockoff copy
  option$copy <- ifelse(copy,      "true", "false")
  # choice on the set of regularization parameters for split LASSO path
  option$lambda <- lambda_seq
  # choice of nu for split knockoffs
  option$nu <- nu
  # choice on whether to estimate the directional effect, 'disabled'/'enabled'
  option$sign <- sign
  option <- option[-1]
  
  SplitKnockoff_Fit <- tryCatch({
    result <- sk.filter(
      X[sample_index2, S_beta, drop = FALSE],
      D_scr,
      matrix(y[sample_index2], ncol = 1),
      option
    )
    result
  }, error = function(e) {
    message("Error in SplitKnockoff: ", conditionMessage(e))
    list(selected = NULL)
  })

  SplitKnockoff_Sel <- S_gamma[SplitKnockoff_Fit$results$data]
  
  return(list(
    Sel       = SplitKnockoff_Sel,
    S_beta    = S_beta,
    S_gamma   = S_gamma,
    D_scr = D_scr,
    Fit       = SplitKnockoff_Fit
  ))
}