###############################################################################
# StatKnock: Stability Transformational Knockoff
# -----------------------------------------------------------------------------
# This script implements the StatKnock method for selection in structural sparsity problems.
#
# Inputs:
#   X (n x p), y (n), D (m x p);
#   L (# of stability splits), q (target FDR),
#   cv/cvfolds/lambda_cr (Lasso tuning), X_N_option/proj (augmentation choices),
#   knockoff_args (e.g., list(fixed_method = "equi")).
#
# Outputs:
#   list(S       = selected set via knockoff,
#        S_plus  = selected set via knockoff+,
#        W       = knockoff statistics,
#        T1, T_plus = thresholds,
#        X_star, y_star = transformed design/response)
###############################################################################
library(knockoff)
library(glmnet)
library(MASS)
library(reticulate)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
StatKnock <- function(X, y, D, L = 100, q = 0.1, 
                      cv = TRUE, cvfolds = 5, lambda_cr = c("min","1se"), 
                      X_N_option = 1, proj = TRUE,
                      knockoff_args = list()){
  lambda_cr <- match.arg(lambda_cr)
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- dim(D)[1]
  qrD <- qr(D)
  r <- qrD$rank
  
  if (r == m){
    D_plus = D
    X_plus = X
  } else {
    message("rank(D) < m, auxiliary noises introduced!")
    Q_full <- qr.Q(qrD, complete = TRUE)
    N <- Q_full[, (r + 1):m, drop = FALSE]
    D_plus <- cbind(D, N)
    
    if (X_N_option == 1) {
      X_N <- matrix(rnorm(n * (m - r)), nrow = n, ncol = m - r)
    } else if (X_N_option == 2) {
      X_N <- matrix(sample(c(-1, 1), n * (m - r), replace = TRUE), nrow = n, ncol = m - r)
    } else if (X_N_option == 3) {
      if ((m - r) > p) stop("Option 3 requires (m - r) <= p")
      X_N <- sapply(sample(ncol(X), m - r, replace = FALSE), function(j) sample(X[, j]))
    } else {
      stop("X_N_option must be 1, 2, or 3.")
    }
    
    if (proj) {
      if (n >= p) {
        P_X <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
        X_N <- P_X %*% X_N
      } else {
        message("Projection skipped: n < p")
      }
    }
    
    X_plus = cbind(X, X_N) 
  }
  
  if (r == p){
    X_star = X_plus %*% solve(D_plus)
    y_star = y
  } else {
    message("rank(D) < p, projection used!")
    E <- pracma::nullspace(D_plus)
    D_tilde <- rbind(D_plus, t(E))
    D_tilde_inverse <- solve(D_tilde)
    D_dagger <- D_tilde_inverse[, 1:m]
    D_star <- D_tilde_inverse[, (m + 1):(m + (p - r))]
    XDstar <- X_plus %*% D_star
    M_XDstar <- diag(n) - XDstar %*% solve(t(XDstar) %*% XDstar) %*% t(XDstar)
    y_star <- M_XDstar %*% y
    X_star <- M_XDstar %*% X_plus %*% D_dagger
  }
  
  res <- generate_knockoff(
    X_star,
    method_args = knockoff_args
  )
  X_org   <- res$X_org
  X_tilde <- res$X_tilde
  X_aug   <- cbind(X_org, X_tilde)

  S_hat <- rep(0, 2*m)
  if(cv == TRUE){
    cv.fit = cv.glmnet(x = X_aug, y = y_star, family = "gaussian", intercept = F, alpha = 1, nfolds = cvfolds, standardize = F)
    if (lambda_cr == "min") {
      best_lambda = cv.fit$lambda.min
    } else {
      best_lambda = cv.fit$lambda.1se
    }
  } else {
    best_lambda = sqrt(log(2*m)/n) / sqrt(n)

  }
  
  for(i in 1:L){
    index <- sample((1:n), floor(n/2), replace = F)
    
    fit1 = glmnet(x = X_aug[index, ], y = y_star[index], family = "gaussian", intercept = F, alpha = 1, lambda = best_lambda, standardize = F)
    bhat1 <- as.vector(stats::coef(fit1)[-1])
    S_index1 <- which(bhat1!=0)

    fit2 = glmnet(x = X_aug[-index, ],y = y_star[-index], family = "gaussian", intercept = F, alpha = 1, lambda = best_lambda, standardize = F)
    bhat2 <- as.vector(stats::coef(fit2)[-1])
    S_index2 <- which(bhat2!=0)
    
    S_index <- intersect(S_index1, S_index2)
    
    S = rep(0, 2*m)
    S[S_index]=1
    
    S_hat = S_hat + S
  }
  S_hat <- S_hat/L

  # compute Knockoff statistics
  W <- S_hat[1:m]-S_hat[(m+1):(2*m)]
  
  # find the knockoff threshold T
  t = sort(unique(abs(W)[abs(W) > 0]))
  ratio = numeric(length(t))
  
  for (j in seq_along(t)) {
    ratio[j] = (sum(W <= -t[j]))/max(1, sum(W >= t[j]))
  }
  id = which(ratio <= q)[1]
  if (is.na(id)) {
    T1 = Inf
  } else {
    T1 = t[id]
  }
  
  # find the knockoff_plus threshold T_plus
  t_plus = sort(unique(abs(W)[abs(W) > 0]))
  ratio_plus = numeric(length(t_plus))
  
  for (j in seq_along(t_plus)) {
    ratio_plus[j] = (1 + sum(W <= -t_plus[j]))/max(1, sum(W >= t_plus[j]))
  }
  id_plus = which(ratio_plus <= q)[1]
  if (is.na(id_plus)) {
    T_plus = Inf
  } else {
    T_plus = t_plus[id_plus]
  }
  
  # set of discovered variables
  S <- which(W >= T1)
  S_plus <- which(W >= T_plus)
  
  return(list(S = S, S_plus = S_plus, W = W, T1 = T1, T_plus = T_plus,
              X_star = X_star, y_star = y_star))
}
###############################################################################
###############################################################################
###############################################################################
### Helper function
generate_knockoff <- function(X_org, method_args = list()) {
  X_org <- scale(X_org, center = FALSE, scale = sqrt(colSums(X_org^2)))[,]
  fixed_method <- if (!is.null(method_args$fixed_method)) {
    method_args$fixed_method
  } else {
    "equi"
  }
  knock <- create.fixed(X_org, method = fixed_method)
  
  list(
    X_org   = knock$X,
    X_tilde = knock$Xk
  )
}