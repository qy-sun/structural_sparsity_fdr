###############################################################################
# TLasso: Transformational Lasso
# -----------------------------------------------------------------------------
# This script implements the TLasso method for structural sparsity recovery.
#
# Inputs:
#   X (n x p), y (n), D (m x p);
#   X_N_option ∈ {1,2,3} for auxiliary design generation; proj ∈ {TRUE,FALSE}.
#
# Outputs:
#   list(
#     gamma        = estimated Dβ,
#     beta         = recovered β on original p coordinates,
#     beta_N       = coefficients on auxiliary columns (if any),
#     gamma_debias = γ debiased by removing null-space component (if r < m),
#     best_lam     = CV-selected λ,
#     X_star, y_star = transformed design/response used for Lasso
#   )
###############################################################################
library(glmnet)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
TLasso <- function(X, y, D, X_N_option = 1, proj = TRUE) {
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- nrow(D)
  qrD <- qr(D)
  r <- qrD$rank

  if (r == m){
    D_plus = D
    X_plus = X
  } else {
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
  
  cv_fit <- cv.glmnet(X_star, y_star, family = "gaussian", alpha = 1, nfolds = 5, standardize = T, intercept = F)
  best_lambda <- cv_fit$lambda.min
  final_fit <- glmnet(X_star, y_star, family = "gaussian", alpha = 1, lambda = best_lambda, standardize = T, intercept = F)
  gamma <- as.vector(final_fit$beta)
  
  # recover beta
  if (r == p) {
    beta_plus <- solve(D_plus) %*% gamma
  } else {
    beta_partial <- D_dagger %*% gamma
    residual_y <- y - X_plus %*% beta_partial
    alpha <- solve(t(X_plus %*% D_star) %*% (X_plus %*% D_star)) %*% t(X_plus %*% D_star) %*% residual_y
    beta_plus <- beta_partial + D_star %*% alpha
  }
  beta = beta_plus[1:p]
  beta_N = beta_plus[(p+1):length(beta_plus)]
  
  # debias gammma (optional)
  if (r == m){
    gamma_debias = gamma
  } else {
    gamma_debias = gamma - N %*% beta_N
  }
  
  return(list(gamma = gamma, 
              beta = beta, 
              beta_N = beta_N, 
              gamma_debias = gamma_debias,
              best_lam = best_lambda,
              X_star = X_star, 
              y_star = y_star))
}