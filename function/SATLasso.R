###############################################################################
# SATLasso: Sequential Adjusted Transformational Lasso
# -----------------------------------------------------------------------------
# This script implements the SATLasso method for screening in structural sparsity problems.
#
# Reference:
#   Luo, S., & Chen, Z. (2014).
#   Sequential Lasso cum EBIC for feature selection with ultra-high dimensional
#   feature space. Journal of the American Statistical Association, 109(507), 1229–1240.
###############################################################################
library(glmnet)
library(glmpath)
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
SATLasso <- function(X,y,D, alpha_grid = 10^seq(-3, 0, length.out = 30),
                     cv_alpha_folds = 5,
                     cv_lambda_folds = 5,
                     stop = NULL, signal_ind = NULL){
  n <- nrow(X)
  p <- ncol(X)
  m <- nrow(D)
  qrD <- qr(D)
  r <- qrD$rank
  
  if (r == m){
    final_res <- STLasso(X, y, D, stop = stop, signal_ind = signal_ind)
    return(final_res)
  }
  
  folds_alpha <- sample(rep(seq_len(cv_alpha_folds), length.out = n))
  
  alpha_mse <- numeric(length(alpha_grid))
  
  Q_full  <- qr.Q(qrD, complete = TRUE)
  Nmat    <- Q_full[, (r + 1):m, drop = FALSE]
  D_plus  <- cbind(D, Nmat)
  
  for (i in seq_along(alpha_grid)) {
    alpha <- alpha_grid[i]
    fold_errors <- numeric(cv_alpha_folds)
    
    for (f in seq_len(cv_alpha_folds)) {
      train_idx <- which(folds_alpha != f)
      test_idx  <- which(folds_alpha == f)
      X_train   <- X[train_idx, , drop = FALSE]
      y_train   <- y[train_idx]
      X_test    <- X[test_idx,  , drop = FALSE]
      y_test    <- y[test_idx]
      n_train   <- length(y_train)
      n_test    <- length(y_test)
      
      X_N_train <- matrix(0, nrow = length(train_idx), ncol = (m - r))
      X_N_test <- matrix(0, nrow = nrow(X_test), ncol = (m - r))
      
      y_aug_train <- rbind(matrix(1/sqrt(2*n_train) * y_train, ncol = 1),
                           matrix(0, m - r, 1))
      X_aug_train <- rbind(
        cbind(1/sqrt(2*n_train) * X_train, 1/sqrt(2*n_train) * X_N_train),
        cbind(matrix(0, m - r, p), sqrt(alpha) * diag(m - r))
      )
      
      y_aug_test <- rbind(matrix(1/sqrt(2*n_test) * y_test, ncol = 1),
                          matrix(0, m - r, 1))
      X_aug_test <- rbind(
        cbind(1/sqrt(2*n_test) * X_test, 1/sqrt(2*n_test) * X_N_test),
        cbind(matrix(0, m - r, p), sqrt(alpha) * diag(m - r))
      )
      
      if (r == p) {
        X_star_train <- X_aug_train %*% MASS::ginv(D_plus)
        y_star_train <- y_aug_train
        X_star_test <- X_aug_test %*% MASS::ginv(D_plus)
        y_star_test <- y_aug_test
      } else {
        E       <- pracma::nullspace(D_plus)
        D_tilde <- rbind(D_plus, t(E))
        D_inv   <- MASS::ginv(D_tilde)
        
        D_dagger <- D_inv[, 1:m]
        D_star   <- D_inv[, (m + 1):(m + (p - r))]
        
        XDstar_train   <- X_aug_train %*% D_star
        Mmat_train     <- diag(nrow(X_aug_train)) -
          XDstar_train %*% MASS::ginv(t(XDstar_train) %*% XDstar_train) %*% t(XDstar_train)
        y_star_train <- Mmat_train %*% y_aug_train
        X_star_train <- Mmat_train %*% X_aug_train %*% D_dagger
        
        XDstar_test   <- X_aug_test %*% D_star
        Mmat_test     <- diag(nrow(X_aug_test)) -
          XDstar_test %*% MASS::ginv(t(XDstar_test) %*% XDstar_test) %*% t(XDstar_test)
        y_star_test <- Mmat_test %*% y_aug_test
        X_star_test <- Mmat_test %*% X_aug_test %*% D_dagger
      }
      
      cv_fit <- cv.glmnet(X_star_train, y_star_train,
                          family = "gaussian", alpha = 1,
                          nfolds = cv_lambda_folds,
                          standardize = FALSE,
                          intercept = FALSE)
      best_lambda <- cv_fit$lambda.min
      final_fit <- glmnet(X_star_train, y_star_train,
                          family = "gaussian", alpha = 1,
                          lambda = best_lambda,
                          standardize = FALSE,
                          intercept = FALSE)
      gamma_hat <- as.vector(final_fit$beta)
      
      y_star_pred <- X_star_test %*% gamma_hat
      fold_errors[f] <- mean( (y_star_test - y_star_pred)^2 )
    }
    
    alpha_mse[i] <- mean(fold_errors)
  }
  
  best_alpha_idx <- which.min(alpha_mse)
  best_alpha     <- alpha_grid[best_alpha_idx]
  print(paste("best alpha =", best_alpha))
  
  final_res <- SATLasso_alpha(X=X, y=y, D=D, 
                              alpha = best_alpha,
                              stop = stop,
                              signal_ind = signal_ind)
  return(final_res)
}
###############################################################################
###############################################################################
###############################################################################
### Helper functions
## (1) STLasso
STLasso <- function(X,y,D, X_N_option = 1, proj = TRUE, stop = NULL, signal_ind = NULL){
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
        P_X <- diag(n) - X %*% MASS::ginv(t(X) %*% X) %*% t(X)
        X_N <- P_X %*% X_N
      } else {
        message("Projection skipped: n < p")
      }
    }
    
    X_plus = cbind(X, X_N) 
  }
  
  if (r == p){
    X_star = X_plus %*% MASS::ginv(D_plus)
    y_star = y
  } else {
    E <- pracma::nullspace(D_plus)
    D_tilde <- rbind(D_plus, t(E))
    D_tilde_inverse <- MASS::ginv(D_tilde)
    D_dagger <- D_tilde_inverse[, 1:m]
    D_star <- D_tilde_inverse[, (m + 1):(m + (p - r))]
    XDstar <- X_plus %*% D_star
    M_XDstar <- diag(n) - XDstar %*% MASS::ginv(t(XDstar) %*% XDstar) %*% t(XDstar)
    y_star <- M_XDstar %*% y
    X_star <- M_XDstar %*% X_plus %*% D_dagger
  }
  
  cb <- NULL
  if (!is.null(signal_ind)) {
    cb <- function(path) {
      length(intersect(path, signal_ind)) == length(signal_ind)
    }
  }
  
  final_res <- SLasso(X_star, y_star, stop = stop, callback = cb)
  return(final_res)
}
###############################################################################
###############################################################################
###############################################################################
## (2) SATLasso_alpha (fix alpha)
SATLasso_alpha <- function(X,y,D, alpha = 0.1, stop = NULL, signal_ind = NULL){
  y <- matrix(y, ncol = 1)
  n <- dim(X)[1]
  p <- dim(X)[2]
  m <- nrow(D)
  qrD <- qr(D)
  r <- qrD$rank
  
  if (r == m){
    final_res <- STLasso(X, y, D, stop = stop, signal_ind = signal_ind)
    return(final_res)
  }
  
  Q_full <- qr.Q(qrD, complete = TRUE)
  N <- Q_full[, (r + 1):m, drop = FALSE]
  D_plus <- cbind(D, N)
  
  X_N <- matrix(0, n, (m-r))
  
  y_aug <- rbind(1/sqrt(2*n) * y, matrix(0, m-r, 1))
  X_aug <- rbind(
    cbind(1/sqrt(2*n) * X,     1/sqrt(2*n) * X_N),
    cbind(matrix(0, m-r, p),  sqrt(alpha) * diag(m-r))
  )
  
  if (r == p) {
    X_star <- X_aug %*% MASS::ginv(D_plus)
    y_star <- y_aug
  } else {
    E        <- pracma::nullspace(D_plus)
    D_tilde  <- rbind(D_plus, t(E))
    D_inv    <- MASS::ginv(D_tilde)
    
    D_dagger <- D_inv[, 1:m]
    D_star   <- D_inv[, (m + 1):(m + (p - r))]
    
    XDstar   <- X_aug %*% D_star
    M        <- diag(n+m-r) - XDstar %*% MASS::ginv(t(XDstar) %*% XDstar) %*% t(XDstar)
    
    y_star <- M %*% y_aug
    X_star <- M %*% X_aug %*% D_dagger
  }
  
  cb <- NULL
  if (!is.null(signal_ind)) {
    cb <- function(path) {
      length(intersect(path, signal_ind)) == length(signal_ind)
    }
  }
  
  final_res <- SLasso(X_star, y_star, stop = stop, callback = cb)
  return(final_res)
}
###############################################################################
###############################################################################
###############################################################################
## (3) SLasso (Luo & Chen, 2014)
SLasso <- function(X,y, stop = NULL, callback = NULL){
  n<-nrow(X)
  p<-ncol(X)
  eta<- 2.1 -log(n)/log(p)
  path<-setdiff(1:p,1:p)
  EBIC<-NULL
  
  ## Initial model
  index.m <- sequential_select_whole(X,y,p,path)
  path    <- c(path,index.m)
  W <- diag(rep(1,n)) - (X[,index.m] %*% t(X[,index.m])) / as.numeric(t(X[,index.m]) %*% X[,index.m])
  Ycheck <- W%*%y
  EBIC[1] <- n*log(t(Ycheck)%*%Ycheck/n) + length(path)*log(n) + eta*lchoose(p,length(path))

  k<-1
  while (k < p){
    ## Early Stop
    if (!is.null(callback) && callback(path)) {
      break
    }

    ## Sequential Select
    Xcheck<-W%*%X
    Ycheck<-W%*%y
    index.m<-sequential_select_whole(Xcheck,Ycheck,p,path)
    path<-c(path,index.m)

    if (!is.null(stop) && k == stop) break

    ## Update Projector
    W <- W %*% (diag(rep(1,n)) - (X[,index.m] %*% t(X[,index.m]) %*% W)
                / as.numeric(t(X[,index.m]) %*% W %*% X[,index.m]))

    Ycheck <- W%*%y
    EBIC[k+1] <- n*log(t(Ycheck)%*%Ycheck/n) + length(path)*log(n) + eta*lchoose(p,length(path))
    # if (EBIC[k+1]>EBIC[k]) break

    k<-k+1
  }
  bic <- which.min(EBIC)
  em  <- path[1:bic]
  
  return(list(es=path,
              em=em,
              bic=EBIC))
}

sequential_select_whole<- function(X, y, p, exist_path)
{
  non_zero_index <-c()
  original_indices <- setdiff(1:ncol(X), exist_path)
  X_new = X[, original_indices]
  model_beta <- glmnet(X_new, y, alpha = 1, standardize = FALSE, intercept = FALSE)$beta
  for (i in 1:dim(model_beta)[2]) {
    current_beta <- model_beta[, i]
    if (any(current_beta != 0)) {
      non_zero_index <- as.numeric(which(current_beta != 0)[1])
      non_zero_index <- original_indices[non_zero_index]
      break
    }
  }
  return(non_zero_index)
}