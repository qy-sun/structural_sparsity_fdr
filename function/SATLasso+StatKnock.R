###############################################################################
# SATLasso + StatKnock
# -----------------------------------------------------------------------------
# This script implements the following high-dimensional two-stage pipeline:
#   1) Screening: SATLasso on (X1, y1) to select constraints (gamma) and
#      recycle connected coefficients (beta) via D.
#   2) Selection: StatKnock on (X2, y2) with the screened design (X2[, S_beta])
#      and D_scr = D[S_gamma, S_beta], controlling FDR.
#
# Inputs:
#   X (n x p), y (n x 1), D (m x p);  L, q, cv, cvfolds, lambda_cr, knockoff_args.
#
# Outputs:
#   list(S = selected indices in original gamma space,
#        S_plus = conservative (knockoff+) selection if available)
###############################################################################
this.dir <- dirname(sys.frame(1)$ofile)
source(file.path(this.dir, "StatKnock.R"))
source(file.path(this.dir, "SATLasso.R"))
###############################################################################
###############################################################################
###############################################################################
SATLasso_StatKnock <- function(X, y, D,
                               L = 100, q = 0.1,
                               cv = TRUE, cvfolds = 10,
                               lambda_cr = c("min", "1se"),
                               knockoff_args = list()) {
  lambda_cr <- match.arg(lambda_cr)
  y <- matrix(y, ncol = 1)
  n <- nrow(X); p <- ncol(X); m <- nrow(D)
  
  ## If low-dimensional, skip screening
  if (n >= 2 * m) {
    message("Low-dimensional setting: running StatKnock directly.")
    return(
      StatKnock(X, y, D, L = L, q = q,
                cv = cv, cvfolds = cvfolds,
                lambda_cr = lambda_cr,
                knockoff_args = knockoff_args)
    )
  }
  
  #############################################
  ## 1) Random split into (X1, y1) and (X2, y2)
  #############################################
  indices <- sample.int(n)
  n1 <- floor(n / 2); n2 <- n - n1
  idx1 <- indices[1:n1]
  idx2 <- indices[(n1 + 1):n]
  X1 <- X[idx1, , drop = FALSE];  y1 <- y[idx1, , drop = FALSE]
  X2 <- X[idx2, , drop = FALSE];  y2 <- y[idx2, , drop = FALSE]
  
  #############################################
  ## 2) Screening on (X1, y1)
  #############################################
  
  # ## Conditions to ensure that StatKnock works
  # max_gamma = floor(n2/2)
  # max_beta  = n2 - 1
  
  # (i) S_gamma: screening gamma via SATLasso
  S_gamma <- SATLasso(X1, y1, D, stop = floor(n2/log(n2)))$em
  
  # (ii) S_beta: recycling beta based on S_gamma (vertices connected to the screened edges)
  S_beta  <- which(colSums(abs(D[S_gamma, , drop = FALSE]) != 0) > 0)
  
  ## (v) Form screened D matrix
  D_scr <- D[S_gamma, S_beta]
  
  # print(dim(D_scr))
  # print(all(which(abs(D %*% Beta0)> 1e-10) %in% S_gamma)) ## Sure screening of gamma
  # print(all(which(abs(Beta0)> 1e-10) %in% S_beta))        ## Sure screening of beta
  # print(which(D %*% Beta0!=0))                             ## Original active set S
  # print(S_gamma[which(D_scr %*% Beta0[S_beta]!=0)])    ## S after screening
  
  #############################################
  ## 3) Selection on (X2, y2)
  #############################################
  if (is.null(dim(D_scr))) {
    Sel <- integer(0)
    Sel_plus <- integer(0)
  } else {
    stat_fit <- StatKnock(X2[, S_beta, drop = FALSE],
                          y2, D = D_scr,
                          L = L, q = q,
                          cv = cv, cvfolds = cvfolds,
                          lambda_cr = lambda_cr,
                          knockoff_args = knockoff_args)
    Sel      <- S_gamma[stat_fit$S]
    Sel_plus <- S_gamma[stat_fit$S_plus]
  }
  
  return(list(S=Sel, S_plus=Sel_plus))
}
