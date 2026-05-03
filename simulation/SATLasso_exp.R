###############################################################################
# Simulation Experiments for SATLasso
# -----------------------------------------------------------------------------
# This script reproduces the results in Table S.1 and Table S.2 in the Supplement
###############################################################################
################################################################################
################################################################################
################################################################################
source("helper/Graph Difference Operator.R")
source("helper/Beta Structure Generation.R")
source("../function/SATLasso.R")
source("../function/other methods/HD_SplitKnockoff.R")
source("../function/other methods/HD_GKnockoff.R")
library(genlasso)
################################################################################
################################################################################
################################################################################
### Helper function
fdp_power_ssr <- function(selected_index, signal_index){
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
mean.se = function(x, ns = 1) { paste0(format(round(mean(x)*100, ns), nsmall=ns), ' (', format(round(100*sd(x)/sqrt(length(x)), ns), nsmall=ns), ')')}
################################################################################
################################################################################
################################################################################
set.seed(2025)
nIterations = 100
################################################################################
################################################################################
################################################################################
### Table S.1-1: piecewise constant over Chain

set.seed(2025)
p = 1000; n = 200; A = 1; rho = 0.2
p1 = 1000
m1 = 20
Sigma= toeplitz(rho^(1:p-1))
D <- create.chain(p)
m <- nrow(D)

MMS_satlasso_ChainConstant = list()
MMS_fusis_ChainConstant = list()
MMS_genLasso_ChainConstant = list()
MMS_splitlasso_ChainConstant = list()

for (iter in 1:nIterations) {
  print(iter)
  ### Generate the Data
  Beta0 = Chain_PiecewiseConstant(p,p1,m1,A,min_interval = 5, edge_buffer = 5)
  X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)
  signal_ind = which(abs(D%*%Beta0) > 1e-10)
  print(length(signal_ind))
  
  ## SATLasso
  SATLasso_Fit <- SATLasso(X, y, D, signal_ind = signal_ind)
  MMS_satlasso_ChainConstant[[iter]] <- length(SATLasso_Fit$es)
  print(MMS_satlasso_ChainConstant[[iter]])
  
  # FuSIS
  fusisSel = FuSIS(X, y, D)
  for (j2 in length(signal_ind):m) {
    if (length(intersect(fusisSel[1:j2],signal_ind))==length(signal_ind)){
      MMS_fusis_ChainConstant[[iter]] = j2
      break
    }
  }
  print(MMS_fusis_ChainConstant[[iter]])
  
  ## GenLasso
  GenLasso_Fit <- genlasso(y, X, D)
  gamma = D %*% coef(GenLasso_Fit, lambda=sqrt(n*log(p)))$beta
  GenLasso_Sel <- order(-abs(gamma))
  for (j3 in length(signal_ind):m) {
    if (length(intersect(GenLasso_Sel[1:j3],signal_ind))==length(signal_ind)){
      MMS_genLasso_ChainConstant[[iter]] = j3
      break
    }
  }
  print(MMS_genLasso_ChainConstant[[iter]])
  
  ## SplitLasso
  SplitLasso_Fit <- SplitLasso(X, y, D)
  SplitLasso_Sel <- order(-abs(SplitLasso_Fit$gammas))
  for (j4 in length(signal_ind):m) {
    if (length(intersect(SplitLasso_Sel[1:j4],signal_ind))==length(signal_ind)){
      MMS_splitlasso_ChainConstant[[iter]] = j4
      break
    }
  }
  print(MMS_splitlasso_ChainConstant[[iter]])
}
quantile(unlist(MMS_satlasso_ChainConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
quantile(unlist(MMS_fusis_ChainConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
quantile(unlist(MMS_genLasso_ChainConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
quantile(unlist(MMS_splitlasso_ChainConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
save(MMS_satlasso_ChainConstant, file = "MMS_satlasso_ChainConstant.RData")
save(MMS_fusis_ChainConstant, file = "MMS_fusis_ChainConstant.RData")
save(MMS_genLasso_ChainConstant, file = "MMS_genLasso_ChainConstant.RData")
save(MMS_splitlasso_ChainConstant, file = "MMS_splitlasso_ChainConstant.RData")
################################################################################
################################################################################
################################################################################
### Table S.1-2: piecewise constant over Lattice

set.seed(2025)
p = 1000; n = 200; A = 1; rho = 0.2
p1 = 1000
s = 25
t = 40
m1 = 20
Sigma= toeplitz(rho^(1:p-1))
D <- create.lattice(s,t)
m <- nrow(D)

MMS_satlasso_LatticeConstant = list()
MMS_fusis_LatticeConstant = list()
MMS_genLasso_LatticeConstant = list()
MMS_splitlasso_LatticeConstant = list()

for (iter in 1:nIterations) {
  print(iter)
  ### Generate the Data
  Beta0 = Lattice_PiecewiseConstant(s, t, p1, m1, A)
  X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)
  signal_ind = which(abs(D%*%Beta0) > 1e-10)
  
  ## SATLasso
  SATLasso_Fit <- SATLasso(X, y, D, signal_ind = signal_ind)
  MMS_satlasso_LatticeConstant[[iter]] <- length(SATLasso_Fit$es)
  print(MMS_satlasso_LatticeConstant[[iter]])
  
  # FuSIS
  fusisSel = FuSIS(X, y, D)
  for (j2 in length(signal_ind):m) {
    if (length(intersect(fusisSel[1:j2],signal_ind))==length(signal_ind)){
      MMS_fusis_LatticeConstant[[iter]] = j2
      break
    }
  }
  print(MMS_fusis_LatticeConstant[[iter]])
  
  ## GenLasso
  GenLasso_Fit <- genlasso(y, X, D)
  gamma = D %*% coef(GenLasso_Fit, lambda=sqrt(n*log(p)))$beta
  GenLasso_Sel <- order(-abs(gamma))
  for (j3 in length(signal_ind):m) {
    if (length(intersect(GenLasso_Sel[1:j3],signal_ind))==length(signal_ind)){
      MMS_genLasso_LatticeConstant[[iter]] = j3
      break
    }
  }
  print(MMS_genLasso_LatticeConstant[[iter]])
  
  ## SplitLasso
  SplitLasso_Fit <- SplitLasso(X, y, D)
  SplitLasso_Sel <- order(-abs(SplitLasso_Fit$gammas))
  for (j4 in length(signal_ind):m) {
    if (length(intersect(SplitLasso_Sel[1:j4],signal_ind))==length(signal_ind)){
      MMS_splitlasso_LatticeConstant[[iter]] = j4
      break
    }
  }
  print(MMS_splitlasso_LatticeConstant[[iter]])
}
quantile(unlist(MMS_satlasso_LatticeConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
quantile(unlist(MMS_fusis_LatticeConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
quantile(unlist(MMS_genLasso_LatticeConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
quantile(unlist(MMS_splitlasso_LatticeConstant),c(0.05,0.25,0.5,0.75,0.95),type = 1)
save(MMS_satlasso_LatticeConstant, file = "MMS_satlasso_LatticeConstant.RData")
save(MMS_fusis_LatticeConstant, file = "MMS_fusis_LatticeConstant.RData")
save(MMS_genLasso_LatticeConstant, file = "MMS_genLasso_LatticeConstant.RData")
save(MMS_splitlasso_LatticeConstant, file = "MMS_splitlasso_LatticeConstant.RData")
################################################################################
################################################################################
################################################################################
### Table S.2: 1D sparse fused Lasso

set.seed(2025)
p = 1000
n = 200
p1 = 20
m1 = 5

A = 1
rho = 0.2
Sigma= toeplitz(rho^(1:p-1))

D_1 <- diag(p)
m_1 <- nrow(D_1)
D_2 <- create.chain(p)
m_2 <- nrow(D_2)
D_3 <- rbind(diag(p), create.chain(p))
m_3 <- nrow(D_3)

fdp_power_list_1D_SFL = sapply(1:nIterations, function(it) {
  ### Generate the Data
  Beta0 = Chain_PiecewiseConstant(p,p1,m1,A, min_interval = 3, edge_buffer = 3)
  X <- matrix(rnorm(n*p),n) %*% chol(Sigma)
  y <- X %*% Beta0 + rnorm(n, sd = 1)
  signal_ind_1 = which(abs(D_1%*%Beta0) > 1e-10)
  signal_ind_2 = which(abs(D_2%*%Beta0) > 1e-10)
  signal_ind_3 = which(abs(D_3%*%Beta0) > 1e-10)
  
  ### SATLasso
  SATLasso_Fit_1 <- SATLasso(X, y, D_1, stop = floor(n/log(n)))
  plot(SATLasso_Fit_1$bic)
  SATLasso_Sel_1 <- SATLasso_Fit_1$em
  fdp_power1_1 = fdp_power_ssr(selected_index=SATLasso_Sel_1,signal_index=signal_ind_1)
  SATLasso_Fit_2 <- SATLasso(X, y, D_2, stop = floor(n/log(n)))
  plot(SATLasso_Fit_2$bic)
  SATLasso_Sel_2 <- SATLasso_Fit_2$em
  fdp_power1_2 = fdp_power_ssr(selected_index=SATLasso_Sel_2,signal_index=signal_ind_2)
  SATLasso_Fit_3 <- SATLasso(X, y, D_3, stop = floor(n/log(n)))
  plot(SATLasso_Fit_3$bic)
  SATLasso_Sel_3 <- SATLasso_Fit_3$em
  fdp_power1_3 = fdp_power_ssr(selected_index=SATLasso_Sel_3,signal_index=signal_ind_3)
  
  ### FuSIS
  fusisSel_1 <- FuSIS(X, y, D_1, stop = floor(n/log(n)))
  fdp_power2_1 = fdp_power_ssr(selected_index=fusisSel_1,signal_index=signal_ind_1)
  fusisSel_2 <- FuSIS(X, y, D_2, stop = floor(n/log(n)))
  fdp_power2_2 = fdp_power_ssr(selected_index=fusisSel_2,signal_index=signal_ind_2)
  fusisSel_3 <- FuSIS(X, y, D_3, stop = floor(n/log(n)))
  fdp_power2_3 = fdp_power_ssr(selected_index=fusisSel_3,signal_index=signal_ind_3)
  
  #### GLasso
  GLasso_Fit_1 <- genlasso(y, X, D_1)
  gamma_genLasso_1 = D_1 %*% coef(GLasso_Fit_1, lambda=sqrt(n*log(p)))$beta
  GLasso_Sel_1 <- which(gamma_genLasso_1!=0)
  fdp_power3_1 = fdp_power_ssr(selected_index=GLasso_Sel_1,signal_index=signal_ind_1)
  GLasso_Fit_2 <- genlasso(y, X, D_2)
  gamma_genLasso_2 = D_2 %*% coef(GLasso_Fit_2, lambda=sqrt(n*log(p)))$beta
  GLasso_Sel_2 <- which(gamma_genLasso_2!=0)
  fdp_power3_2 = fdp_power_ssr(selected_index=GLasso_Sel_2,signal_index=signal_ind_2)
  GLasso_Fit_3 <- genlasso(y, X, D_3)
  gamma_genLasso_3 = D_3 %*% coef(GLasso_Fit_3, lambda=sqrt(n*log(p)))$beta
  GLasso_Sel_3 <- which(gamma_genLasso_3!=0)
  fdp_power3_3 = fdp_power_ssr(selected_index=GLasso_Sel_3,signal_index=signal_ind_3)
  
  ### SplitLasso
  SplitLasso_Fit_1 <- SplitLasso(X, y, D_1)
  SplitLasso_Sel_1 <- SplitLasso_Fit_1$S
  fdp_power4_1 = fdp_power_ssr(selected_index=SplitLasso_Sel_1,signal_index=signal_ind_1)
  SplitLasso_Fit_2 <- SplitLasso(X, y, D_2)
  SplitLasso_Sel_2 <- SplitLasso_Fit_2$S
  fdp_power4_2 = fdp_power_ssr(selected_index=SplitLasso_Sel_2,signal_index=signal_ind_2)
  SplitLasso_Fit_3 <- SplitLasso(X, y, D_3)
  SplitLasso_Sel_3 <- SplitLasso_Fit_3$S
  fdp_power4_3 = fdp_power_ssr(selected_index=SplitLasso_Sel_3,signal_index=signal_ind_3)
  
  
  r = c(fdp_power1_1$fdp, fdp_power1_1$power, fdp_power1_1$ssr,
        fdp_power1_2$fdp, fdp_power1_2$power, fdp_power1_2$ssr,
        fdp_power1_3$fdp, fdp_power1_3$power, fdp_power1_3$ssr,
        fdp_power2_1$fdp, fdp_power2_1$power, fdp_power2_1$ssr,
        fdp_power2_2$fdp, fdp_power2_2$power, fdp_power2_2$ssr,
        fdp_power2_3$fdp, fdp_power2_3$power, fdp_power2_3$ssr,
        fdp_power3_1$fdp, fdp_power3_1$power, fdp_power3_1$ssr,
        fdp_power3_2$fdp, fdp_power3_2$power, fdp_power3_2$ssr,
        fdp_power3_3$fdp, fdp_power3_3$power, fdp_power3_3$ssr,
        fdp_power4_1$fdp, fdp_power4_1$power, fdp_power4_1$ssr,
        fdp_power4_2$fdp, fdp_power4_2$power, fdp_power4_2$ssr,
        fdp_power4_3$fdp, fdp_power4_3$power, fdp_power4_3$ssr
  )
})
save(fdp_power_list_1D_SFL, file = "SATLasso_1D_SFL.RData")