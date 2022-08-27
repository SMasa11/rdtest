#'
#' return max statistics for covariates balance
#'
#' @inheritParams return_result_balance
#'

compute_stat_ma <- function(df_data,
                             int_dim_Z,
                             int_J = 3,
                             bool_max_test = FALSE,
                             bool_max_test_inv_V = FALSE,
                             bool_L2_std = FALSE)
{
  real_mean_effective_N_Z <- 0

  vec_h_L <- rep(0,int_dim_Z)
  vec_h_R <- rep(0,int_dim_Z)
  vec_V_robust_L <- rep(0,int_dim_Z)
  vec_V_robust_R <- rep(0,int_dim_Z)

  mat_beta_L <- matrix(rep(0,int_dim_Z*3),int_dim_Z,3)
  mat_beta_R <- matrix(rep(0,int_dim_Z*3),int_dim_Z,3)

  vec_X_L <- df_data$vec_X[df_data$vec_X<0,drop=FALSE]
  vec_X_R <- df_data$vec_X[df_data$vec_X>=0,drop=FALSE]
  mat_dif_in_X_L <- matrix(rep(t(vec_X_L),length(vec_X_L)),length(vec_X_L),length(vec_X_L))
  mat_dif_in_X_L <- abs(mat_dif_in_X_L - t(mat_dif_in_X_L)) + diag(Inf,length(vec_X_L),length(vec_X_L))

  mat_index_close_L <- matrix(rep(0,length(vec_X_L)*int_J),length(vec_X_L),int_J)
  for (j in 1:int_J) {
    mat_index_CloseL[,j] <- max.col(-mat.difInXL)
    for (i in 1:length(vec.XL)) {mat.difInXL[i,mat.indexCloseL[i,j]] <- Inf}
  }
  # vec.indexClose1stL <- max.col(-mat.difInXL)
  # for (i in 1:length(vec.XL)) {mat.difInXL[i,vec.indexClose1stL[i]] <- Inf}
  # vec.indexClose2ndL <- max.col(-mat.difInXL)
  # for (i in 1:length(vec.XL)) {mat.difInXL[i,vec.indexClose2ndL[i]] <- Inf}
  # vec.indexClose3rdL <- max.col(-mat.difInXL)

  mat.difInXR <- matrix(rep(t(vec.XR),length(vec.XR)),length(vec.XR),length(vec.XR))
  mat.difInXR <- abs(mat.difInXR - t(mat.difInXR)) + diag(Inf,length(vec.XR),length(vec.XR))

  mat.indexCloseR <- matrix(rep(0,length(vec.XR)*int_J),length(vec.XR),int_J)
  for (j in 1:int_J) {
    mat.indexCloseR[,j] <- max.col(-mat.difInXR)
    for (i in 1:length(vec.XR)) {mat.difInXR[i,mat.indexCloseR[i,j]] <- Inf}
  }
  # vec.indexClose1stR <- max.col(-mat.difInXR)
  # for (i in 1:length(vec.XR)) {mat.difInXR[i,vec.indexClose1stR[i]] <- Inf}
  # vec.indexClose2ndR <- max.col(-mat.difInXR)
  # for (i in 1:length(vec.XR)) {mat.difInXR[i,vec.indexClose2ndR[i]] <- Inf}
  # vec.indexClose3rdR <- max.col(-mat.difInXR)

  # Arrays to store
  mat.XhL <- matrix(rep(0,length(vec.XL)*int_dim_Z),length(vec.XL),int_dim_Z)
  mat.XhR <- matrix(rep(0,length(vec.XR)*int_dim_Z),length(vec.XR),int_dim_Z)

  arr.PsiL <- array(rep(0,3*length(vec.XL)*int_dim_Z),dim=c(length(vec.XL),3,int_dim_Z))
  arr.PsiR <- array(rep(0,3*length(vec.XR)*int_dim_Z),dim=c(length(vec.XR),3,int_dim_Z))

  arr.GinvL <- array(rep(0,3^2*int_dim_Z),dim=c(3,3,int_dim_Z))
  arr.GinvR <- array(rep(0,3^2*int_dim_Z),dim=c(3,3,int_dim_Z))

  mat.residualNnL <- matrix(rep(0,length(vec.XL)*int_dim_Z),length(vec.XL),int_dim_Z)
  mat.residualNnR <- matrix(rep(0,length(vec.XR)*int_dim_Z),length(vec.XR),int_dim_Z)
  ###

  vec.statTZ <- rep(0,int_dim_Z)
  vec.seZ <- rep(0,int_dim_Z)


  for (d in 1:int_dim_Z)
  {
    # enforce the bandwidths are the same for both sides
    # Note: density bands choose different bands by default, but the value seems to be sensitive to the option. I'll keep them as their default.
    eval(parse(text = paste0(paste0("list.resultRDRobustZ <- rdrobust(y=df_data$vec.Z.",d),",x=df_data$vec.X,rho=1,bwselect='mserd')")))

    # Zstat
    real.statTZ <- list.resultRDRobustZ$z[3]
    vec.statTZ[d] <- real.statTZ
    vec.seZ[d] <- list.resultRDRobustZ$se[3]
    vec.hL[d] <- list.resultRDRobustZ$bws[1,1]
    vec.hR[d] <- list.resultRDRobustZ$bws[1,2]
    if (is.null(fun.adjustBand) == FALSE) {
      vec.hL[d] <- vec.hL[d]*fun.adjustBand(int_dim_Z)
      vec.hR[d] <- vec.hR[d]*fun.adjustBand(int_dim_Z)
      eval(parse(text = paste0(paste0("list.resultRDRobustZ <- rdrobust(y=df_data$vec.Z.",d),",x=df_data$vec.X,rho=1,h=c(vec.hL[d],vec.hR[d]))")))
    }
    real.effectiveNMeanZ <- real.effectiveNMeanZ + sum(list.resultRDRobustZ$N_h)/int_dim_Z

    eval(parse(text = paste0(paste0("list.resultRDRobustZBias <- rdrobust(y=df_data$vec.Z.",d),",x=df_data$vec.X,p=2,rho=1,h=c(vec.hL[d],vec.hR[d]))")))

    vec.VRobustL[d] <- list.resultRDRobustZ$V_rb_l[1,1]
    vec.VRobustR[d] <- list.resultRDRobustZ$V_rb_r[1,1]
    if (abs(real.statTZ) > qnorm(1- 0.025)) {
      bool.rejectNaiveNull <- TRUE
    }
    if (abs(real.statTZ) > qnorm(1 - 0.025/(int_dim_Z+bool.joint))){
      bool.rejectBonferroniNull <- TRUE
    }

    # For each d, compute X(h_d)
    vec.XhL <- df_data$vec.X[df_data$vec.X< 0,drop=FALSE]/vec.hL[d]
    vec.XhR <- df_data$vec.X[df_data$vec.X>=0,drop=FALSE]/vec.hR[d]

    mat.XhLinL <- cbind(rep(1,length(vec.XhL)),vec.XhL)
    mat.XhLinR <- cbind(rep(1,length(vec.XhR)),vec.XhR)

    mat.XhSqL <- cbind(rep(1,length(vec.XhL)),vec.XhL,vec.XhL^2)
    mat.XhSqR <- cbind(rep(1,length(vec.XhR)),vec.XhR,vec.XhR^2)

    mat.KhL <- diag(((1-abs(vec.XhL))*(abs(vec.XhL)<=1))/vec.hL[d])
    mat.KhR <- diag(((1-abs(vec.XhR))*(abs(vec.XhR)<=1))/vec.hR[d])

    mat.GinvLinL <- chol2inv(chol(crossprod(sqrt(mat.KhL)%*%mat.XhLinL)))*length(vec.XhL)
    mat.GinvLinR <- chol2inv(chol(crossprod(sqrt(mat.KhR)%*%mat.XhLinR)))*length(vec.XhR)

    mat.GinvSqL <- chol2inv(chol(crossprod(sqrt(mat.KhL)%*%mat.XhSqL)))*length(vec.XhL)
    mat.GinvSqR <- chol2inv(chol(crossprod(sqrt(mat.KhR)%*%mat.XhSqR)))*length(vec.XhR)

    arr.GinvL[,,d] <- mat.GinvSqL
    arr.GinvR[,,d] <- mat.GinvSqR

    vec.ZL <- eval(parse(text = paste0(paste0("df_data$vec.Z.",d),"[df_data$vec.X<0]")))
    vec.ZR <- eval(parse(text = paste0(paste0("df_data$vec.Z.",d),"[df_data$vec.X>=0]")))

    mat.betaL[d,] <- list.resultRDRobustZBias$beta_p_l
    mat.betaR[d,] <- list.resultRDRobustZBias$beta_p_r

    # This residual should be the raw conditional mean projection errors of Z, not normalized by bandwidths
    vec.residualL <- vec.ZL - (mat.betaL[d,1] + mat.betaL[d,2]*vec.XL + mat.betaL[d,3]*vec.XL^2)
    vec.residualR <- vec.ZR - (mat.betaR[d,1] + mat.betaR[d,2]*vec.XR + mat.betaR[d,3]*vec.XR^2)

    # Initialize the containers of near-neighbor estimates of residuals
    vec.residualNnL <- vec.residualL
    vec.residualNnR <- vec.residualR
    for (i in 1:length(vec.XL))
    {
      #vec.residualNnL[i] <- vec.residualNnL[i] - (vec.residualL[vec.indexClose1stL[i]] + vec.residualL[vec.indexClose2ndL[i]] + vec.residualL[vec.indexClose3rdL[i]])/3
      for (j in 1:int_J)
      {
        vec.residualNnL[i] <- vec.residualNnL[i] - vec.residualL[mat.indexCloseL[i,j]]/int_J
      }
    }
    for (i in 1:length(vec.XR))
    {
      #vec.residualNnR[i] <- vec.residualNnR[i] - (vec.residualR[vec.indexClose1stR[i]] + vec.residualR[vec.indexClose2ndR[i]] + vec.residualR[vec.indexClose3rdR[i]])/3
      for (j in 1:int_J)
      {
        vec.residualNnR[i] <- vec.residualNnR[i] - vec.residualR[mat.indexCloseR[i,j]]/int_J
      }
    }

    arr.PsiL[,,d] <- diag(vec.residualNnL) %*% (mat.KhL %*% mat.XhSqL)
    arr.PsiR[,,d] <- diag(vec.residualNnR) %*% (mat.KhR %*% mat.XhSqR)

  }
  real.meanStatTZ <- mean(vec.statTZ)
  real.medianStatTZ <- median(vec.statTZ)
  real.maxAbsStatTZ <- max(abs(vec.statTZ))

  mat.CovL <- matrix(rep(0,int_dim_Z^2),int_dim_Z,int_dim_Z)
  mat.CovR <- matrix(rep(0,int_dim_Z^2),int_dim_Z,int_dim_Z)
  for (d1 in 1:int_dim_Z)
  {
    for (d2 in 1:int_dim_Z)
    {
      if (d2 >= d1)
      {
        mat.PsiL <- (int_J/(int_J+1))*t(arr.PsiL[,,d1]) %*% arr.PsiL[,,d2] / length(vec.XL)
        mat.PsiR <- (int_J/(int_J+1))*t(arr.PsiR[,,d1]) %*% arr.PsiR[,,d2] / length(vec.XR)

        mat.VarL <- arr.GinvL[,,d1] %*% mat.PsiL %*% arr.GinvL[,,d2]/length(vec.XL)
        mat.VarR <- arr.GinvR[,,d1] %*% mat.PsiR %*% arr.GinvR[,,d2]/length(vec.XR)

        # Multiply the hjk = (hj * hk)^(1/2)
        mat.CovL[d1,d2] <- mat.VarL[1,1] * sqrt(vec.hL[d1] * vec.hL[d2])
        mat.CovR[d1,d2] <- mat.VarR[1,1] * sqrt(vec.hR[d1] * vec.hR[d2])
      } else {
        mat.CovL[d1,d2] <- mat.CovL[d2,d1]
        mat.CovR[d1,d2] <- mat.CovR[d2,d1]
      }
    }
  }
  mat.Cov <- mat.CovL + mat.CovR
  # Z.stat is standardized vector: multiply the se back
  # also, multiply hk^(1/2): we take vec.hL here, but not that bandwidths are chosen to be the same for both sides.
  for (l in 1:int_dim_Z) {vec.statTZ[l] <- vec.statTZ[l]*vec.seZ[l]*(vec.hL[l]^(1/2))/sqrt(mat.Cov[l,l])}

  mat.Cor <- cov2cor(mat.Cov)

  if (bool_max_test_inv_V) {
    if (bool.knownV == TRUE) {
      mat.Cor <- diag(int_dim_Z)
      for (l1 in 1:int_dim_Z) { for (l2 in 1:int_dim_Z) {if (l1 != l2) {mat.Cor[l1,l2] = real.covZ} }}
    }
    U <- svd(mat.Cor)$u
    V <- svd(mat.Cor)$v
    D2inv <- sqrt(chol2inv(chol(diag(sqrt(svd(mat.Cor)$d)))))
    mat.CorSqrtInv <- U %*% D2inv %*% t(V)
    # take it to the identity matrix
    mat.Cor <- diag(int_dim_Z)
    real.statMax <- max((mat.CorSqrtInv %*% vec.statTZ)^2)
  } else {
    real.statMax <- max(vec.statTZ^2)
  }
  list.resultCovariatesStatMax <- list(real.statMax=real.statMax,
                                       mat.Cor = mat.Cor,
                                       bool.rejectBonferroniNull=bool.rejectBonferroniNull,
                                       bool.rejectNaiveNull=bool.rejectNaiveNull,
                                       real.meanStatTZ=real.meanStatTZ,
                                       real.medianStatTZ=real.medianStatTZ,
                                       real.maxAbsStatTZ=real.maxAbsStatTZ,
                                       real.effectiveNMeanZ=real.effectiveNMeanZ)
  return(list.resultCovariatesStatMax)
}
