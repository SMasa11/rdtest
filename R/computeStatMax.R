####
#' return max statistics for covariates balance
#'
#' @param df.data data.frame of data to evaluate
#' @param int.dimZ  integer dimension of covariates
#' @param bool.joint Boolean TRUE if joint test, FALSE if covariate balance test only
#' @param int.J integer of nearest neighbor for the variance estimation <= 3
#' @param bool.maxTestVinv boolean option for using inv V for max test
#' @param bool.knownV boolean option for setting known V instead of estimated Vinv
#' @param real.covZ real value of correlation in V, nullified if bool.knownV is FALSE or in default
#' @param fun.adjustBand real-valued function which returns a fraction of reduced bandwidth
#'
#' @export
####

computeStatMax <- function(df.data,int.dimZ,bool.joint,int.J,bool.maxTestVinv,bool.knownV=FALSE,real.covZ=NULL,fun.adjustBand=NULL)
{
  # booleans which return TRUE if tests are rejected by covariates only
  bool.rejectNaiveNull <- FALSE
  bool.rejectBonferroniNull <- FALSE

  real.effectiveNMeanZ <- 0

  vec.hL <- rep(0,int.dimZ)
  vec.hR <- rep(0,int.dimZ)
  vec.VRobustL <- rep(0,int.dimZ)
  vec.VRobustR <- rep(0,int.dimZ)

  mat.betaL <- matrix(rep(0,int.dimZ*3),int.dimZ,3)
  mat.betaR <- matrix(rep(0,int.dimZ*3),int.dimZ,3)

  vec.XL <- df.data$vec.X[df.data$vec.X<0,drop=FALSE]
  vec.XR <- df.data$vec.X[df.data$vec.X>=0,drop=FALSE]
  mat.difInXL <- matrix(rep(t(vec.XL),length(vec.XL)),length(vec.XL),length(vec.XL))
  mat.difInXL <- abs(mat.difInXL - t(mat.difInXL)) + diag(Inf,length(vec.XL),length(vec.XL))

  mat.indexCloseL <- matrix(rep(0,length(vec.XL)*int.J),length(vec.XL),int.J)
  for (j in 1:int.J) {
    mat.indexCloseL[,j] <- max.col(-mat.difInXL)
    for (i in 1:length(vec.XL)) {mat.difInXL[i,mat.indexCloseL[i,j]] <- Inf}
  }
  # vec.indexClose1stL <- max.col(-mat.difInXL)
  # for (i in 1:length(vec.XL)) {mat.difInXL[i,vec.indexClose1stL[i]] <- Inf}
  # vec.indexClose2ndL <- max.col(-mat.difInXL)
  # for (i in 1:length(vec.XL)) {mat.difInXL[i,vec.indexClose2ndL[i]] <- Inf}
  # vec.indexClose3rdL <- max.col(-mat.difInXL)

  mat.difInXR <- matrix(rep(t(vec.XR),length(vec.XR)),length(vec.XR),length(vec.XR))
  mat.difInXR <- abs(mat.difInXR - t(mat.difInXR)) + diag(Inf,length(vec.XR),length(vec.XR))

  mat.indexCloseR <- matrix(rep(0,length(vec.XR)*int.J),length(vec.XR),int.J)
  for (j in 1:int.J) {
    mat.indexCloseR[,j] <- max.col(-mat.difInXR)
    for (i in 1:length(vec.XR)) {mat.difInXR[i,mat.indexCloseR[i,j]] <- Inf}
  }
  # vec.indexClose1stR <- max.col(-mat.difInXR)
  # for (i in 1:length(vec.XR)) {mat.difInXR[i,vec.indexClose1stR[i]] <- Inf}
  # vec.indexClose2ndR <- max.col(-mat.difInXR)
  # for (i in 1:length(vec.XR)) {mat.difInXR[i,vec.indexClose2ndR[i]] <- Inf}
  # vec.indexClose3rdR <- max.col(-mat.difInXR)

  # Arrays to store
  mat.XhL <- matrix(rep(0,length(vec.XL)*int.dimZ),length(vec.XL),int.dimZ)
  mat.XhR <- matrix(rep(0,length(vec.XR)*int.dimZ),length(vec.XR),int.dimZ)

  arr.PsiL <- array(rep(0,3*length(vec.XL)*int.dimZ),dim=c(length(vec.XL),3,int.dimZ))
  arr.PsiR <- array(rep(0,3*length(vec.XR)*int.dimZ),dim=c(length(vec.XR),3,int.dimZ))

  arr.GinvL <- array(rep(0,3^2*int.dimZ),dim=c(3,3,int.dimZ))
  arr.GinvR <- array(rep(0,3^2*int.dimZ),dim=c(3,3,int.dimZ))

  mat.residualNnL <- matrix(rep(0,length(vec.XL)*int.dimZ),length(vec.XL),int.dimZ)
  mat.residualNnR <- matrix(rep(0,length(vec.XR)*int.dimZ),length(vec.XR),int.dimZ)
  ###

  vec.statTZ <- rep(0,int.dimZ)
  vec.seZ <- rep(0,int.dimZ)


  for (d in 1:int.dimZ)
  {
    # enforce the bandwidths are the same for both sides
    # Note: density bands choose different bands by default, but the value seems to be sensitive to the option. I'll keep them as their default.
    eval(parse(text = paste0(paste0("list.resultRDRobustZ <- rdrobust(y=df.data$vec.Z.",d),",x=df.data$vec.X,rho=1,bwselect='mserd')")))

    # Zstat
    real.statTZ <- list.resultRDRobustZ$z[3]
    vec.statTZ[d] <- real.statTZ
    vec.seZ[d] <- list.resultRDRobustZ$se[3]
    vec.hL[d] <- list.resultRDRobustZ$bws[1,1]
    vec.hR[d] <- list.resultRDRobustZ$bws[1,2]
    if (is.null(fun.adjustBand) == FALSE) {
      vec.hL[d] <- vec.hL[d]*fun.adjustBand(int.dimZ)
      vec.hR[d] <- vec.hR[d]*fun.adjustBand(int.dimZ)
      eval(parse(text = paste0(paste0("list.resultRDRobustZ <- rdrobust(y=df.data$vec.Z.",d),",x=df.data$vec.X,rho=1,h=c(vec.hL[d],vec.hR[d]))")))
    }
    real.effectiveNMeanZ <- real.effectiveNMeanZ + sum(list.resultRDRobustZ$N_h)/int.dimZ

    eval(parse(text = paste0(paste0("list.resultRDRobustZBias <- rdrobust(y=df.data$vec.Z.",d),",x=df.data$vec.X,p=2,rho=1,h=c(vec.hL[d],vec.hR[d]))")))

    vec.VRobustL[d] <- list.resultRDRobustZ$V_rb_l[1,1]
    vec.VRobustR[d] <- list.resultRDRobustZ$V_rb_r[1,1]
    if (abs(real.statTZ) > qnorm(1- 0.025)) {
      bool.rejectNaiveNull <- TRUE
    }
    if (abs(real.statTZ) > qnorm(1 - 0.025/(int.dimZ+bool.joint))){
      bool.rejectBonferroniNull <- TRUE
    }

    # For each d, compute X(h_d)
    vec.XhL <- df.data$vec.X[df.data$vec.X< 0,drop=FALSE]/vec.hL[d]
    vec.XhR <- df.data$vec.X[df.data$vec.X>=0,drop=FALSE]/vec.hR[d]

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

    vec.ZL <- eval(parse(text = paste0(paste0("df.data$vec.Z.",d),"[df.data$vec.X<0]")))
    vec.ZR <- eval(parse(text = paste0(paste0("df.data$vec.Z.",d),"[df.data$vec.X>=0]")))

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
      for (j in 1:int.J)
      {
        vec.residualNnL[i] <- vec.residualNnL[i] - vec.residualL[mat.indexCloseL[i,j]]/int.J
      }
    }
    for (i in 1:length(vec.XR))
    {
      #vec.residualNnR[i] <- vec.residualNnR[i] - (vec.residualR[vec.indexClose1stR[i]] + vec.residualR[vec.indexClose2ndR[i]] + vec.residualR[vec.indexClose3rdR[i]])/3
      for (j in 1:int.J)
      {
        vec.residualNnR[i] <- vec.residualNnR[i] - vec.residualR[mat.indexCloseR[i,j]]/int.J
      }
    }

    arr.PsiL[,,d] <- diag(vec.residualNnL) %*% (mat.KhL %*% mat.XhSqL)
    arr.PsiR[,,d] <- diag(vec.residualNnR) %*% (mat.KhR %*% mat.XhSqR)

  }
  real.meanStatTZ <- mean(vec.statTZ)
  real.medianStatTZ <- median(vec.statTZ)
  real.maxAbsStatTZ <- max(abs(vec.statTZ))

  mat.CovL <- matrix(rep(0,int.dimZ^2),int.dimZ,int.dimZ)
  mat.CovR <- matrix(rep(0,int.dimZ^2),int.dimZ,int.dimZ)
  for (d1 in 1:int.dimZ)
  {
    for (d2 in 1:int.dimZ)
    {
      if (d2 >= d1)
      {
        mat.PsiL <- (int.J/(int.J+1))*t(arr.PsiL[,,d1]) %*% arr.PsiL[,,d2] / length(vec.XL)
        mat.PsiR <- (int.J/(int.J+1))*t(arr.PsiR[,,d1]) %*% arr.PsiR[,,d2] / length(vec.XR)

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
  for (l in 1:int.dimZ) {vec.statTZ[l] <- vec.statTZ[l]*vec.seZ[l]*(vec.hL[l]^(1/2))/sqrt(mat.Cov[l,l])}

  mat.Cor <- cov2cor(mat.Cov)

  if (bool.maxTestVinv) {
    if (bool.knownV == TRUE) {
      mat.Cor <- diag(int.dimZ)
      for (l1 in 1:int.dimZ) { for (l2 in 1:int.dimZ) {if (l1 != l2) {mat.Cor[l1,l2] = real.covZ} }}
    }
    U <- svd(mat.Cor)$u
    V <- svd(mat.Cor)$v
    D2inv <- sqrt(chol2inv(chol(diag(sqrt(svd(mat.Cor)$d)))))
    mat.CorSqrtInv <- U %*% D2inv %*% t(V)
    # take it to the identity matrix
    mat.Cor <- diag(int.dimZ)
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
