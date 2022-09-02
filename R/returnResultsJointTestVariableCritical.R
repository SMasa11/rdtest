####
#' return joint test result
#'
#' @param df.data data.frame of data to evaluate
#' @param int.dimZ  integer dimension of covariates
#' @param int.J integer of nearest neighbor for the variance estimation <= 3
#' @param bool.maxTest boolean option to use max test rather than chi-squared test
#' @param bool.maxTestVinv boolean option for using inv V for max test
#' @param bool.knownV boolean option for setting known V instead of estimated Vinv
#' @param real.covZ real value of correlation in V, nullified if bool.knownV is FALSE or in default
#' @param fun.adjustBand real-valued function which returns a fraction of reduced bandwidth
#' @param real.criticalAlt real-valued critical value (in nominal coverage) to adjust the size of the joint test
#'
#' @export
####

returnResultsJointTestVariableCritical <- function(df.data,int.dimZ,int.J=3,bool.maxTest=FALSE,bool.maxTestVinv=FALSE,bool.knownV=FALSE,real.covZ=NULL,fun.adjustBand=NULL,real.criticalAlt=0.95)
{
  library("rdrobust")
  library("rddensity")
  library("R.utils")
  int.numSimulDraws <- 3000

  # Zstat
  if (bool.maxTest == FALSE) {
    list.resultCovariates <- computeStatChiSquared(df.data=df.data,int.dimZ=int.dimZ,bool.joint=TRUE,int.J=int.J,bool.knownV=bool.knownV,real.covZ=real.covZ,fun.adjustBand=fun.adjustBand)
  } else {
    list.resultCovariates <- computeStatMax(df.data=df.data,int.dimZ=int.dimZ,bool.joint=TRUE,int.J=int.J,bool.maxTestVinv=bool.maxTestVinv,bool.knownV=bool.knownV,real.covZ=real.covZ,fun.adjustBand=fun.adjustBand)
  }

  # Xstat
  list.resultX <- rddensity(X=df.data$vec.X)
  if (is.null(fun.adjustBand) == FALSE) {
    real.hLTemp <- list.resultX$h$left * fun.adjustBand(int.dimZ)
    real.hRTemp <- list.resultX$h$right * fun.adjustBand(int.dimZ)
    vec.hAdjusted <- c(real.hLTemp,real.hRTemp)
    list.resultX <- rddensity(X=df.data$vec.X,h=vec.hAdjusted)
  }
  # t stat
  real.statTX <- list.resultX$test$t_jk

  bool.rejectNaiveNull <- list.resultCovariates$bool.rejectNaiveNull
  if (abs(real.statTX) > qnorm(1- 0.025)) {
    bool.rejectNaiveNull <- TRUE
  }

  bool.rejectBonferroniNull <- list.resultCovariates$bool.rejectBonferroniNull
  if (abs(real.statTX) > qnorm(1 - 0.025/(int.dimZ+1))) {
    bool.rejectBonferroniNull <- TRUE
  }

  if (bool.maxTest==FALSE) {
    bool.rejectJointNull <- (list.resultCovariates$real.statChiSq + real.statTX^2 >= qchisq(real.criticalAlt, df=(int.dimZ+1)))
  } else {
    if (bool.maxTestVinv == TRUE) {
      # standardized max stats
      real.statMaxJoint <- max(list.resultCovariates$real.statMax,real.statTX^2)
      bool.rejectJointNull <- (real.statMaxJoint >= qchisq(real.criticalAlt^(1/(int.dimZ+1)),1))
      real.criticalValueJoint <- qchisq(real.criticalAlt^(1/(int.dimZ+1)),1)
    } else {
      real.statMaxJoint <- max(list.resultCovariates$real.statMax,real.statTX^2)
      vec.randomDrawStats <- rep(0,int.numSimulDraws)
      #Generates the same random normal vector with the seed of 373568, restoring the currently loading seed
      mat.randomNorm <- matrix(withSeed({rnorm((int.dimZ+1)*int.numSimulDraws,0,1)},seed=373568),int.dimZ+1,int.numSimulDraws)
      mat.Cor <- matrix(rep(0,(int.dimZ+1)^2),(int.dimZ+1),(int.dimZ+1))
      mat.Cor[1:int.dimZ,1:int.dimZ] <- list.resultCovariates$mat.Cor
      mat.Cor[(int.dimZ+1),(int.dimZ+1)] <- 1
      U <- svd(mat.Cor)$u
      V <- svd(mat.Cor)$v
      D <- diag(sqrt(svd(mat.Cor)$d))
      mat.CorSqrt <- U %*% D %*% t(V)
      for (l in seq(1:int.numSimulDraws)) {
        # Replicating the (signed) vector of statistics by multiplying the corr matrix, then take square and max.
        vec.randomDrawStats[l] <- max((mat.CorSqrt %*% mat.randomNorm[,l])^2)
      }
      real.pvalue <- mean((vec.randomDrawStats >= real.statMaxJoint))
      bool.rejectJointNull <- (real.pvalue <= (1 - real.criticalAlt))
      real.criticalValueJoint <- quantile(vec.randomDrawStats,real.criticalAlt)
    }
  }


  real.covZ1Z2 <- NULL
  if (int.dimZ > 1) {real.covZ1Z2 <- cor(df.data$vec.Z.1[df.data$vec.X <= quantile(abs(df.data$vec.X),0.1)],df.data$vec.Z.2[df.data$vec.X <= quantile(abs(df.data$vec.X),0.1)])} #cor(df.data$vec.Z.1,df.data$vec.Z.2)}
  real.effNx <- list.resultX$N$eff_left + list.resultX$N$eff_right

  if (bool.maxTest == FALSE) {
    list.resultJointTest <- list(real.statChiSq = list.resultCovariates$real.statChiSq,
                                 real.statTX = real.statTX,
                                 real.meanStatTZ = list.resultCovariates$real.meanStatTZ,
                                 real.medianStatTZ = list.resultCovariates$real.medianStatTZ,
                                 real.maxAbsStatTZ = list.resultCovariates$real.maxAbsStatTZ,
                                 real.effectiveNMeanZ = list.resultCovariates$real.effectiveNMeanZ,
                                 real.effNx = real.effNx,
                                 real.covZ1Z2 = real.covZ1Z2,
                                 bool.rejectJointNull = bool.rejectJointNull,
                                 bool.rejectNaiveNull = bool.rejectNaiveNull,
                                 bool.rejectBonferroniNull = bool.rejectBonferroniNull,
                                 real.criticalValueJoint = NA)
  } else {
    list.resultJointTest <- list(real.statMaxJoint = real.statMaxJoint,
                                 real.statTX = real.statTX,
                                 real.meanStatTZ = list.resultCovariates$real.meanStatTZ,
                                 real.medianStatTZ = list.resultCovariates$real.medianStatTZ,
                                 real.maxAbsStatTZ = list.resultCovariates$real.maxAbsStatTZ,
                                 real.effectiveNMeanZ = list.resultCovariates$real.effectiveNMeanZ,
                                 real.effNx = real.effNx,
                                 real.covZ1Z2 = real.covZ1Z2,
                                 bool.rejectJointNull = bool.rejectJointNull,
                                 bool.rejectNaiveNull = bool.rejectNaiveNull,
                                 bool.rejectBonferroniNull = bool.rejectBonferroniNull,
                                 real.criticalValueJoint = real.criticalValueJoint)
  }

  return(list.resultJointTest)
}


