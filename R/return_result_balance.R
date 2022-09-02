####
#' return covariates balance test results
#'
#'
#' @param df.data data.frame of data to evaluate
#' @param int.dimZ  integer dimension of covariates
#' @param int.J integer of nearest neighbor for the variance estimation <= 3
#' @param bool.maxTest boolean option to use max test rather than chi-squared test
#'
####

return_result_balance <- function(df_data,
                                  int_dimZ,
                                  int_J = 3,
                                  bool_maxTest = FALSE)
{
  if (int.dimZ < 2) {stop("For the balance test with a single covariate, use rdrobust package.")}

  int.numSimulDraws <- 3000

  real.covZ1Z2 <- cor(df.data$vec.Z.1,df.data$vec.Z.2)

  if (bool.maxTest == FALSE) {
    list.resultCovariatesPart <- compute_stat_L2(df.data=df.data,int.dimZ=int.dimZ,bool.joint=FALSE,int.J=int.J)
  } else {
    list.resultCovariatesPart <- compute_stat_max(df.data=df.data,int.dimZ=int.dimZ,bool.joint=FALSE,int.J=int.J)
  }

  if (bool.maxTest==FALSE) {
    bool.rejectJointNull <- (list.resultCovariatesPart$real.statChiSq >= qchisq(.95, df=int.dimZ))
  } else {
    real.statMaxJoint <- list.resultCovariatesPart$real.statMax
    vec.randomDrawStats <- rep(0,int.numSimulDraws)
    #Generates the same random normal vector with the seed of 373568, restoring the currently loading seed
    mat.randomNorm <- matrix(R.utils::withSeed({rnorm((int.dimZ+1)*int.numSimulDraws,0,1)},seed=373568),int.dimZ,int.numSimulDraws)
    mat.Cor <- matrix(rep(0,(int.dimZ)^2),(int.dimZ),(int.dimZ))
    mat.Cor[1:int.dimZ,1:int.dimZ] <- list.resultCovariatesPart$mat.Cor
    U <- svd(mat.Cor)$u
    V <- svd(mat.Cor)$v
    D <- diag(sqrt(svd(mat.Cor)$d))
    mat.CorSqrt <- U %*% D %*% t(V)
    for (l in seq(1:int.numSimulDraws)) {
      # Replicating the (signed) vector of statistics by multiplying the corr matrix, then take square and max.
      vec.randomDrawStats[l] <- max((mat.CorSqrt %*% mat.randomNorm[,l])^2)
    }
    real.pvalue <- mean((vec.randomDrawStats >= real.statMaxJoint))
    bool.rejectJointNull <- (real.pvalue <= 0.05)
    real.criticalValueJoint <- quantile(vec.randomDrawStats,0.95)
  }

  if (bool.maxTest == FALSE) {
    list.resultCovariateBalance <- list(bool.rejectNaiveNull=list.resultCovariatesPart$bool.rejectNaiveNull,
                                      bool.rejectBonferroniNull=list.resultCovariatesPart$bool.rejectBonferroniNull,
                                      bool.rejectJointNull=bool.rejectJointNull,
                                      real.statChiSq=list.resultCovariatesPart$real.statChiSq,
                                      real.meanStatTZ=list.resultCovariatesPart$real.meanStatTZ,
                                      real.medianStatTZ=list.resultCovariatesPart$real.medianStatTZ,
                                      real.maxAbsStatTZ=list.resultCovariatesPart$real.maxAbsStatTZ,
                                      real.covZ1Z2=real.covZ1Z2,
                                      real.effectiveNMeanZ=list.resultCovariatesPart$real.effectiveNMeanZ)
  } else {
    list.resultCovariateBalance <- list(real.statMaxJoint = real.statMaxJoint,
                               real.meanStatTZ = list.resultCovariatesPart$real.meanStatTZ,
                               real.medianStatTZ = list.resultCovariatesPart$real.medianStatTZ,
                               real.maxAbsStatTZ = list.resultCovariatesPart$real.maxAbsStatTZ,
                               real.effectiveNMeanZ = list.resultCovariatesPart$real.effectiveNMeanZ,
                               real.covZ1Z2 = real.covZ1Z2,
                               bool.rejectJointNull = bool.rejectJointNull,
                               bool.rejectNaiveNull = list.resultCovariatesPart$bool.rejectNaiveNull,
                               bool.rejectBonferroniNull = list.resultCovariatesPart$bool.rejectBonferroniNull,
                               real.criticalValueJoint = real.criticalValueJoint)
  }
  return(list.resultCovariateBalance)
}
