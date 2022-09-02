####
#' return joint test result
#'
#' @param df_data data.frame of data to evaluate
#' @param int_dimZ  integer dimension of covariates
#' @param int_J integer of nearest neighbor for the variance estimation <= 3
#' @param bool_maxTest boolean option to use max test rather than chi-squared test
#' @param bool_maxTestVinv boolean option for using inv V for max test
#' @param bool_knownV boolean option for setting known V instead of estimated Vinv
#' @param real_covZ real value of correlation in V, nullified if bool.knownV is FALSE or in default
#' @param fun_adjustBand real-valued function which returns a fraction of reduced bandwidth
#' @param bool_chisqStd boolean option for using standardized t stat to conduct L2 test
#'
####

return_result_joint <- function(df_data,
                                   int_dimZ,
                                   int_J = 3,
                                   bool_maxTest = FALSE,
                                   bool_maxTestVinv = FALSE,
                                   bool_knownV = FALSE,
                                   real_covZ = NULL,
                                   fun_adjustBand = NULL,
                                   bool_chisqStd = FALSE)
{
  int_numSimulDraws <- 3000

  # Zstat
  if (bool_maxTest == FALSE) {
    list_resultCovariates <-
      compute_stat_L2(df.data = df_data,
                            int.dimZ = int_dimZ,
                            bool.joint = TRUE,
                            int.J = int_J,
                            bool.knownV = bool_knownV,
                            real.covZ = real_covZ,
                            fun.adjustBand = fun_adjustBand,
                            bool.chisqStd = bool_chisqStd)
  } else {
    list_resultCovariates <-
      compute_stat_max(df.data = df_data,
                     int.dimZ = int_dimZ,
                     bool.joint = TRUE,
                     int.J = int_J,
                     bool.maxTestVinv = bool_maxTestVinv,
                     bool.knownV = bool_knownV,
                     real.covZ = real_covZ,
                     fun.adjustBand = fun_adjustBand)
  }

  # Xstat
  list_resultX <- rddensity::rddensity(X=df_data$vec.X)
  if (is.null(fun_adjustBand) == FALSE) {
    real_hLTemp <- list_resultX$h$left * fun_adjustBand(int_dimZ)
    real_hRTemp <- list_resultX$h$right * fun_adjustBand(int_dimZ)
    vec_hAdjusted <- c(real_hLTemp,real_hRTemp)
    list_resultX <- rddensity::rddensity(X=df_data$vec_X,
                                         h=vec_hAdjusted)
  }
  # t stat
  real_statTX <- list_resultX$test$t_jk

  bool_rejectNaiveNull <- list_resultCovariates$bool.rejectNaiveNull
  if (abs(real_statTX) > qnorm(1- 0.025)) {
    bool_rejectNaiveNull <- TRUE
  }

  bool_rejectBonferroniNull <- list_resultCovariates$bool.rejectBonferroniNull
  if (abs(real_statTX) > qnorm(1 - 0.025/(int_dimZ+1))) {
    bool_rejectBonferroniNull <- TRUE
  }

  if (bool_maxTest==FALSE) {
    if (bool_chisqStd) {
      real_statChisqStdJoint <-
        list_resultCovariates$real.statChiSq + real_statTX^2
      vec_randomDrawStats <- rep(0,int_numSimulDraws)
      mat_Cor <- matrix(rep(0,(int_dimZ+1)^2),(int_dimZ+1),(int_dimZ+1))
      mat_Cor[1:int_dimZ,1:int_dimZ] <- list_resultCovariates$mat.Cor
      mat_Cor[(int_dimZ+1),(int_dimZ+1)] <- 1
      mat_randomDrawStats <-
        return_random_stat(int_num_simul = int_numSimulDraws,
                           int_dim = int_dimZ+1,
                           mat_Cor = mat_Cor)
      real_pvalue <- 0
      for (l in seq(1:int_numSimulDraws)) {
        # Replicating the (signed) vector of statistics
        # by multiplying the corr matrix, then take square and max.
        vec_randomDrawStats[l] <- sum(mat_randomDrawStats[,l]^2)
        real_pvalue <- real_pvalue +
          (vec_randomDrawStats[l] >= real_statChisqStdJoint)/int_numSimulDraws
      }
      bool_rejectJointNull <- (real_pvalue <= 0.05)
      real_criticalValueJoint <- quantile(vec_randomDrawStats,0.95)
    } else {
      bool_rejectJointNull <- (list_resultCovariates$real.statChiSq +
                                 real_statTX^2 >= qchisq(.95, df=(int_dimZ+1)))
    }

  } else {
    if (bool_maxTestVinv == TRUE) {
      # standardized max stats
      real_statMaxJoint <-
        max(list_resultCovariates$real_statMax,real_statTX^2)
      bool_rejectJointNull <-
        (real_statMaxJoint >= qchisq(0.95^(1/(int_dimZ+1)),1))
      real_criticalValueJoint <- qchisq(0.95^(1/(int_dimZ+1)),1)
    } else {
      real_statMaxJoint <- max(list_resultCovariates$real.statMax,real_statTX^2)
      vec_randomDrawStats <- rep(0,int_numSimulDraws)
      mat_Cor <- matrix(rep(0,(int_dimZ+1)^2),(int_dimZ+1),(int_dimZ+1))
      mat_Cor[1:int_dimZ,1:int_dimZ] <- list_resultCovariates$mat.Cor
      mat_Cor[(int_dimZ+1),(int_dimZ+1)] <- 1
      mat_randomDrawStats <-
        return_random_stat(int_num_simul = int_numSimulDraws,
                           int_dim = int_dimZ+1,
                           mat_Cor = mat_Cor)
      for (l in seq(1:int_numSimulDraws)) {
        # Replicating the (signed) vector of statistics
        # by multiplying the corr matrix, then take square and max.
        vec_randomDrawStats[l] <- max(mat_randomDrawStats[,l]^2)
      }
      real_pvalue <- mean((vec_randomDrawStats >= real_statMaxJoint))
      bool_rejectJointNull <- (real_pvalue <= 0.05)
      real_criticalValueJoint <- quantile(vec_randomDrawStats,0.95)
    }
  }


  real_covZ1Z2 <- NULL
  if (int_dimZ > 1) {
    real_covZ1Z2 <-
      cor(df_data$vec.Z.1[df_data$vec.X <= quantile(abs(df_data$vec.X),0.1)],
          df_data$vec.Z.2[df_data$vec.X <= quantile(abs(df_data$vec.X),0.1)])}
  real_effNx <- list_resultX$N$eff_left + list_resultX$N$eff_right

  if (bool_maxTest == FALSE) {
    if (bool_chisqStd == FALSE) {
      real_criticalValueJoint = NA
    }
    list_resultJointTest <-
      list(real.statChiSq = list_resultCovariates$real.statChiSq,
           real.statTX = real_statTX,
           real.meanStatTZ = list_resultCovariates$real.meanStatTZ,
           real.medianStatTZ = list_resultCovariates$real.medianStatTZ,
           real.maxAbsStatTZ = list_resultCovariates$real.maxAbsStatTZ,
           real.effectiveNMeanZ = list_resultCovariates$real.effectiveNMeanZ,
           real.effNx = real_effNx,
           real.covZ1Z2 = real_covZ1Z2,
           bool.rejectJointNull = bool_rejectJointNull,
           bool.rejectNaiveNull = bool_rejectNaiveNull,
           bool.rejectBonferroniNull = bool_rejectBonferroniNull,
           real.criticalValueJoint = real_criticalValueJoint)
  } else {
    list_resultJointTest <-
      list(real.statMaxJoint = real_statMaxJoint,
           real.statTX = real_statTX,
           real.meanStatTZ = list_resultCovariates$real.meanStatTZ,
           real.medianStatTZ = list_resultCovariates$real.medianStatTZ,
           real.maxAbsStatTZ = list_resultCovariates$real.maxAbsStatTZ,
           real.effectiveNMeanZ = list_resultCovariates$real.effectiveNMeanZ,
           real.effNx = real_effNx,
           real.covZ1Z2 = real_covZ1Z2,
           bool.rejectJointNull = bool_rejectJointNull,
           bool.rejectNaiveNull = bool_rejectNaiveNull,
           bool.rejectBonferroniNull = bool_rejectBonferroniNull,
           real.criticalValueJoint = real_criticalValueJoint)
  }

  return(list_resultJointTest)
}


