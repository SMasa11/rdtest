####
#' MC joint test of coverage / power check using simulated data
#'
#' Example code for simulation excercise
#'
#' @param int_ns Interger, number of simulations to run
#'  default is 300.
#' @param jump Real, z jump size.
#' @param x_jump Real, x jump size.
#' @param a_2 Real, gap in heteroskedasticity.
#' @param dim Integer, number of z variables.
#' @param n Integer, sample size.
#' @param cov Real, correlation of Z variables.
#' @param frac_jump String, either "None", "Half", "All", "Last".
#' @param bool_mutePrint Boolean, option for muting outputs.
#'  default is FALSE
#' @param int_J Integer, nearest neighbor for the variance estimation <= 3
#'  default is 3.
#' @param bool_maxTest Boolean, use max test instead of L2 test,
#'  default is false.
#' @param bool_maxTestVinv Boolean, option for using inverse V for max test.
#' @param bool_knownV Boolean, option for setting V as known.
#' @param fun.adjustBand Function, returns a fraction of reduced bandwidth.
#' @param bool_chisqStd Boolean, for using standardized t stat for L2 test.
#'
####

returnResultsMCJointTest <- function(int_ns = 300,
                                     x_jump,
                                     jump,
                                     a_2,
                                     dim,
                                     n,
                                     cov_z,
                                     frac_jump = "Half",
                                     bool_mutePrint = FALSE,
                                     int_J = 3,
                                     bool_maxTest = FALSE,
                                     bool_maxTestVinv = FALSE,
                                     bool_knownV = FALSE,
                                     fun_adjustBand = NULL,
                                     bool_chisqStd = FALSE)
{
  set.seed(52622)

  if (jump != 0)
  {
    if (frac_jump == 0)
    {
      warning(
        "WARNING:
        JUMP IS SPECIFIED BUT NOT APPLIED AS frac.jump IS NOT SPECIFIED.")
    }
  }

  ns <- int_ns
  naive.num.reject <- 0
  bonfe.num.reject <- 0
  joint.num.reject <- 0
  chi.stat.vec <- rep(0,length(ns))
  vec.statMaxJoint <- rep(0,length(ns))
  effN.z.vec <- c(1:ns)
  effN.x.vec <- c(1:ns)
  covZ1Z2.vec <- c(1:ns)
  mean.jump.z.vec <- c(1:ns)
  median.jump.z.vec <- c(1:ns)
  max.jump.z.vec <- c(1:ns)
  vec.criticalValue <- c(1:ns)

  #fun.mu <- function(x){return(x^2)}
  fun_mu <- function(x){
    return(
      0.48
      + (1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5)*(x < 0)
      + (0.84*x - 3.00*x^2 + 7.99*x^3 - 9.01*x^4 + 3.56*x^5)*(x >= 0))}

  options <- list(int.DGP = 3,
                  int.n = n,
                  real.jumpX = x_jump,
                  real.jumpZ = jump,
                  real.a2Param = a_2,
                  int.dimZ = dim,
                  real.covZ = cov_z,
                  str.fracJumpZ = frac_jump,
                  fun.mu = fun_mu)

  for (loop in 1:ns)
  {
    data <- simulateDGPs(options)
    list.resultJointTest <- return_result_joint(
      df_data = data,
      int_dimZ = options$int.dimZ,
      int_J = int_J,
      bool_maxTest = bool_maxTest,
      bool_maxTestVinv = bool_maxTestVinv,
      bool_knownV = bool_knownV,
      real_covZ = cov_z,
      fun_adjustBand = fun_adjustBand,
      bool_chisqStd = bool_chisqStd)

    naive.num.reject <-
      naive.num.reject + list.resultJointTest$bool.rejectNaiveNull
    bonfe.num.reject <-
      bonfe.num.reject + list.resultJointTest$bool.rejectBonferroniNull
    joint.num.reject <-
      joint.num.reject + list.resultJointTest$bool.rejectJointNull

    if (bool_maxTest == FALSE) {
      chi.stat.vec[loop] <- list.resultJointTest$real.statChiSq
    } else {
      vec.statMaxJoint[loop]  <- list.resultJointTest$real.statMaxJoint
    }
    mean.jump.z.vec[loop] <- list.resultJointTest$real.meanStatTZ
    median.jump.z.vec[loop] <- list.resultJointTest$real.medianStatTZ
    max.jump.z.vec[loop] <- list.resultJointTest$real.maxAbsStatTZ
    effN.z.vec[loop] <- list.resultJointTest$real.effectiveNMeanZ
    effN.x.vec[loop] <- list.resultJointTest$real.effNx
    if (dim > 1) {
      covZ1Z2.vec[loop] <- list.resultJointTest$real.covZ1Z2
    }

    vec.criticalValue[loop] <- list.resultJointTest$real.criticalValueJoint

    if (loop %% 100 == 0)
    {
      if (!bool_mutePrint) {message(loop)}
    }
  }
  if (!bool_mutePrint) {
    print("average effective sample size")
    print(mean(effN.z.vec))

    print("average correlation of Z1 and Z2")
    print(mean(covZ1Z2.vec))

    print("mean z stat")
    print(mean(mean.jump.z.vec))

    print("median z stat")
    print(mean(median.jump.z.vec))

    print("max z stat")
    print(mean(max.jump.z.vec))


    print("naive rejection rate")
    print(naive.num.reject/ns)
    print("bonferroni rejection rate")
    print(bonfe.num.reject/ns)
    print("joint rejection rate")
    if (bool_maxTest == FALSE) {joint.num.reject <- joint.num.reject[1,1]}
    print(joint.num.reject/ns)
  }

  summary <- list(vec.statMaxJoint = vec.statMaxJoint,
                  chi.stat.vec = chi.stat.vec,
                  effN.mean = mean(effN.z.vec),
                  effNx.mean = mean(effN.x.vec),
                  cor.z = mean(covZ1Z2.vec),
                  naive = naive.num.reject/ns,
                  bonfe = bonfe.num.reject/ns,
                  joint = joint.num.reject/ns,
                  maxCritical = vec.criticalValue)
  return(summary)
}


