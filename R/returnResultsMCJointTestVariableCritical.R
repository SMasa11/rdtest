####
#' MC joint test of coverage / power check using simulated data
#'
#' Example code for simulation excercise
#'
#' @param int.ns number of simulations to run
#' @param jump z jump size
#' @param x.jump x jump size
#' @param a.2 gap in heteroskedasticity
#' @param dim number of z variables
#' @param n sample size
#' @param cov correlation of Z variables
#' @param frac.jump Either None, Half, All, Last
#' @param bool.mutePrint boolean option for muting outputs
#' @param int.J integer of nearest neighbor for the variance estimation <= 3
#' @param bool.maxTest boolean option for using the max test instead of chisqaured test, default is false
#' @param bool.maxTestVinv boolean option for using inverse V for max test
#' @param bool.knownV boolean option for setting V as known (instead of inverse V)
#' @param fun.adjustBand real-valued function which returns a fraction of reduced bandwidth
#' @param real.criticalAlt real-valued critical value (in nominal coverage) to adjust the size of the joint test
#'
#' @export
####

returnResultsMCJointTestVariableCritical <- function(int.ns = 300, x.jump,jump,a.2,dim,n,cov.z,frac.jump = "Half",bool.mutePrint=FALSE,int.J = 3,bool.maxTest=FALSE,bool.maxTestVinv=FALSE,bool.knownV=FALSE,fun.adjustBand=NULL,real.criticalAlt=0.95)
{
  library("rdrobust")
  library("rddensity")
  library("ggplot2")
  library("grid")
  set.seed(52622)

  if (jump != 0)
  {
    if (frac.jump == 0)
    {
      warning("WARNING: JUMP IS SPECIFIED BUT NOT APPLIED AS frac.jump IS NOT SPECIFIED.")
    }
  }

  ns <- int.ns
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
  fun.mu <- function(x){return(0.48
                                 + (1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5)*(x < 0)
                                 + (0.84*x - 3.00*x^2 + 7.99*x^3 - 9.01*x^4 + 3.56*x^5)*(x >= 0))}

  options <- list(int.DGP = 3, int.n = n,real.jumpX = x.jump,real.jumpZ = jump, real.a2Param = a.2,int.dimZ = dim,real.covZ = cov.z,str.fracJumpZ=frac.jump,fun.mu=fun.mu)

  for (loop in 1:ns)
  {
    data <- RDTests::simulateDGPs(options)
    list.resultJointTest <- returnResultsJointTestVariableCritical(df.data=data,int.dimZ=options$int.dimZ,int.J=int.J,bool.maxTest=bool.maxTest,bool.maxTestVinv=bool.maxTestVinv,bool.knownV=bool.knownV,real.covZ=cov.z,fun.adjustBand=fun.adjustBand,real.criticalAlt=real.criticalAlt)
    #
    # list.resultCovariatesPart <- computeStatChiSquared(df.data=data,int.dimZ=options$int.dimZ,bool.joint=TRUE)
    #
    # # Xstat
    # Xresult <- rddensity(X=data$vec.X)
    # effN.x.vec[loop] <- Xresult$N$eff_left + Xresult$N$eff_right
    # # t stat
    # stat.X <- Xresult$test$t_jk
    #
    # naive.null.reject <- list.resultCovariatesPart$bool.rejectNaiveNull
    # if (abs(stat.X) > qnorm(1- 0.025)) {
    #   naive.null.reject <- TRUE
    # }
    #
    # bonfe.null.reject <- list.resultCovariatesPart$bool.rejectBonferroniNull
    # if (abs(stat.X) > qnorm(1 - 0.025/(options$int.dimZ+1))) {
    #   bonfe.null.reject <- TRUE
    # }

    #
    #     joint.null.reject <- ((list.resultCovariatesPart$real.statChiSq >= qchisq(.95 + 0.05/(options$int.dimZ+1), df=(options$int.dimZ))) | (stat.X >= qnorm(p=(1 - 0.05/(options$int.dimZ+1)))))
    naive.num.reject <- naive.num.reject + list.resultJointTest$bool.rejectNaiveNull
    bonfe.num.reject <- bonfe.num.reject + list.resultJointTest$bool.rejectBonferroniNull
    joint.num.reject <- joint.num.reject + list.resultJointTest$bool.rejectJointNull

    if (bool.maxTest == FALSE) {
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
      if (!bool.mutePrint) {message(loop)}
    }
  }
  if (!bool.mutePrint) {
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
    if (bool.maxTest == FALSE) {joint.num.reject <- joint.num.reject[1,1]}
    print(joint.num.reject/ns)
  }

  summary <- list(vec.statMaxJoint=vec.statMaxJoint,chi.stat.vec=chi.stat.vec,effN.mean = mean(effN.z.vec),effNx.mean=mean(effN.x.vec),cor.z=mean(covZ1Z2.vec),naive=naive.num.reject/ns,bonfe=bonfe.num.reject/ns,joint=joint.num.reject/ns,maxCritical=vec.criticalValue)
  return(summary)
}


