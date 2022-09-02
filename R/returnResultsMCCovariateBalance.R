####
#' MC covarage / power check using simulated data
#'
#' Example code for simulation excercise
#'
#' @param jump z jump size
#' @param dim number of z variables
#' @param n sample size
#' @param cov correlation of Z variables
#' @param frac.jump fraction of z with jumps specified as jump. 1 if only the last one, 0 if none, 0 < frac < 1 if fraction of vector, default: null. 0
#' @param bool.mutePrint Boolean option for muting outputs
#'
#' @export
####

returnResultsMCCovariateBalance <- function(jump,dim,n,cov.z,frac.jump = 0,bool.mutePrint=FALSE)
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

  ns <- 300
  naive.num.reject <- 0
  bonfe.num.reject <- 0
  joint.num.reject <- 0
  chi.stat.vec <- rep(0,length(ns))
  effN.z.vec <- c(1:ns)
  covZ1Z2.vec <- c(1:ns)
  mean.jump.z.vec <- c(1:ns)
  median.jump.z.vec <- c(1:ns)
  max.jump.z.vec <- c(1:ns)

  options <- list(int.DGP = 2, int.n = n,real.pParam=0.8,real.jumpZ = jump,int.dimZ = dim,real.covZ = cov.z,real.fracJumpZ=frac.jump,fun.mu=function(x){return(x^2)})

  for (loop in 1:ns)
  {
    data <- RDTests::simulateDGPs(options)
    list.resultCovariateBalance <- returnResultsBalanceTest(df.data=data,int.dimZ=options$int.dimZ)
    covZ1Z2.vec[loop] <- list.resultCovariateBalance$real.covZ1Z2

    mean.jump.z.vec[loop] <- list.resultCovariateBalance$real.meanStatTZ
    median.jump.z.vec[loop] <- list.resultCovariateBalance$real.medianStatTZ
    max.jump.z.vec[loop] <- list.resultCovariateBalance$real.maxAbsStatTZ
    effN.z.vec[loop] <- list.resultCovariateBalance$real.effectiveNMeanZ
    chi.stat.vec[loop] <- list.resultCovariateBalance$real.statChiSq
    naive.num.reject <- naive.num.reject + list.resultCovariateBalance$bool.rejectNaiveNull
    bonfe.num.reject <- bonfe.num.reject + list.resultCovariateBalance$bool.rejectBonferroniNull
    joint.num.reject <- joint.num.reject + list.resultCovariateBalance$bool.rejectJointNull

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
    print("chi-sq rejection rate")
    print(joint.num.reject[1,1]/ns)

    message("average effective sample size")
    message(mean(effN.z.vec))

    message("average correlation of Z1 and Z2")
    message(mean(covZ1Z2.vec))

    message("mean z stat")
    message(mean(mean.jump.z.vec))

    message("median z stat")
    message(mean(median.jump.z.vec))

    message("max z stat")
    message(mean(max.jump.z.vec))

    message("naive rejection rate")
    message(naive.num.reject/ns)
    message("bonferroni rejection rate")
    message(bonfe.num.reject/ns)
    message("chi-sq rejection rate")
    message(joint.num.reject[1,1]/ns)
  }
  summary <- list(chi.stat.vec=chi.stat.vec,effN.mean = mean(effN.z.vec),cor.z=mean(covZ1Z2.vec),naive=naive.num.reject/ns,bonfe=bonfe.num.reject/ns,joint=joint.num.reject[1,1]/ns)
  return(summary)
}


