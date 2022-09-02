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
#' @param vec.critRange vector of critical values (in nominal size) to find the size adjustable critical value
#'
#' @export
####

returnResultsMCJointTestAdjustedSize <- function(int.ns = 300, x.jump,jump,a.2,dim,n,cov.z,frac.jump = "Half",bool.mutePrint=FALSE,int.J = 3,bool.maxTest=FALSE,bool.maxTestVinv=FALSE,bool.knownV=FALSE,fun.adjustBand=NULL,vec.critRange=vec.critRange)
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
  mat.num.reject <- rep(0,length(vec.critRange))

  #fun.mu <- function(x){return(x^2)}
  fun.mu <- function(x){return(0.48
                                 + (1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5)*(x < 0)
                                 + (0.84*x - 3.00*x^2 + 7.99*x^3 - 9.01*x^4 + 3.56*x^5)*(x >= 0))}

  options <- list(int.DGP = 3, int.n = n,real.jumpX = x.jump,real.jumpZ = jump, real.a2Param = a.2,int.dimZ = dim,real.covZ = cov.z,str.fracJumpZ=frac.jump,fun.mu=fun.mu)

  for (loop in 1:ns)
  {
    data <- RDTests::simulateDGPs(options)
    vec.rejectJointNull <- returnResultsJointTestAdjustedSize(df.data=data,int.dimZ=options$int.dimZ,int.J=int.J,bool.maxTest=bool.maxTest,bool.maxTestVinv=bool.maxTestVinv,bool.knownV=bool.knownV,real.covZ=cov.z,fun.adjustBand=fun.adjustBand,vec.critRange=vec.critRange)
    mat.num.reject <- mat.num.reject + vec.rejectJointNull
  }
  print("result")
  print(vec.critRange)
  print(mat.num.reject/ns)
  real.pval <- vec.critRange[which.min(abs(mat.num.reject/ns - 0.05))]
  print(real.pval)

  return(real.pval)
}


