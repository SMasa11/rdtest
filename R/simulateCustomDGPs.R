####
#' Simulating customized DGP using utilRDMC
#'
#'
#' @param int.typeDGP an integer representing the type of custom DGP
#' @param int.dim an integer of covariates dimensions
#' @param int.sampleSize an integer of sample size
#' @param real.TE real value of the discontinuity gap
#' @param real.corr real value of correlation coefficient across dep. vars
#' @export
####

simulateCustomDGPs <- function(int.typeDGP,int.dim,int.sampleSize,real.TE,real.corr)
{
  bool.flagKnown <- FALSE
  if (int.typeDGP == 1) {
    func.eval <- function(vec.X) {return(func.mu1(vec.X = vec.X,real.TE=real.TE))}
    df.DGP <- utilRDMC::generateMC(str.typeMu = "User",func.mu=func.eval,int.dim = int.dim,int.sampleSize = int.sampleSize,real.corr = real.corr)
    bool.flagKnown <- TRUE
  }
  if (int.typeDGP == 2) {
    func.eval <- function(vec.X) {return(func.mu2(vec.X = vec.X,real.TE=real.TE))}
    df.DGP <- utilRDMC::generateMC(str.typeMu = "User",func.mu=func.eval,int.dim = int.dim,int.sampleSize = int.sampleSize,real.corr = real.corr)
    bool.flagKnown <- TRUE
  }
  if (bool.flagKnown == FALSE) {stop(paste0(paste0("DGP type ",int.typeDGP)," is not specified."))}

  # Returning a data.frame ---
  ## mat.Z is converted into a list of vectors, vec.Z.1, vec.Z.2,...
  #df.DGPSimulated <- data.frame(vec.Z=mat.Z,vec.X=vec.X,vec.M=vec.M,vec.D=vec.D,vec.Xstar=vec.Xstar)
  #return(df.DGPSimulated)
  return (df.DGP)
}

func.mu1 <- function(vec.X,real.TE) {return((vec.X > 0)*(real.TE + vec.X) + (vec.X <= 0)*0) }

func.mu2 <- function(vec.X,real.TE) {return((vec.X > 0)*(real.TE + utilRDMC::specMuCCT(vec.X = vec.X)) + (vec.X <= 0)*0)}

