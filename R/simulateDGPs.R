####
#' Simulating DGP
#'
#' Simulation code for DGP type 1
#'
#' @param list.options list
#' @export
####

simulateDGPs = function(list.options)
{
  # Reading options ----
  list.parameter <- loadOptionsSimulation(list.options)

  if ((list.parameter$int.DGP != 1) & (list.parameter$int.DGP != 2) & (list.parameter$int.DGP != 3)) {stop("Currently DGP specifications are either 1, 2 or 3.")}

  # DGP1
  if (list.parameter$int.DGP == 1)
  {
    # Generating the base varaibles ----
    vec.Xtil <- rnorm(list.parameter$int.n)
    vec.Ztil <- rnorm(list.parameter$int.n,0,3)

    # Z shock matrix with covariance
    mat.ZtilShock <- matrix(rnorm(list.parameter$int.n*list.parameter$int.dimZ,0,1),list.parameter$int.n,list.parameter$int.dimZ)
    mat.ZSigma <- diag(1,nrow=list.parameter$int.dimZ,ncol=list.parameter$int.dimZ)
    for (i in 1:list.parameter$int.dimZ)
    {
      for (j in 1:list.parameter$int.dimZ)
      {
        if (i != j) {mat.ZSigma[i,j] <- list.parameter$real.covZ}
      }
    }

    vec.Xstar <- rnorm(list.parameter$int.n)
    vec.eps <- rnorm(list.parameter$int.n)


    # Indicator of manipulation
    vec.M <- 1*((vec.Xstar <= 0) & (list.parameter$real.aParam + (list.parameter$real.bParam)*vec.Ztil + vec.eps <= 0))
    # If manipulated, then Xtil enters
    vec.X <- (1 - vec.M)*vec.Xstar + vec.M*abs(vec.Xtil)

    # Set the vector of jump volumes in mean Z
    if (list.parameter$real.fracJumpZ == 0) {vec.b <- rep(0,list.parameter$int.dimZ)}
    # fracjump == 1: only the last value gets 1, and others 0
    if (list.parameter$real.fracJumpZ == 1) {
      vec.b <- rep(0,list.parameter$int.dimZ)
      vec.b[length(vec.b)] <- 1
    }
    # 0 < fracjump < 1: round(dim*frac)+1 to the end is 1
    ## 0.5, dim=3 -> 2+1:3 -> 3:3
    ## 0.5, dim=10 -> 5+1:10 -> 6:10
    ## 0.0001, dim < 5000 -> 0+1:dim -> full
    if ((list.parameter$real.fracJumpZ > 0) & (list.parameter$real.fracJumpZ < 1)) {
      vec.b <- rep(0,list.parameter$int.dimZ)
      vec.b[(round(list.parameter$real.fracJumpZ*length(vec.b))+1):length(vec.b)] <- 1
    }
    # Reshape the vector b into a matrix
    mat.B <- t(matrix(rep(vec.b,list.parameter$int.n),list.parameter$int.dimZ,list.parameter$int.n))
    # Construct the base vector of Z
    vec.Ztil <- mat.B*matrix(rep(vec.Ztil,list.parameter$int.dimZ),list.parameter$int.n,list.parameter$int.dimZ)
    # Adding the correlated shocks in Z variable
    #mat.ratioMuX <- t(matrix(rep(rnorm(list.parameter$int.dimZ,1,5),list.parameter$int.n),list.parameter$int.dimZ,list.parameter$int.n))
    mat.Z <-  list.options$fun.mu(2*qbeta(pnorm(vec.X),2,2)) + (1 + list.parameter$real.a2Param*(vec.X >= 0) + abs(vec.X)) * (vec.Ztil  + mat.ZtilShock %*% chol(mat.ZSigma))
    # D vector by the sharp RD definition
    vec.D <- (vec.X >= 0)
  }

  # DGP2
  if (list.parameter$int.DGP == 2)
  {
    # Generating the variables ----
    mat.Ztil <- matrix(rnorm(list.parameter$int.n*list.parameter$int.dimZ),list.parameter$int.n,list.parameter$int.dimZ)
    mat.ZSigma <- diag(1,nrow=list.parameter$int.dimZ,ncol=list.parameter$int.dimZ)
    for (i in 1:list.parameter$int.dimZ)
    {
      for (j in 1:list.parameter$int.dimZ)
      {
        if (i != j) {mat.ZSigma[i,j] <- list.parameter$real.covZ}
      }
    }

    vec.U <- runif(list.parameter$int.n)
    vec.Xstar <- qnorm(vec.U)
    vec.Xtil <- qnorm((vec.U + 0.5)*(vec.U < 0.5) + (vec.U - 0.5)*(vec.U >= 0.5))

    vec.M <- (runif(list.parameter$int.n) <= list.parameter$real.pParam)
    vec.X <- (1 - vec.M)*vec.Xstar + vec.M*vec.Xtil
    if (list.parameter$real.fracJumpZ == 0) {vec.a <- rep(0,list.parameter$int.dimZ)}
    if (list.parameter$real.fracJumpZ == 1) {
      vec.a <- rep(0,list.parameter$int.dimZ)
      vec.a[length(vec.a)] <- list.parameter$real.aParam
    }
    if ((list.parameter$real.fracJumpZ > 0) & (list.parameter$real.fracJumpZ < 1)) {
      vec.a <- rep(0,list.parameter$int.dimZ)
      vec.a[(round(list.parameter$real.fracJumpZ*length(vec.a))+1):length(vec.a)] <- list.parameter$real.aParam
    }

    #seq(0,parameter$a,length=(parameter$dim.z+1))[2:(parameter$dim.z+1)]
    mat.A <- t(matrix(rep(vec.a,list.parameter$int.n),list.parameter$int.dimZ,list.parameter$int.n))
    mat.Z <- mat.A*matrix(rep((vec.U - 0.5),list.parameter$int.dimZ),list.parameter$int.n,list.parameter$int.dimZ) + mat.Ztil %*% chol(mat.ZSigma)
    vec.D <- (vec.X >= 0)

  }


  # DGP3
  if (list.parameter$int.DGP == 3)
  {
    # Generating the variables --
    ## X is a mixture of truncated normal p*N(0,0.12;[0,1])+(1-p)*N(0,0.12;[-1,0)) (mocking the variance of 2*B(2,4)-1),
    # D = 1 then RHS, D = 0 then LHS, pParam = 1/2 then no jump. 1 > pParam > 1/2 means a jump.
    vec.D <- (runif(list.parameter$int.n) <= list.parameter$real.pParam)
    vec.X <- vec.D*truncnorm::rtruncnorm(n=list.parameter$int.n,a=0,b=1,mean=0,sd=0.12) + (1 - vec.D)*truncnorm::rtruncnorm(n=list.parameter$int.n,a=-1,b=1e-30,mean=0,sd=0.12)

    ## Z is a*(X >= 0) Sigma Ztil
    mat.Ztil <- matrix(rnorm(list.parameter$int.n*list.parameter$int.dimZ),list.parameter$int.n,list.parameter$int.dimZ)
    mat.ZSigma <- diag(1,nrow=list.parameter$int.dimZ,ncol=list.parameter$int.dimZ)
    for (i in 1:list.parameter$int.dimZ)
    {
      for (j in 1:list.parameter$int.dimZ)
      {
        if (i != j) {mat.ZSigma[i,j] <- list.parameter$real.covZ}
      }
    }
    bool.fracEvalFlag <- FALSE
    if (list.parameter$str.fracJumpZ == "None") {
      vec.a <- rep(0,list.parameter$int.dimZ)
      bool.fracEvalFlag <- TRUE
    }
    if (list.parameter$str.fracJumpZ == "All") {
      vec.a <- rep(1,list.parameter$int.dimZ)*list.parameter$real.aParam/list.parameter$int.dimZ
      bool.fracEvalFlag <- TRUE
    }
    if (list.parameter$str.fracJumpZ == "Half") {
      vec.a <- rep(0,list.parameter$int.dimZ)
      vec.a[(round(0.5*length(vec.a))+1):length(vec.a)] <- list.parameter$real.aParam*(2/list.parameter$int.dimZ)
      bool.fracEvalFlag <- TRUE
    }
    if (list.parameter$str.fracJumpZ == "Last") {
      vec.a <- rep(0,list.parameter$int.dimZ)
      vec.a[length(vec.a)] <- list.parameter$real.aParam
      bool.fracEvalFlag <- TRUE
    }
    if (bool.fracEvalFlag == FALSE) {stop("str.fracJumpZ is taking an undefined value.")}

    mat.A <- t(matrix(rep(vec.a,list.parameter$int.n),list.parameter$int.dimZ,list.parameter$int.n))
    mat.Z <- list.options$fun.mu(vec.X) + mat.A*vec.D + mat.Ztil %*% chol(mat.ZSigma)
    vec.M <- rep(0,list.parameter$int.n)
    vec.Xstar <- rep(0,list.parameter$int.n)
  }

  # Returning a data.frame ---
  ## mat.Z is converted into a list of vectors, vec.Z.1, vec.Z.2,...
  if (list.parameter$int.dimZ > 1) {
    df.DGPSimulated <- data.frame(vec.Z=mat.Z,vec.X=vec.X,vec.M=vec.M,vec.D=vec.D,vec.Xstar=vec.Xstar)
  } else {
    df.DGPSimulated <- data.frame(vec.Z.1=mat.Z,vec.X=vec.X,vec.M=vec.M,vec.D=vec.D,vec.Xstar=vec.Xstar)
  }

  return(df.DGPSimulated)
}
