#' MC orthogonality check
#'
#' MC example of simulation study for orthogonality
#'
#' @param jump x jump size
#' @param a.2 heteroskedasticity added at cutoff
#' @param n sample size
#' @param bool.mutePrint boolean option to mute all the outputs
#'
#' @export

returnResultsMCOrthogonalStats <- function(jump,a.2,n,bool.mutePrint=FALSE)
{
  library("rdrobust")
  library("rddensity")
  library("ggplot2")
  set.seed(76276)

  ns <- 5000
  stat.Z.vec <- c(1:ns)
  stat.X.vec <- c(1:ns)
  effN.z <- c(1:ns)
  opth.z <- c(1:ns)
  effN.x <- c(1:ns)
  opth.x <- c(1:ns)
  var_l.z <-  c(1:ns)
  var_r.z <-  c(1:ns)

  # case1 DGP1 x.jump = 0.2, z.jump = 1
  # case2 DGP1 x.jump = 0.01, z.jump = 0
  # case3 DGP2 p = 0.001, z.jump = 0

  for (i in 1:ns) {

    # Draw simulation
    options <- list(int.DGP = 1, int.n = n,real.pParam=0.001,real.jumpX = jump,real.jumpZ = 0,real.fracJumpZ=0,int.dimZ = 1, real.a2Param = a.2,fun.mu=function(x){return(x^2)})
    data <- RDTests::simulateDGPs(options)

    Zresult <- rdrobust(y=data$vec.Z,x=data$vec.X,rho=1)
    effN.z[i] <- sum(Zresult$N_h)
    opth.z[i] <- Zresult$bws[1,1]
    # Zstat
    stat.Z <- Zresult$z[3]
    stat.Z.vec[i] <- stat.Z
    var_l.z[i] <- var(data$vec.Z[data$vec.X < 0])
    var_r.z[i] <- var(data$vec.Z[data$vec.X >= 0])

    Xresult <- rddensity(X=data$vec.X)
    effN.x[i] <- Xresult$N$eff_left + Xresult$N$eff_right
    # t stat
    stat.X <- Xresult$test$t_jk
    stat.X.vec[i] <- stat.X

    #if (i %% 250 == 0)
    #{
    #  message(i)
    #}
    if ((i == 100) | (i == 500) | (i == 1000) | (i == 2000) | (i == 5000))
    {
      #check <- data.frame(X=stat.X.vec[1:i],Z=stat.Z.vec[1:i])
      #print(summary(lm(X~Z,data=check)))
      if (!bool.mutePrint) {message(i)}
    }
  }

  check <- data.frame(X=stat.X.vec,Z=stat.Z.vec)
  if (!bool.mutePrint) {print(summary(lm(X~Z,data=check)))}

  plot1 <- ggplot(check,aes(x=X,y=Z))+
    stat_density2d(aes(fill=..level..), geom="polygon")
  #print(ggplot(check, aes(x = X)) + stat_density())
  #print(ggplot(check, aes(x = Z)) + stat_density())
  #print("Simulated histogram of joint test")
  #print(plot1)

  Xn <- rnorm(ns)
  Zn <- rnorm(ns)
  checkn <- data.frame(X = Xn, Z = Zn)

  #plot2 <- ggplot(checkn,aes(x=X,y=Z))+
  #  stat_density2d(aes(fill=..level..), geom="polygon")
  #print("Random sample of joint normal")
  #print(plot2)

  if (!bool.mutePrint) {
    print("average effective sample: x")
    print(mean(effN.x))
    print("average effective sample: z")
    print(mean(effN.z))

    print("mean variance of Z variable, X < c")
    print(mean(var_l.z))
    print("mean variance of Z variable, X >= c")
    print(mean(var_r.z))
  }

  list.testResult5 <- returnResultsChiIndependence(df.check=check,int.numPartition = 5)
  if (!bool.mutePrint) {print(list.testResult5)}

  list.testResult10 <- returnResultsChiIndependence(df.check=check,int.numPartition = 10)
  if (!bool.mutePrint) {print(list.testResult10)}

  list.testResult30 <- returnResultsChiIndependence(df.check=check,int.numPartition = 30)
  if (!bool.mutePrint) {print(list.testResult30)}

  list.resultMCOrthogonal <- list(plot1=plot1,list.testResult5=list.testResult5,list.testResult10=list.testResult10,list.testResult30=list.testResult30)
  return(list.resultMCOrthogonal)
}
