####
#' Option loaders for simulations
#'
#' Loading options for simulations tuned up for specified DGPs
#'
#' @param list.options list
#' @export
####

loadOptionsSimulation = function(list.options)
{
  # options must be a list containing
  ##  n: simulation sample size
  ##  DGP: choice of DGP
  ##  x.jump: [used only for DGP1] parameter relating the size of density jump
  ##  z.jump: parameter relating the size of mean function jump
  ##  p: [used only for DGP2] fraction of M = 1
  ##  dim.z: number of Z variables

  ## Create a container list ----
  list.optionValues <- list()
  list.optionValues$int.n <- list.options$int.n
  list.optionValues$int.DGP <- list.options$int.DGP

  ## DGP selector ----
  ### DGP 1 ----
  if (list.options$int.DGP == 1) {
    # parameter a determines x.jump solely
    if (list.options$real.jumpX > dnorm(0)) {
      stop("error! jumpX in DGP1 cannot exceed the phi(0) value. \n")
    }
    if (list.options$real.jumpX <= 0){
      stop("error! jumpX cannot take exactly 0. \n")
    }
    list.optionValues$real.aParam <- qnorm(pnorm(0) + 0.5*(1 - list.options$real.jumpX/dnorm(0)))
    # parameter b relates z.jump in complicated way, so set as is
    list.optionValues$real.bParam <- list.options$real.jumpZ
    # degree of heteroskedasticity at the cut-off
    list.optionValues$real.a2Param <- list.options$real.a2Param

    if (is.null(list.options$real.fracJumpZ) == FALSE) {list.optionValues$real.fracJumpZ <- list.options$real.fracJumpZ}
    if (is.null(list.options$int.dimZ) == FALSE) {list.optionValues$int.dimZ <- list.options$int.dimZ}
    if (is.null(list.options$real.covZ) == FALSE){list.optionValues$real.covZ <- list.options$real.covZ}

    # DGP1 accepts a smooth function mu
    list.optionValues$fun.mu <- list.options$fun.mu
  }
  ### DGP 2 ----
  if (list.options$int.DGP == 2){
    list.optionValues$real.aParam <- list.options$real.jumpZ / list.options$real.pParam
    list.optionValues$real.fracJumpZ <- list.options$real.fracJumpZ
    list.optionValues$real.pParam <- list.options$real.pParam
    list.optionValues$int.dimZ <- list.options$int.dimZ
    list.optionValues$real.covZ <- list.options$real.covZ
  }

  ### DGP 3----
  if (list.options$int.DGP == 3) {
    list.optionValues$real.aParam <- list.options$real.jumpZ
    list.optionValues$str.fracJumpZ <- list.options$str.fracJumpZ
    list.optionValues$real.pParam <- 0.5 + list.options$real.jumpX/2
    list.optionValues$int.dimZ <- list.options$int.dimZ
    list.optionValues$real.covZ <- list.options$real.covZ
    list.optionValues$fun.mu <- list.options$fun.mu
  }
  # returning the parameters -----
  return(list.optionValues)
}
