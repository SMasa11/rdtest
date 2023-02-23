####
#' Option loaders for simulations
#'
#' Loading options for simulations tuned up for specified DGPs
#'
#' @param list_option list
####

load_option_simulation = function(list_option)
{
  # options must be a list containing
  ##  n: simulation sample size
  ##  DGP: choice of DGP
  ##  x.jump: [used only for DGP1] parameter relating the size of density jump
  ##  z.jump: parameter relating the size of mean function jump
  ##  p: [used only for DGP2] fraction of M = 1
  ##  dim.z: number of Z variables

  ## Create a container list ----
  list_option_value <- list()
  list_option_value$int_n <- list_option$int_n
  list_option_value$int_DGP <- list_option$int_DGP

  ## DGP selector ----
  ### DGP 1 ----
  if (list_option$int_DGP == 1) {
    # parameter a determines x.jump solely
    if (list_option$real_jump_X > stats::dnorm(0)) {
      stop("error! jumpX in DGP1 cannot exceed the phi(0) value. \n")
    }
    if (list_option$real_jump_X <= 0){
      stop("error! jumpX cannot take exactly 0. \n")
    }
    list_option_value$real_a_param <-
      stats::qnorm(stats::pnorm(0) +
                     0.5 * (1 - list_option$real_jump_X/stats::dnorm(0)))
    # parameter b relates z.jump in complicated way, so set as is
    list_option_value$real_b_param <- list_option$real_jump_Z
    # degree of heteroskedasticity at the cut-off
    list_option_value$real_a2_param <- list_option$real_a2_param

    if (is.null(list_option$real_frac_jump_Z) == FALSE) {
      list_option_value$real_frac_jump_Z <- list_option$real_frac_jump_Z}
    if (is.null(list_option$int_dim_Z) == FALSE) {
      list_option_value$int_dim_Z <- list_option$int_dim_Z}
    if (is.null(list_option$real_cov_Z) == FALSE){
      list_option_value$real_cov_Z <- list_option$real_cov_Z}

    # DGP1 accepts a smooth function mu
    list_option_value$fun_mu <- list_option$fun_mu
  }
  ### DGP 2 ----
  if (list_option$int_DGP == 2){
    list_option_value$real_a_param <-
      list_option$real_jump_Z / list_option$real_p_param
    list_option_value$real_frac_jump_Z <- list_option$real_frac_jump_Z
    list_option_value$real_p_param <- list_option$real_p_param
    list_option_value$int_dim_Z <- list_option$int_dim_Z
    list_option_value$real_cov_Z <- list_option$real_cov_Z
  }

  ### DGP 3----
  if (list_option$int_DGP == 3) {
    list_option_value$real_a_param <- list_option$real_jump_Z
    list_option_value$str_frac_jump_Z <- list_option$str_frac_jump_Z
    list_option_value$real_p_param <- 0.5 + list_option$real_jump_X / 2
    list_option_value$int_dim_Z <- list_option$int_dim_Z
    list_option_value$real_cov_Z <- list_option$real_cov_Z
    list_option_value$fun_mu <- list_option$fun_mu
  }
  # returning the parameters -----
  return(list_option_value)
}
