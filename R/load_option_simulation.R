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
  # deplicated
  list_option_value$int_DGP <- list_option$int_DGP

  list_option_value$real_a_param <- list_option$real_jump_Z
  list_option_value$str_frac_jump_Z <- list_option$str_frac_jump_Z
  list_option_value$real_p_param <- 0.5 + list_option$real_jump_X / 2
  list_option_value$int_dim_Z <- list_option$int_dim_Z
  list_option_value$real_cov_Z <- list_option$real_cov_Z
  list_option_value$fun_mu <- list_option$fun_mu

  # returning the parameters -----
  return(list_option_value)
}
