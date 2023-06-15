#' Simulating DGP
#'
#' Simulation code for DGP
#' @param list_option List, list of options that specify the DGP.
#'

simulate_DGP = function(list_option)
{
  # Reading options ----
  list_parameter <- load_option_simulation(list_option)

  # DGP3: this is the only dgp used
  if (list_parameter$int_DGP == 3)
  {
    # Generating the variables --
    ## X is a mixture of truncated normal
    # p*N(0,0.12;[0,1])+(1-p)*N(0,0.12;[-1,0))
    # (mocking the variance of 2*B(2,4)-1),
    # D = 1 then RHS, D = 0 then LHS, pParam = 1/2
    # then no jump. 1 > pParam > 1/2 means a jump.
    vec_D <- (stats::runif(list_parameter$int_n) <= list_parameter$real_p_param)
    vec_X <- vec_D *
      truncnorm::rtruncnorm(n = list_parameter$int_n,
                            a = 0,
                            b = 1,
                            mean = 0,
                            sd = 0.12) +
      (1 - vec_D) *
      truncnorm::rtruncnorm(n = list_parameter$int_n,
                            a = -1,
                            b = 1e-30,
                            mean = 0,
                            sd = 0.12)

    ## Z is a*(X >= 0) Sigma Ztil
    mat_Z_til <- matrix(
      stats::rnorm(list_parameter$int_n * list_parameter$int_dim_Z),
      list_parameter$int_n,
      list_parameter$int_dim_Z)
    mat_Z_Sigma <- diag(1,nrow = list_parameter$int_dim_Z,
                        ncol = list_parameter$int_dim_Z)
    for (i in 1:list_parameter$int_dim_Z)
    {
      for (j in 1:list_parameter$int_dim_Z)
      {
        if (i != j) {
          if (list_parameter$real_cov_Z >= 0) {
            mat_Z_Sigma[i,j] <- list_parameter$real_cov_Z
          } else {
            if (i == 1 | j == 1) {
              mat_Z_Sigma[i,j] <- list_parameter$real_cov_Z
            } else {
              mat_Z_Sigma[i,j] <- -list_parameter$real_cov_Z
            }
          }
        }
      }
    }
    bool_frac_eval_flag <- FALSE
    if (list_parameter$str_frac_jump_Z == "None") {
      vec_a <- rep(0,list_parameter$int_dim_Z)
      bool_frac_eval_flag <- TRUE
    }
    if (list_parameter$str_frac_jump_Z == "All") {
      vec_a <- rep(1,list_parameter$int_dim_Z) *
        list_parameter$real_a_param / list_parameter$int_dim_Z
      bool_frac_eval_flag <- TRUE
    }
    if (list_parameter$str_frac_jump_Z == "Half") {
      vec_a <- rep(0,list_parameter$int_dim_Z)
      vec_a[
        (round(0.5 * length(vec_a))+1):length(vec_a)
        ] <- list_parameter$real_a_param * (2 / list_parameter$int_dim_Z)
      bool_frac_eval_flag <- TRUE
    }
    if (list_parameter$str_frac_jump_Z == "Last") {
      vec_a <- rep(0,list_parameter$int_dim_Z)
      vec_a[length(vec_a)] <- list_parameter$real_a_param
      bool_frac_eval_flag <- TRUE
    }
    # for development debugging
    # if (bool_frac_eval_flag == FALSE) {
    #   stop("str_frac_jump_Z is taking an undefined value.")}

    mat_A <- t(matrix(rep(vec_a,list_parameter$int_n),
                      list_parameter$int_dim_Z,
                      list_parameter$int_n))
    mat_Z <- list_option$fun_mu(vec_X) +
      mat_A * vec_D + mat_Z_til %*% chol(mat_Z_Sigma)
    vec_M <- rep(0,list_parameter$int_n)
    vec_X_star <- rep(0,list_parameter$int_n)
  }

  # Returning a data.frame ---
  ## mat.Z is converted into a list of vectors, vec.Z.1, vec.Z.2,...
  if (list_parameter$int_dim_Z > 1) {
    df_DGP_simulated <- data.frame(vec_Z = mat_Z,
                                   vec_X = vec_X,
                                   vec_M = vec_M,
                                   vec_D = vec_D,
                                   vec_X_star = vec_X_star)
  } else {
    df_DGP_simulated <- data.frame(vec_Z.1 = mat_Z,
                                   vec_X = vec_X,
                                   vec_M = vec_M,
                                   vec_D = vec_D,
                                   vec_X_star = vec_X_star)
  }

  return(df_DGP_simulated)
}
