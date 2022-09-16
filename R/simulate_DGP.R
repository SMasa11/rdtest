####
#' Simulating DGP
#'
#' Simulation code for DGP
#'
#' @param list_option list
#' @export
####

simulate_DGP = function(list_option)
{
  # Reading options ----
  list_parameter <- load_option_simulation(list_option)

  if ((list_parameter$int_DGP != 1) &
      (list_parameter$int_DGP != 2) &
      (list_parameter$int_DGP != 3))
    {stop("Currently DGP specifications are either 1, 2 or 3.")}

  # THESE DGPs are not used.
  # # DGP1
  # if (list_parameter$int_DGP == 1)
  # {
  #   # Generating the base varaibles ----
  #   vec_X_til <- stats::rnorm(list_parameter$int_n)
  #   vec_Z_til <- stats::rnorm(list_parameter$int_n,0,3)
  #
  #   # Z shock matrix with covariance
  #   mat_Z_til_shock <-
  #     matrix(stats::rnorm(list_parameter$int_n * list_parameter$int_dim_Z,0,1),
  #            list_parameter$int_n,
  #            list_parameter$int_dim_Z)
  #   mat_Z_Sigma <- diag(1,nrow = list_parameter$int_dim_Z,
  #                       ncol = list_parameter$int_dim_Z)
  #   for (i in 1:list_parameter$int_dim_Z)
  #   {
  #     for (j in 1:list_parameter$int_dim_Z)
  #     {
  #       if (i != j) {mat_Z_Sigma[i,j] <- list_parameter$real_cov_Z}
  #     }
  #   }
  #
  #   vec_X_star <- stats::rnorm(list_parameter$int_n)
  #   vec_eps <- stats::rnorm(list_parameter$int_n)
  #
  #
  #   # Indicator of manipulation
  #   vec_M <- 1*((vec_X_star <= 0) &
  #                 (list_parameter$real_a_param +
  #                    (list_parameter$real_b_param) * vec_Z_til + vec_eps <= 0))
  #   # If manipulated, then Xtil enters
  #   vec_X <- (1 - vec_M) * vec_X_star + vec_M * abs(vec_X_til)
  #
  #   # Set the vector of jump volumes in mean Z
  #   if (list_parameter$real_frac_jump_Z == 0) {
  #     vec_b <- rep(0,list_parameter$int_dim_Z)}
  #   # fracjump == 1: only the last value gets 1, and others 0
  #   if (list_parameter$real_frac_jump_Z == 1) {
  #     vec_b <- rep(0,list_parameter$int_dim_Z)
  #     vec_b[length(vec_b)] <- 1
  #   }
  #   # 0 < fracjump < 1: round(dim*frac)+1 to the end is 1
  #   ## 0.5, dim=3 -> 2+1:3 -> 3:3
  #   ## 0.5, dim=10 -> 5+1:10 -> 6:10
  #   ## 0.0001, dim < 5000 -> 0+1:dim -> full
  #   if ((list_parameter$real_frac_jump_Z > 0) &
  #       (list_parameter$real_frac_jump_Z < 1)) {
  #     vec_b <- rep(0,list_parameter$int_dim_Z)
  #     vec_b[
  #       (round(list_parameter$real_frac_jump_Z * length(vec_b))+1):length(vec_b)
  #       ] <- 1
  #   }
  #   # Reshape the vector b into a matrix
  #   mat_B <- t(matrix(rep(vec_b,list_parameter$int_n),
  #                     list_parameter$int_dim_Z,
  #                     list_parameter$int_n))
  #   # Construct the base vector of Z
  #   vec_Z_til <- mat_B * matrix(rep(vec_Z_til,list_parameter$int_dim_Z),
  #                               list_parameter$int_n,
  #                               list_parameter$int_dim_Z)
  #   # Adding the correlated shocks in Z variable
  #   mat_Z <-  list_option$fun_mu(2 * stats::qbeta(stats::pnorm(vec_X),2,2)) +
  #     (1 + list_parameter$real_a2_param*(vec_X >= 0) + abs(vec_X)) *
  #     (vec_Z_til  + mat_Z_til_shock %*% chol(mat_Z_Sigma))
  #   # D vector by the sharp RD definition
  #   vec_D <- (vec_X >= 0)
  # }
  #
  # # DGP2
  # if (list_parameter$int_DGP == 2)
  # {
  #   # Generating the variables ----
  #   mat_Z_til <- matrix(
  #     stats::rnorm(list_parameter$int_n * list_parameter$int_dim_Z),
  #     list_parameter$int_n,
  #     list_parameter$int_dim_Z)
  #   mat_Z_Sigma <- diag(1,nrow = list_parameter$int_dim_Z,
  #                       ncol = list_parameter$int_dim_Z)
  #   for (i in 1:list_parameter$int_dim_Z)
  #   {
  #     for (j in 1:list_parameter$int_dim_Z)
  #     {
  #       if (i != j) {mat_Z_Sigma[i,j] <- list_parameter$real_cov_Z}
  #     }
  #   }
  #
  #   vec_U <- stats::runif(list_parameter$int_n)
  #   vec_X_star <- stats::qnorm(vec_U)
  #   vec_X_til <- stats::qnorm(
  #     (vec_U + 0.5)*(vec_U < 0.5) + (vec_U - 0.5)*(vec_U >= 0.5))
  #
  #   vec_M <-
  #     (stats::runif(list_parameter$int_n) <= list_parameter$real_p_param)
  #   vec_X <- (1 - vec_M) * vec_X_star + vec_M * vec_X_til
  #   if (list_parameter$real_frac_jump_Z == 0) {
  #     vec_a <- rep(0,list_parameter$int_dim_Z)}
  #   if (list_parameter$real_frac_jump_Z == 1) {
  #     vec_a <- rep(0,list_parameter$int_dim_Z)
  #     vec_a[length(vec_a)] <- list_parameter$real_a_param
  #   }
  #   if ((list_parameter$real_frac_jump_Z > 0) &
  #       (list_parameter$real_frac_jump_Z < 1)) {
  #     vec_a <- rep(0,list_parameter$int_dim_Z)
  #     vec_a[
  #       (round(list_parameter$real_frac_jump_Z * length(vec_a))+1):length(vec_a)
  #       ] <- list_parameter$real_a_param
  #   }
  #
  #   #seq(0,parameter$a,length=(parameter$dim.z+1))[2:(parameter$dim.z+1)]
  #   mat_A <- t(matrix(rep(vec_a,list_parameter$int_n),
  #                     list_parameter$int_dim_Z,
  #                     list_parameter$int_n))
  #   mat_Z <- mat_A*
  #     matrix(rep((vec_U - 0.5),list_parameter$int_dim_Z),
  #            list_parameter$int_n,
  #            list_parameter$int_dim_Z) +
  #     mat_Z_til %*% chol(mat_Z_Sigma)
  #   vec_D <- (vec_X >= 0)
  # }

  # DGP3
  if (list_parameter$int_DGP == 3)
  {
    # Generating the variables --
    ## X is a mixture of truncated normal
    # p*N(0,0.12;[0,1])+(1-p)*N(0,0.12;[-1,0))
    # (mocking the variance of 2*B(2,4)-1),
    # D = 1 then RHS, D = 0 then LHS, pParam = 1/2
    # then no jump. 1 > pParam > 1/2 means a jump.
    vec_D <- (runif(list_parameter$int_n) <= list_parameter$real_p_param)
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
        if (i != j) {mat_Z_Sigma[i,j] <- list_parameter$real_cov_Z}
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
    if (bool_frac_eval_flag == FALSE) {
      stop("str_frac_jump_Z is taking an undefined value.")}

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
