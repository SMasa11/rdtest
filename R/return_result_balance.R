#'
#' Return result of joint covariates balance test
#'
#' @param df_data A data.frame of data to evaluate.
#' @param int_dim_Z  An integer dimension of covariates.
#' @param int_J An integer of nearest neighbor for the variance estimation <= 3.
#' @param bool_max_test A boolean option to return max test rather than L2 test.
#' @param bool_max_test_inv_V A boolean option for using inv V for max test.
#' @param bool_L2_std A boolean option for using standardized t stat to conduct
#'   L2 test.
#'

return_result_balance <- function(df_data,
                                  int_dim_Z,
                                  int_J = 3,
                                  bool_max_test = FALSE,
                                  bool_max_test_inv_V = FALSE,
                                  bool_L2_std = FALSE)
{
  if (int_dim_Z < 2) {stop("Error: For the balance test with a single covariate, use rdrobust package directly.")}

  int_num_simulation_draw <- 3000

  real_cov_Z1_Z2 <- cor(df_data$vec_Z.1,df_data$vec_Z.2)

  if (bool_max_test == FALSE) {
    list_result_balance_part <- compute_stat_L2(df_data = df_data,
                                                int_dim_Z = int_dim_Z,
                                                bool_joint = FALSE,
                                                int_J = int_J)
  } else {
    list_result_balance_part <- compute_stat_max(df_data = df_data,
                                                  int_dim_Z = int_dim_Z,
                                                  bool_joint = FALSE,
                                                  int_J = int_J)
  }

  if (bool_max_test==FALSE) {
    real_stat_balance <- list_result_balance_part$real_stat_balance_L2
    real_critical_value_balance <- stats::qchisq(.95, df = int_dim_Z)
    bool_reject_balance <- (real_stat_balance >= real_critical_value_balance)
  } else {
    real_stat_balance <- list_result_balance_part$real_stat_balance_max
    vec_random_draw_stats <- rep(0,int_num_simulation_draw)
    # Generates the same random normal vector with the seed of 373568,
    # restoring the currently loading seed
    mat_random_Norm <- matrix(
        R.utils::withSeed(
          {rnorm((int_dim_Z+1)*int_num_simulation_draw,0,1)},seed=373568
          ),
        int_dim_Z,
        int_num_simulation_draw
      )
    mat_cor <- matrix(rep(0,(int_dim_Z)^2),(int_dim_Z),(int_dim_Z))
    mat_cor[1:int_dim_Z,1:int_dim_Z] <- list_result_balance_part$mat_cor
    U <- svd(mat_cor)$u
    V <- svd(mat_cor)$v
    D <- diag(sqrt(svd(mat_Cor)$d))
    mat_cor_sqrt <- U %*% D %*% t(V)
    for (l in seq(1:int_num_simulation_draw)) {
      # Replicating the (signed) vector of statistics by multiplying the corr matrix, then take square and max.
      vec_random_draw_stats[l] <- max((mat_cor_sqrt %*% mat_random_norm[,l])^2)
    }
    real_pvalue <- mean((vec_random_draw_stats >= real_stat_balance))
    bool_reject_balance <- (real_pvalue <= 0.05)
    real_critical_value_balance <- stats::quantile(vec_random_draw_stats,0.95)
  }


  list_result_balance <- list(bool_reject_balance = bool_reject_balance,
                              real_stat_balance = real_stat_balance,
                              real_mean_tstat_Z = list_result_balance_part$real_mean_tstat_Z,
                              real_median_tstat_Z = list_result_balance_part$real_median_tstat_Z,
                              real_max_abs_tstat_Z = list_result_balance_part$real_max_abs_tstat_Z,
                              real_cov_Z1_Z2 = real_cov_Z1_Z2,
                              real_mean_effective_N_Z = list_result_balance_part$real_mean_effective_N_Z,
                              real_critical_value_balance = real_critical_value_balance)
  return(list_result_balance)
}
