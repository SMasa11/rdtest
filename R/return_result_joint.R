#'
#' return joint test result
#'
#' @inheritParams return_result_balance
#'
####

return_result_joint <- function(df_data,
                                int_dim_Z,
                                int_J = 3,
                                bool_max_test = FALSE,
                                bool_max_test_inv_V = FALSE,
                                bool_L2_std = FALSE)
{
  int_num_simul_draw <- 3000

  # Zstat
  list_result_balance <- return_result_balance(df_data = df_data,
                                               int_dim_Z = int_dim_Z,
                                               int_J = int_J,
                                               bool_max_test = bool_max_test,
                                               bool_max_test_inv_V = bool_max_test_inv_V,
                                               bool_L2_std = bool_L2_std)
  # Xstat
  list_result_X <- rddensity::rddensity(X=df_data$vec_X)
  # t stat
  real_tstat_X <- list_result_X$test$t_jk

  if (bool_max_test==FALSE) {
    real_stat_joint <- list_result_balance$real_stat_balance + real_tstat_X^2
    if (bool_L2_std) {
      vec_random_draw_stats <- rep(0,int_num_simulation_draw)
      # Generates the same random normal vector with the seed of 373568,
      # restoring the currently loading seed
      mat_random_Norm <- matrix(
        R.utils::withSeed(
          {rnorm((int_dim_Z+1)*int_num_simulation_draw,0,1)},seed=373568
        ),
        int_dim_Z+1,
        int_num_simulation_draw
      )
      mat_cor <- matrix(rep(0,(int_dim_Z+1)^2),(int_dim_Z+1),(int_dim_Z+1))
      mat_cor[1:int_dim_Z,1:int_dim_Z] <- list_result_balance_part$mat_cor
      mat_cor[(int.dimZ+1),(int.dimZ+1)] <- 1
      U <- svd(mat_cor)$u
      V <- svd(mat_cor)$v
      D <- diag(sqrt(svd(mat_Cor)$d))
      mat_cor_sqrt <- U %*% D %*% t(V)
      real_pvalue <- 0
      for (l in seq(1:int_num_simul_draw)) {
        # Replicating the (signed) vector of statistics
        # by multiplying the corr matrix, then take square and max.
        vec_random_draw_stat[l] <- sum((mat_cor_sqrt %*% mat_random_norm[,l])^2)
        real_pvalue <- real_pvalue + (vec_random_draw_stat[l] >= real_stat_joint)/int_num_simul_draw
      }
      bool_reject_joint <- (real_pvalue <= 0.05)
      real_critical_value_joint <- stats::quantile(vec_random_draw_stat,0.95)
    } else {
      real_critical_value_joint <- stats::qchisq(.95, df=(int.dimZ+1))
      bool_reject_joint <- (real_stat_joint >= real_critical_value_joint)
    }

  } else {
    if (bool.maxTestVinv == TRUE) {
      # standardized max stats
      real_stat_joint <- max(list_result_balance$real_stat_balance,real_tstat_X^2)
      real.criticalValueJoint <- stats::qchisq(0.95^(1/(int_dim_Z+1)),1)
      bool_reject_joint <- (real_stat_joint >= real.criticalValueJoint)
    } else {
      real_stat_joint <- max(list_result_balance$real_stat_balance,real_tstat_X^2)
      vec_random_draw_stats <- rep(0,int_num_simulation_draw)
      # Generates the same random normal vector with the seed of 373568,
      # restoring the currently loading seed
      mat_random_Norm <- matrix(
        R.utils::withSeed(
          {rnorm((int_dim_Z+1)*int_num_simulation_draw,0,1)},seed=373568
        ),
        int_dim_Z+1,
        int_num_simulation_draw
      )
      mat_cor <- matrix(rep(0,(int_dim_Z+1)^2),(int_dim_Z+1),(int_dim_Z+1))
      mat_cor[1:int_dim_Z,1:int_dim_Z] <- list_result_balance_part$mat_cor
      mat_cor[(int.dimZ+1),(int.dimZ+1)] <- 1
      U <- svd(mat_cor)$u
      V <- svd(mat_cor)$v
      D <- diag(sqrt(svd(mat_Cor)$d))
      mat_cor_sqrt <- U %*% D %*% t(V)
      for (l in seq(1:int_num_simul_draw)) {
        # Replicating the (signed) vector of statistics
        # by multiplying the corr matrix, then take square and max.
        vec_random_draw_stat[l] <- max((mat_cor_sqrt %*% mat_random_norm[,l])^2)
      }
      real_pvalue <- mean((vec_random_draw_stats >= real_stat_joint))
      bool_reject_joint <- (real_pvalue <= 0.05)
      real_critical_value_joint <- stats::quantile(vec_random_draw_stats,0.95)
    }
  }


  real_cov_Z1_Z2 <- NULL
  if (int_dim_Z > 1) {
      real_cov_Z1_Z2 <- stats::cor(
        df_data$vec_Z.1[df_data$vec_X <= stats::quantile(abs(df_data$vec_X),0.1)],
        df_data$vec_Z.2[df_data$vec_X <= stats::quantile(abs(df_data$vec_X),0.1)]
      )
    }
  real_effective_N_X <- list_result_X$N$eff_left + list_result_X$N$eff_right


  list_resultJointTest <- list(real_stat_jonit = real_stat_joint,
                               real_tstat_X = real_tstat_X,
                               real_mean_tstat_Z = list_result_balance$real_mean_tstat_Z,
                               real_median_tstat_Z = list_result_balance$real_median_tstat_Z,
                               real_maxAbsStatTZ = list_result_balance$real_max_abs_tstat_Z,
                               real_mean_effective_N_Z = list_result_balance$real_mean_effective_N_Z,
                               real_effective_N_X = real_effective_N_X,
                               real_cov_Z1_Z2 = real_cov_Z1_Z2,
                               bool_reject_joint = bool_reject_joint,
                               real_critical_joint = real_critical_joint)

  return(list.resultJointTest)
}


