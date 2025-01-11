#' return joint test result
#'
#' @param df_data Data.frame, data to evaluate.
#' @param int_dim_Z Integer, dimension of covariates.
#' @param int_J Integer, nearest neighbor for the variance estimation <= 3
#'  default is 3.
#' @param bool_max_test Boolean, use max test instead of L2 test,
#'  default is false.
#' @param bool_L2_std Boolean, for using standardized t stat for L2 test.
#' @param bool_joint Boolean, for returning joint test, instead of balance only.

return_result_joint <- function(df_data,
                               int_dim_Z,
                               int_J = 3,
                               bool_max_test = FALSE,
                               bool_L2_std = TRUE,
                               bool_joint = TRUE)
{
  int_num_simul_draw <- 3000

  # Zstat
  list_result_covariate <-
    compute_test_stat(df_data = df_data,
                    int_dim_Z = int_dim_Z,
                    bool_joint = bool_joint,
                    int_J = int_J,
                    bool_L2_std = bool_L2_std,
                    bool_max_test = bool_max_test)


  bool_reject_naive_null <-
    list_result_covariate$bool_reject_naive_null
  bool_reject_bonferroni_null <-
    list_result_covariate$bool_reject_bonferroni_null

  int_dim_total <- int_dim_Z + bool_joint
  # Xstat
  tryCatch({
    list_result_X <- rddensity::rddensity(X = df_data$vec_X)},
    error = function(e) {
      message(paste0("ERROR at rddensity"))
      message(e)
    }
  )

  # t stat
  real_tstat_X <- list_result_X$test$t_jk

  if (abs(real_tstat_X) > stats::qnorm(1- 0.025)) {
    bool_reject_naive_null <- TRUE
  }
  if (abs(real_tstat_X) > stats::qnorm(1 - 0.025/int_dim_total)) {
    bool_reject_bonferroni_null <- TRUE
  }

  if (bool_max_test==FALSE) {
    if (bool_L2_std) {
      if (bool_joint) {
        real_stat_L2_std_joint <-
          list_result_covariate$real_stat + real_tstat_X^2
      } else {
        real_stat_L2_std_joint <-
          list_result_covariate$real_stat
      }
      vec_random_draw_stat <- rep(0,int_num_simul_draw)
      mat_cor <- matrix(rep(0,int_dim_total^2),int_dim_total,int_dim_total)
      mat_cor[1:int_dim_Z,1:int_dim_Z] <- list_result_covariate$mat_cor
      if (bool_joint) {
        mat_cor[int_dim_total,int_dim_total] <- 1
      }
      mat_random_draw_stat <-
        return_random_stat(int_num_simul = int_num_simul_draw,
                           int_dim = int_dim_total,
                           mat_cor = mat_cor)
      real_pvalue <- 0
      for (l in seq(1:int_num_simul_draw)) {
        # Replicating the (signed) vector of statistics
        # by multiplying the corr matrix, then take square and max.
        vec_random_draw_stat[l] <- sum(mat_random_draw_stat[,l]^2)
        real_pvalue <- real_pvalue +
          (vec_random_draw_stat[l] >=
             real_stat_L2_std_joint)/int_num_simul_draw
      }
      bool_reject_joint_null <- (real_pvalue <= 0.05)
      real_critical_value_joint <- stats::quantile(vec_random_draw_stat,0.95)
      real_stat_joint <- real_stat_L2_std_joint

    } else {
      real_critical_value_joint <- stats::qchisq(.95, df=int_dim_total)
      if (bool_joint) {
        bool_reject_joint_null <-
          (list_result_covariate$real_stat + real_tstat_X^2
           >= real_critical_value_joint)
        real_stat_joint <- list_result_covariate$real_stat + real_tstat_X^2
        real_pvalue <- 1-stats::pchisq(real_stat_joint,df=int_dim_total)
      } else {
        bool_reject_joint_null <-
          (list_result_covariate$real_stat
           >= real_critical_value_joint)
        real_stat_joint <- list_result_covariate$real_stat
        real_pvalue <- 1-stats::pchisq(real_stat_joint,df=int_dim_total)
      }

    }
  } else {
    if (bool_joint) {
      real_stat_max_joint <-
        max(list_result_covariate$real_stat,real_tstat_X^2)
    } else {
      real_stat_max_joint <- max(list_result_covariate$real_stat)
    }
    vec_random_draw_stat <- rep(0,int_num_simul_draw)
    mat_cor <- matrix(rep(0,(int_dim_total)^2),int_dim_total,int_dim_total)
    mat_cor[1:int_dim_Z,1:int_dim_Z] <- list_result_covariate$mat_cor
    if (bool_joint) {
      mat_cor[int_dim_total,int_dim_total] <- 1
    }
    mat_random_draw_stat <-
      return_random_stat(int_num_simul = int_num_simul_draw,
                         int_dim = int_dim_total,
                         mat_cor = mat_cor)
    for (l in seq(1:int_num_simul_draw)) {
      # Replicating the (signed) vector of statistics
      # by multiplying the corr matrix, then take square and max.
      vec_random_draw_stat[l] <- max(mat_random_draw_stat[,l]^2)
    }
    real_stat_joint <- real_stat_max_joint
    real_pvalue <- mean((vec_random_draw_stat >= real_stat_max_joint))
    bool_reject_joint_null <- (real_pvalue <= 0.05)
    real_critical_value_joint <- stats::quantile(vec_random_draw_stat,0.95)
  }


  real_cov_Z1_Z2 <- NULL
  if (int_dim_Z > 1) {
    real_cov_Z1_Z2 <-
      stats::cor(df_data$vec_Z.1[df_data$vec_X <= stats::quantile(abs(df_data$vec_X),0.1)],
          df_data$vec_Z.2[df_data$vec_X <= stats::quantile(abs(df_data$vec_X),0.1)])}
  real_eff_N_x <- list_result_X$N$eff_left + list_result_X$N$eff_right

  list_result_joint_test <-
    list(real_stat = list_result_covariate$real_stat,
         real_tstat_X = real_tstat_X,
         real_mean_tstat_Z = list_result_covariate$real_mean_tstat_Z,
         real_median_tstat_Z = list_result_covariate$real_median_tstat_Z,
         real_max_abs_tstat_Z = list_result_covariate$real_max_abs_tstat_Z,
         real_effective_N_mean_Z = list_result_covariate$real_effective_N_mean_Z,
         real_eff_N_x = real_eff_N_x,
         real_cov_Z1_Z2 = real_cov_Z1_Z2,
         bool_reject_joint_null = bool_reject_joint_null,
         bool_reject_naive_null = bool_reject_naive_null,
         bool_reject_bonferroni_null = bool_reject_bonferroni_null,
         real_critical_value_joint = real_critical_value_joint,
         vec_tstat_Z_raw = list_result_covariate$vec_tstat_Z_raw,
         real_pvalue = real_pvalue,
         real_stat_joint = real_stat_joint,
         int_largest_stat = list_result_covariate$int_largest_stat,
         real_largest_stat = list_result_covariate$real_largest_stat)
  return(list_result_joint_test)
}


