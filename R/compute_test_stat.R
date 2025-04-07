#' Return statistics for covariates balance or placebo test.
#'
#' @inherit return_result_joint

compute_test_stat <- function(df_data,
                            int_dim_Z,
                            int_J = 3,
                            bool_max_test = FALSE,
                            bool_L2_std = FALSE,
                            bool_joint = TRUE,
                            bool_equivalence = FALSE,
                            real_epsilon = NULL)
{
  # booleans which return TRUE if tests are rejected by covariates only
  bool_reject_naive_null <- FALSE
  bool_reject_bonferroni_null <- FALSE

  # equivalence test is valid only for the max test form
  if (bool_max_test == FALSE) {
    bool_equivalence = FALSE
  }
  if (bool_equivalence & !is.null(real_epsilon)) {
    max_equivalence_statistic <- -Inf
    min_equivalence_statistic <- Inf
  } else {
    max_equivalence_statistic <- NA
    min_equivalence_statistic <- NA
  }
 
  real_effective_N_mean_Z <- 0

  vec_h_L <- rep(0,int_dim_Z)
  vec_h_R <- rep(0,int_dim_Z)
  vec_V_robust_L <- rep(0,int_dim_Z)
  vec_V_robust_R <- rep(0,int_dim_Z)

  mat_beta_L <- matrix(rep(0,int_dim_Z*3),int_dim_Z,3)
  mat_beta_R <- matrix(rep(0,int_dim_Z*3),int_dim_Z,3)

  vec_X_L <- df_data$vec_X[df_data$vec_X<0,drop=FALSE]
  vec_X_R <- df_data$vec_X[df_data$vec_X>=0,drop=FALSE]
  mat_dif_in_X_L <- matrix(rep(t(vec_X_L),length(vec_X_L)),
                           length(vec_X_L),
                           length(vec_X_L))
  mat_dif_in_X_L <- abs(mat_dif_in_X_L - t(mat_dif_in_X_L)) +
    diag(Inf,length(vec_X_L),length(vec_X_L))

  mat_index_close_L <- matrix(rep(0,length(vec_X_L)*int_J),
                              length(vec_X_L),int_J)
  for (j in 1:int_J) {
    mat_index_close_L[,j] <- max.col(-mat_dif_in_X_L)
    for (i in 1:length(vec_X_L)) {
      mat_dif_in_X_L[i,mat_index_close_L[i,j]] <- Inf
    }
  }

  mat_dif_in_X_R <- matrix(rep(t(vec_X_R),length(vec_X_R)),
                        length(vec_X_R),
                        length(vec_X_R))
  mat_dif_in_X_R <- abs(mat_dif_in_X_R - t(mat_dif_in_X_R)) +
    diag(Inf,length(vec_X_R),length(vec_X_R))

  mat_index_close_R <- matrix(rep(0,length(vec_X_R)*int_J),
                              length(vec_X_R),
                              int_J)
  for (j in 1:int_J) {
    mat_index_close_R[,j] <- max.col(-mat_dif_in_X_R)
    for (i in 1:length(vec_X_R)) {
      mat_dif_in_X_R[i,mat_index_close_R[i,j]] <- Inf
    }
  }

  # Arrays to store
  mat_X_h_L <- matrix(rep(0,length(vec_X_L)*int_dim_Z),
                      length(vec_X_L),
                      int_dim_Z)
  mat_X_h_R <- matrix(rep(0,length(vec_X_R)*int_dim_Z),
                      length(vec_X_R),
                      int_dim_Z)

  arr_Psi_L <- array(rep(0,3*length(vec_X_L)*int_dim_Z),
                    dim=c(length(vec_X_L),3,int_dim_Z))
  arr_Psi_R <- array(rep(0,3*length(vec_X_R)*int_dim_Z),
                    dim=c(length(vec_X_R),3,int_dim_Z))

  arr_G_inv_L <- array(rep(0,3^2*int_dim_Z),dim=c(3,3,int_dim_Z))
  arr_G_inv_R <- array(rep(0,3^2*int_dim_Z),dim=c(3,3,int_dim_Z))

  mat_residual_nn_L <- matrix(rep(0,length(vec_X_L)*int_dim_Z),
                            length(vec_X_L),
                            int_dim_Z)
  mat_residual_nn_R <- matrix(rep(0,length(vec_X_R)*int_dim_Z),
                            length(vec_X_R),
                            int_dim_Z)

  vec_tstat_Z <- rep(0,int_dim_Z)
  vec_se_Z <- rep(0,int_dim_Z)

  list_result_rdrobust_Z <- list()
  list_result_rdrobust_Z_bias <- list()
  for (d in 1:int_dim_Z)
  {
    # enforce the bandwidths are the same for both sides
    # Note: density bands choose different bands by default,
    # but the value seems to be sensitive to the option.
    # I'll keep them as their default.
    if (int_dim_Z == 1) {
      tryCatch({
        list_result_rdrobust_Z <- rdrobust::rdrobust(
          y=df_data$vec_Z.1,
          x = df_data$vec_X,
          rho=1,
          bwselect='mserd'
        )},
        error = function(e) {
          message("ERROR at rdrobust for vec_Z.1")
          message(e)
        }
      )
    } else {
      tryCatch({
        eval(parse(text = paste0(paste0(
          "list_result_rdrobust_Z <- rdrobust::rdrobust(y=df_data$vec_Z.",d),
          ",x = df_data$vec_X,rho=1,bwselect='mserd')"
        )))},
        error = function(e) {
          message(paste0("ERROR at rdrobust for vec_Z.",d))
          message(e)
        }
      )
    }

    # Zstat
    real_tstat_Z <- list_result_rdrobust_Z$z[3]
    vec_tstat_Z[d] <- real_tstat_Z
    vec_se_Z[d] <- list_result_rdrobust_Z$se[3]
    vec_h_L[d] <- list_result_rdrobust_Z$bws[1,1]
    vec_h_R[d] <- list_result_rdrobust_Z$bws[1,2]
    real_effective_N_mean_Z <- real_effective_N_mean_Z +
      sum(list_result_rdrobust_Z$N_h)/int_dim_Z

    tryCatch({
      eval(parse(text = paste0(paste0(
        "list_result_rdrobust_Z_bias <- rdrobust::rdrobust(y=df_data$vec_Z.",d),
        ",x=df_data$vec_X,p=2,rho=1,h=c(vec_h_L[d],vec_h_R[d]))"
        )))},
      error = function(e) {
        message(paste0("ERROR at calculating bias rdrobust for vec_Z.",d))
        message(e)
      }
    )

    vec_V_robust_L[d] <- list_result_rdrobust_Z$V_rb_l[1,1]
    vec_V_robust_R[d] <- list_result_rdrobust_Z$V_rb_r[1,1]
    if (abs(real_tstat_Z) > stats::qnorm(1- 0.025)) {
      bool_reject_naive_null <- TRUE
    }
    if (abs(real_tstat_Z) > stats::qnorm(1 - 0.025/(int_dim_Z+bool_joint))){
      bool_reject_bonferroni_null <- TRUE
    }

    # For each d, compute X(h_d)
    vec_X_h_L <- df_data$vec_X[df_data$vec_X< 0,drop=FALSE]/vec_h_L[d]
    vec_X_h_R <- df_data$vec_X[df_data$vec_X>=0,drop=FALSE]/vec_h_R[d]

    mat_X_h_L_in_L <- cbind(rep(1,length(vec_X_h_L)),vec_X_h_L)
    mat_X_h_L_in_R <- cbind(rep(1,length(vec_X_h_R)),vec_X_h_R)

    mat_X_h_sq_L <- cbind(rep(1,length(vec_X_h_L)),vec_X_h_L,vec_X_h_L^2)
    mat_X_h_sq_R <- cbind(rep(1,length(vec_X_h_R)),vec_X_h_R,vec_X_h_R^2)

    mat_K_h_L <- diag(((1-abs(vec_X_h_L))*(abs(vec_X_h_L)<=1))/vec_h_L[d])
    mat_K_h_R <- diag(((1-abs(vec_X_h_R))*(abs(vec_X_h_R)<=1))/vec_h_R[d])

    mat_G_inv_L_in_L <-
      chol2inv(chol(crossprod(sqrt(mat_K_h_L) %*% mat_X_h_L_in_L))) *
      length(vec_X_h_L)
    mat_G_inv_L_in_R <-
      chol2inv(chol(crossprod(sqrt(mat_K_h_R) %*% mat_X_h_L_in_R))) *
      length(vec_X_h_R)

    mat_G_inv_sq_L <-
      chol2inv(chol(crossprod(sqrt(mat_K_h_L) %*% mat_X_h_sq_L))) *
      length(vec_X_h_L)
    mat_G_inv_sq_R <-
      chol2inv(chol(crossprod(sqrt(mat_K_h_R) %*% mat_X_h_sq_R))) *
      length(vec_X_h_R)

    arr_G_inv_L[,,d] <- mat_G_inv_sq_L
    arr_G_inv_R[,,d] <- mat_G_inv_sq_R

    vec_Z_L <- eval(parse(text = paste0(
      paste0("df_data$vec_Z.",d),"[df_data$vec_X < 0]")))
    vec_Z_R <- eval(parse(text = paste0(
      paste0("df_data$vec_Z.",d),"[df_data$vec_X >= 0]")))

    mat_beta_L[d,] <- list_result_rdrobust_Z_bias$beta_Y_p_l
    mat_beta_R[d,] <- list_result_rdrobust_Z_bias$beta_Y_p_r

    # This residual should be the raw conditional mean projection errors of Z,
    # not normalized by bandwidths
    vec_residual_L <- vec_Z_L -
      (mat_beta_L[d,1] + mat_beta_L[d,2]*vec_X_L + mat_beta_L[d,3]*vec_X_L^2)
    vec_residual_R <- vec_Z_R -
      (mat_beta_R[d,1] + mat_beta_R[d,2]*vec_X_R + mat_beta_R[d,3]*vec_X_R^2)

    # Initialize the containers of near-neighbor estimates of residuals
    vec_residual_nn_L <- vec_residual_L
    vec_residual_nn_R <- vec_residual_R
    for (i in 1:length(vec_X_L))
    {
      for (j in 1:int_J)
      {
        vec_residual_nn_L[i] <-
          vec_residual_nn_L[i] - vec_residual_L[mat_index_close_L[i,j]]/int_J
      }
    }
    for (i in 1:length(vec_X_R))
    {
      for (j in 1:int_J)
      {
        vec_residual_nn_R[i] <-
          vec_residual_nn_R[i] - vec_residual_R[mat_index_close_R[i,j]]/int_J
      }
    }
    arr_Psi_L[,,d] <- diag(vec_residual_nn_L) %*% (mat_K_h_L %*% mat_X_h_sq_L)
    arr_Psi_R[,,d] <- diag(vec_residual_nn_R) %*% (mat_K_h_R %*% mat_X_h_sq_R)

  }
  vec_tstat_Z_raw <- vec_tstat_Z
  real_mean_tstat_Z <- mean(vec_tstat_Z)
  real_median_tstat_Z <- stats::median(vec_tstat_Z)
  real_max_abs_tstat_Z <- max(abs(vec_tstat_Z))

  mat_cov_L <- matrix(rep(0,int_dim_Z^2),int_dim_Z,int_dim_Z)
  mat_cov_R <- matrix(rep(0,int_dim_Z^2),int_dim_Z,int_dim_Z)
  for (d1 in 1:int_dim_Z)
  {
    for (d2 in 1:int_dim_Z)
    {
      if (d2 >= d1)
      {
        mat_Psi_L <-
          (int_J/(int_J+1)) *
          t(arr_Psi_L[,,d1]) %*% arr_Psi_L[,,d2] / length(vec_X_L)
        mat_Psi_R <-
          (int_J/(int_J+1)) *
          t(arr_Psi_R[,,d1]) %*% arr_Psi_R[,,d2] / length(vec_X_R)

        mat_Var_L <-
          arr_G_inv_L[,,d1] %*% mat_Psi_L %*% arr_G_inv_L[,,d2]/length(vec_X_L)
        mat_Var_R <-
          arr_G_inv_R[,,d1] %*% mat_Psi_R %*% arr_G_inv_R[,,d2]/length(vec_X_R)

        # Multiply the hjk = (hj * hk)^(1/2)
        mat_cov_L[d1,d2] <- mat_Var_L[1,1] * sqrt(vec_h_L[d1] * vec_h_L[d2])
        mat_cov_R[d1,d2] <- mat_Var_R[1,1] * sqrt(vec_h_R[d1] * vec_h_R[d2])
      } else {
        mat_cov_L[d1,d2] <- mat_cov_L[d2,d1]
        mat_cov_R[d1,d2] <- mat_cov_R[d2,d1]
      }
    }
  }
  mat_cov <- mat_cov_L + mat_cov_R
  mat_cor <- stats::cov2cor(mat_cov)
  # Z.stat is standardized vector: multiply the se back
  # also, multiply hk^(1/2): we take vec_hL here,
  # but not that bandwidths are chosen to be the same for both sides.

  if (bool_max_test) {
    for (l in 1:int_dim_Z) {
      vec_tstat_Z[l] <-
        vec_tstat_Z[l] * vec_se_Z[l] * (vec_h_L[l]^(1/2))/sqrt(mat_cov[l,l])
      if (bool_equivalence & !is.null(real_epsilon)) {
        real_tau_minus <- 
          (vec_tstat_Z[l]*vec_se_Z[l] - real_epsilon) * (vec_h_L[l]^(1/2))/sqrt(mat_cov[l,l])
        real_tau_plus <-
          (vec_tstat_Z[l]*vec_se_Z[l] + real_epsilon) * (vec_h_L[l]^(1/2))/sqrt(mat_cov[l,l])
        max_equivalence_statistic <-
          max(max_equivalence_statistic,real_tau_minus)
        min_equivalence_statistic <-
          min(min_equivalence_statistic,real_tau_plus)
      }
    }
    real_stat <- max(vec_tstat_Z^2)
  } else {
    ## Why are we multiplying the bandwidth here?
    # because we are allowing for different bands by l,
    # but doing so through Var matrix adjustment.
    for (l in 1:int_dim_Z) {
      vec_tstat_Z[l] <- vec_tstat_Z[l] * vec_se_Z[l] * (vec_h_L[l]^(1/2))
    }

    if (bool_L2_std == TRUE) {
      # standardize the statistics by the SE of the constructed Cov matrix
      vec_tstat_Z_std <- vec_tstat_Z
      for (l in 1:int_dim_Z) {
        vec_tstat_Z_std[l] <- vec_tstat_Z_std[l]/sqrt(mat_cov[l,l])
      }
      # this does not follow chisq, hence we need to simulate the critical value.
      real_stat <- t(vec_tstat_Z_std) %*% vec_tstat_Z_std
    } else {
      mat_cor <- NA
      real_stat <-
        t(vec_tstat_Z) %*% chol2inv(chol(mat_cov)) %*% vec_tstat_Z
    }
  }

  list_result_covariate_part <- list(
    real_stat = real_stat,
    mat_cor = mat_cor,
    bool_L2_std = bool_L2_std,
    bool_reject_bonferroni_null = bool_reject_bonferroni_null,
    bool_reject_naive_null = bool_reject_naive_null,
    real_mean_tstat_Z = real_mean_tstat_Z,
    real_median_tstat_Z = real_median_tstat_Z,
    real_max_abs_tstat_Z = real_max_abs_tstat_Z,
    real_effective_N_mean_Z=real_effective_N_mean_Z,
    vec_tstat_Z_raw = vec_tstat_Z_raw,
    max_equivalence_statistic = max_equivalence_statistic,
    min_equivalence_statistic = min_equivalence_statistic)
  return(list_result_covariate_part)
}
