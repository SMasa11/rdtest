#' MC joint test of coverage / power check using simulated data
#'
#' Example code for simulation excercise
#'
#' @param int_ns Interger, number of simulations to run
#'  default is 300.
#' @param jump Real, z jump size.
#' @param x_jump Real, x jump size.
#' @param a_2 Real, gap in heteroskedasticity.
#' @param dim Integer, number of z variables.
#' @param n Integer, sample size.
#' @param cov_z Real, correlation of Z variables.
#' @param frac_jump String, either "None", "Half", "All", "Last".
#' @param bool_mutePrint Boolean, option for muting outputs.
#'  default is FALSE
#' @param int_J Integer, nearest neighbor for the variance estimation <= 3
#'  default is 3.
#' @param bool_max_test Boolean, use max test instead of L2 test,
#'  default is false.
#' @param bool_L2_std Boolean, for using standardized t stat for L2 test.
#' @param bool_joint Boolean, for running joint test with density or not.
#' @param bool_skip_rwolf Boolean, for skipping rwolf2 call,
#'  which is time consuming and requires Stata to run. Default is FALSE,
#'  but effective only when bool_joint is FALSE as well.
#'
#' @export

return_result_MC_joint <- function(int_ns = 300,
                                     x_jump,
                                     jump,
                                     a_2,
                                     dim,
                                     n,
                                     cov_z,
                                     frac_jump = "Half",
                                     bool_mutePrint = FALSE,
                                     int_J = 3,
                                     bool_max_test = FALSE,
                                     bool_L2_std = FALSE,
                                     bool_joint = TRUE,
                                     bool_skip_rwolf = FALSE)
{
  set.seed(52622)

  # if (jump != 0)
  # {
  #   if (frac_jump == 0)
  #   {
  #     warning(
  #       "WARNING:
  #       JUMP IS SPECIFIED BUT NOT APPLIED AS frac.jump IS NOT SPECIFIED.")
  #   }
  # }

  ns <- int_ns
  if (bool_joint == FALSE & bool_skip_rwolf == FALSE) {
    rwolf.num.reject <- 0
  } else {
    rwolf.num.reject <- NA
  }
  naive.num.reject <- 0
  bonfe.num.reject <- 0
  joint.num.reject <- 0
  chi.stat.vec <- rep(0,length(ns))
  vec.statMaxJoint <- rep(0,length(ns))
  effN.z.vec <- c(1:ns)
  effN.x.vec <- c(1:ns)
  covZ1Z2.vec <- c(1:ns)
  mean.jump.z.vec <- c(1:ns)
  median.jump.z.vec <- c(1:ns)
  max.jump.z.vec <- c(1:ns)
  vec.criticalValue <- c(1:ns)

  #fun.mu <- function(x){return(x^2)}
  fun_mu <- function(x){
    return(
      0.48
      + (1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5)*(x < 0)
      + (0.84*x - 3.00*x^2 + 7.99*x^3 - 9.01*x^4 + 3.56*x^5)*(x >= 0))}

  option <- list(int_DGP = 3,
                  int_n = n,
                  real_jump_X = x_jump,
                  real_jump_Z = jump,
                  real_a_2_param = a_2,
                  int_dim_Z = dim,
                  real_cov_Z = cov_z,
                  str_frac_jump_Z = frac_jump,
                  fun_mu = fun_mu)

  for (loop in 1:ns)
  {
    if (bool_mutePrint == FALSE) {
      message("loop: ", loop)
    }
    data <- simulate_DGP(option)

    if (option$int_dim_Z == 1) {
      list_result_joint <- rdtest(Z = data.frame(vec_Z.1=data$vec_Z.1),
                                  vec_X = data$vec_X,
                                  bool_joint = bool_joint,
                                  int_J = int_J,
                                  real_cutoff = 0,
                                  bool_max_test = bool_max_test,
                                  bool_L2_std = bool_L2_std)
    } else {
      list_result_joint <- rdtest(Z = data[,1:option$int_dim_Z],
                                  vec_X = data[,option$int_dim_Z+1],
                                  bool_joint = bool_joint,
                                  int_J = int_J,
                                  real_cutoff = 0,
                                  bool_max_test = bool_max_test,
                                  bool_L2_std = bool_L2_std)
    }

    if (bool_joint == FALSE) {
      if (bool_skip_rwolf == FALSE) {
        # run rwolf2 through stata to get the data.frame of p-values
        if (option$int_dim_Z == 3) {
          df_pvalues <-
            RStata::stata(
              "caller_rwolf2_3.do",
              data.in=data.frame(
                z1 = data[,1],
                z2 = data[,2],
                z3 = data[,3],
                x = data[,4]
              ),
              data.out=TRUE
            )
            min_pval <- min(
              c(df_pvalues$rw_pval_1,
                df_pvalues$rw_pval_2,
                df_pvalues$rw_pval_3)
            )
        } else {
          if (option$int_dim_Z == 5) {
            df_pvalues <-
              RStata::stata(
                "caller_rwolf2_5.do",
                data.in=data.frame(
                  z1 = data[,1],
                  z2 = data[,2],
                  z3 = data[,3],
                  z4 = data[,4],
                  z5 = data[,5],
                  x = data[,6]
                ),
                data.out=TRUE
              )
            min_pval <- min(
              c(df_pvalues$rw_pval_1,
                df_pvalues$rw_pval_2,
                df_pvalues$rw_pval_3,
                df_pvalues$rw_pval_4,
                df_pvalues$rw_pval_5)
            )
          } else {
            stop("Only dim = 3 or 5 is supported for the current set-up using rwolf2")
          }
        }
        print(min_pval)
        print((min_pval < 0.05))
        rwolf.num.reject <-
          rwolf.num.reject + (min_pval < 0.05)
      } else {
        rwolf.num.reject <- NA
      }
    }

    naive.num.reject <-
      naive.num.reject + list_result_joint$bool_reject_naive_null
    bonfe.num.reject <-
      bonfe.num.reject + list_result_joint$bool_reject_bonferroni_null
    joint.num.reject <-
      joint.num.reject + list_result_joint$bool_reject_joint_null

    if (bool_mutePrint == FALSE) {
      message("naive: ", naive.num.reject/loop)
      message("bonfe: ", bonfe.num.reject/loop)
      message("joint: ", joint.num.reject/loop)
      message("rwolf: ", rwolf.num.reject/loop)
    }

    if (bool_max_test == FALSE) {
      chi.stat.vec[loop] <- list_result_joint$real_stat
    } else {
      vec.statMaxJoint[loop]  <- list_result_joint$real_stat
    }
    mean.jump.z.vec[loop] <- list_result_joint$real_mean_tstat_Z
    median.jump.z.vec[loop] <- list_result_joint$real_median_tstat_Z
    max.jump.z.vec[loop] <- list_result_joint$real_max_abs_tstat_Z
    effN.z.vec[loop] <- list_result_joint$real_effective_N_mean_Z
    effN.x.vec[loop] <- list_result_joint$real_eff_N_x
    if (dim > 1) {
      covZ1Z2.vec[loop] <- list_result_joint$real_cov_Z1_Z2
    }

    vec.criticalValue[loop] <- list_result_joint$real_critical_value_joint

    # for development debugging
    # if (loop %% 100 == 0)
    # {
    #   if (!bool_mutePrint) {message(loop)}
    # }
  }

  # for development debugging
  # if (!bool_mutePrint) {
  #   print("average effective sample size")
  #   print(mean(effN.z.vec))
  #
  #   print("average correlation of Z1 and Z2")
  #   print(mean(covZ1Z2.vec))
  #
  #   print("mean z stat")
  #   print(mean(mean.jump.z.vec))
  #
  #   print("median z stat")
  #   print(mean(median.jump.z.vec))
  #
  #   print("max z stat")
  #   print(mean(max.jump.z.vec))
  #
  #
  #   print("naive rejection rate")
  #   print(naive.num.reject/ns)
  #   print("bonferroni rejection rate")
  #   print(bonfe.num.reject/ns)
  #   print("joint rejection rate")
  #   if (bool_max_test == FALSE) {joint.num.reject <- joint.num.reject[1,1]}
  #   print(joint.num.reject/ns)
  # }

  summary <- list(vec.statMaxJoint = vec.statMaxJoint,
                  chi.stat.vec = chi.stat.vec,
                  effN.mean = mean(effN.z.vec),
                  effNx.mean = mean(effN.x.vec),
                  cor.z = mean(covZ1Z2.vec),
                  naive = naive.num.reject/ns,
                  bonfe = bonfe.num.reject/ns,
                  joint = joint.num.reject/ns,
                  rwolf = rwolf.num.reject/ns,
                  maxCritical = vec.criticalValue)
  return(summary)
}


