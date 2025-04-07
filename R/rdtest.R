# Copyright (c) 2023 Masayuki Sawada
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

#'
#' rdtest: A joint diagnosis for regression discontinuity designs.
#'
#' The rdtest package offers a joint test of commonly used diagnostic tests for
#' regression discontinuity designs.
#'
#' @param Z Data.frame or Vector, pre-determined covariates.
#' @param vec_X Vector, running variable X.
#' @param bool_joint Boolean, TRUE if joint test instead of balance test only,
#'   default is TRUE.
#' @param int_J Integer, the number or nearest neighbor for the variance
#'   estimation, default is 3.
#' @param real_cutoff Real scalar, cutoff value, default = 0.
#' @param bool_L2_std Boolean, use of standardized Wald test, default is TRUE.
#' @param bool_max_test Boolean, use of max test, instead of Wald tests, default
#'   is FALSE.
#' @param bool_equivalence Boolean, for returning equivalence test.
#' @param real_epsilon Numeric value, for equivalence test. Default NULL.
#'
#' @examples
#' # Prepare a mock dataset
#' library(rdtest)
#' set.seed(1)
#' N <- 1000
#' Z <- data.frame(
#'  var1 = rnorm(N),
#'  var2 = rnorm(N),
#'  var3 = rnorm(N),
#'  var4 = rnorm(N),
#'  var5 = rnorm(N))
#'  vec_X <- rnorm(N)
#' # Returning result with sWald test statistic
#' res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X)
#' ans <- summary(res_rd)
#' # Returning result with max test statistic
#' res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X, bool_max_test = TRUE)
#' ans <- summary(res_rd)
#' @export
#'

rdtest <- function(Z,
                   vec_X,
                   bool_joint = TRUE,
                   int_J = 3,
                   real_cutoff = 0,
                   bool_max_test = FALSE,
                   bool_L2_std = TRUE,
                   bool_equivalence = FALSE,
                   real_epsilon = NULL)
{

  #! CHECK IF int.dimZ > 0
  if (length(Z) == 0) {stop("Z must not be empty")}
  if (length(vec_X) == 0) {stop("vec_X must not be empty")}

  if (!(length(vec_X) == length(t(vec_X)))) {
    stop("vec_X should be a vector of a single variable")}
  if (!is.data.frame(Z)) {
    if (!is.numeric(Z)) {
      stop("Z must be either a data.frame
            or a vector of pre-determined covariates")
    } else {
        Z = data.frame(vec_Z = Z)
    }
  }
  if (max(is.nan(vec_X)) | max(is.na(vec_X))) {
    stop("vec_X should not contain NaN or NA")}
  for (i in length(Z)) {
    if (max(is.nan(Z[,i])) | max(is.na(Z[,i]))) {
      stop("Z should not contain NaN or NA")}
    }

  #! CHECK IF real.cutoff is scalar real value
  if (length(real_cutoff) == 0) {stop("real_cutoff must not be empty")}
  if (length(real_cutoff) > 1) {stop("real_cutoff must be a scalar")}

  #! CHECK IF X is a vector of a single variable
  if (!(length(vec_X) == length(Z[,1]))) {
    stop("Number of observations should match for X and Z.")}

  #! CHECK IF bool_joint is boolean
  if (!(bool_joint == 0) & !(bool_joint == 1)) {
    stop("bool_joint need to be either 0 or 1")}

  # Retrieve the colnames of Z
  vec_col_name_Z <- colnames(Z)
  int_dim_Z <- length(vec_col_name_Z)
  colnames(Z) <- NULL

  # Normalize X to have 0 at the cutoff
  vec_X <- vec_X - real_cutoff

  # making Z a data.frame by repeating vector, but keeping int_dim_Z = 1
  if (int_dim_Z == 1) {
    # the second argument will be discarded
    df_data <- data.frame(vec_Z.1 = Z, vec_X = vec_X)
  } else {
    # Load the dataset
    # as the colnames removed, Z is labeled as vec.Z.1, vec.Z.2 ,...
    df_data <- data.frame(vec_Z = Z,vec_X = vec_X)
  }

  if (bool_equivalence == TRUE) {
    # constructing the equivalence CI
    # increase real_eps_candid from 1 by multiplying 2 until list_result_CI$bool_equivalence_reject == TRUE
    real_eps_candid <- 1
    while (TRUE) {
      list_result_CI <- return_result_joint(
        df_data = df_data,
        int_dim_Z = int_dim_Z,
        int_J = int_J,
        bool_max_test = bool_max_test,
        bool_L2_std = bool_L2_std,
        bool_joint = bool_joint,
        bool_equivalence = bool_equivalence,
        real_epsilon = real_eps_candid)
      if (list_result_CI$bool_equivalence_reject == TRUE) {
        break
      }
      real_eps_candid <- real_eps_candid * 2
    }
    # use binary search to find the smallest epsilon that rejects the null
    real_eps_upper <- real_eps_candid
    real_eps_lower <- 0
    while (real_eps_upper - real_eps_lower > 1e-6) {
      real_eps_candid <- (real_eps_upper + real_eps_lower) / 2
      list_result_CI <- return_result_joint(
        df_data = df_data,
        int_dim_Z = int_dim_Z,
        int_J = int_J,
        bool_max_test = bool_max_test,
        bool_L2_std = bool_L2_std,
        bool_joint = bool_joint,
        bool_equivalence = bool_equivalence,
        real_epsilon = real_eps_candid)
      if (list_result_CI$bool_equivalence_reject == TRUE) {
        real_eps_upper <- real_eps_candid
      } else {
        real_eps_lower <- real_eps_candid
      }
    }
    real_eps_CI <- c(-real_eps_upper, real_eps_upper)
  } else {
    real_eps_CI <- NA
  }

  list_result <- return_result_joint(
    df_data = df_data,
    int_dim_Z = int_dim_Z,
    int_J = int_J,
    bool_max_test = bool_max_test,
    bool_L2_std = bool_L2_std,
    bool_joint = bool_joint,
    bool_equivalence = bool_equivalence,
    real_epsilon = real_epsilon)
  list_result$call <- match.call()
  list_result$bool_max_test <- bool_max_test
  list_result$bool_L2_std <- bool_L2_std
  list_result$bool_joint <- bool_joint
  list_result$vec_col_name_Z <- vec_col_name_Z
  list_result$int_dim_Z <- int_dim_Z
  list_result$bool_equivalence <- bool_equivalence
  list_result$real_epsilon <- real_epsilon
  list_result$real_eps_CI <- real_eps_CI

  class(list_result) <- "rdtest"

  return(list_result)
}



#'
#' Summary function of the main caller.
#'
#' @param object Object rdtest, result returned from rdtest function.
#' @param ... Other options.
#' @export
summary.rdtest <- function (object, ...)
{
  ans <- object

  cat("Call: ")
  print(ans$call)

  ans$tstats <- matrix(NA, 0L, 3L,
                       dimnames = list(NULL,
                                       c("Variable",
                                         "t-stat",
                                         "(Naive) p-value")))

  if (ans$bool_joint) {
    ans$tstats <-
      rbind(ans$tstats,
      cbind(Variable = "density",
            `t-stat` = round(ans$real_tstat_X,3),
            `(Naive) p-value` =
              round((1-stats::pnorm(abs(ans$real_tstat_X))),3)))
  }
  for (i in seq(1:ans$int_dim_Z)) {
    ans$tstats <-
      rbind(ans$tstats,
            cbind(Variable = ans$vec_col_name_Z[i],
                  `t-stat` = round(ans$vec_tstat_Z_raw[i],3),
                  `(Naive) p-value`
                  = round((1-stats::pnorm(abs(ans$vec_tstat_Z_raw[i]))),3)))
  }

  cat("\n Naive test results:")
  print(knitr::kable(ans$tstats,"rst"))
  cat("\n")

  if (ans$bool_joint) {
    cat("\n Joint test result: covariates and density.")
  } else {
    cat("\n Joint test result: covariates.")
  }

  ans$joint_result <- matrix(NA, 0L, 4L,
                       dimnames = list(NULL,
                                       c("Test Type",
                                         "Test Statistic",
                                         "Critical Value",
                                         "Joint p-value")))

  if (ans$bool_max_test) {
    test_type = "Max test"
  } else {
    if (ans$bool_L2_std) {
      test_type = "standardized Wald test"
    } else {
      test_type = "(non-standardized) Wald test"
    }
  }
  ans$joint_result <-
    rbind(ans$joint_result,
          cbind(
            `Test Type` = test_type,
            `Test Statistic` = round(ans$real_stat_joint,3),
            `Critical Value` = round(ans$real_critical_value_joint,3),
            `Joint p-alue` = round(ans$real_pvalue,3)
          ))
  print(knitr::kable(ans$joint_result,"rst"))

  if (ans$bool_equivalence) {
    if (is.null(ans$real_epsilon) == FALSE) {
      print(paste0("Equivalence test at epsilon = " ,ans$real_epsilon,": rejection = ",ans$bool_equivalence_reject))
    }
    print(paste0(
      "Equivalence test can reject the null hypothesis of min(stat) <= ",
      round(ans$real_eps_CI[1],3),
      " or max(stat) >= ",
      round(ans$real_eps_CI[2],3),
      " at the 5% level.")
    )
  }

  return(ans)
}
