#'
#' Main function of the joint RD tests
#'
#' rdtest is ...
#'
#' @param df_Z A data frame of pre-determined covariates.
#' @param vec_X A vector of the running variable X.
#' @param bool_balance_only A boolean option, if TRUE, then balance test of z
#'   runs
#' @param int_J An integer of nearest neighbor for the variance estimation, the
#'   default value is 3 and the value should be <= 3.
#' @param real_cutoff A scalar real of cutoff value, default = 0.
#' @param bool_max_test A boolean option, if True, then max test runs instead of
#'   L2 test.
#'
#' @return A list of list_result containing test results.
#'
#'
#' @export
#'

rdtest <- function(df_Z,
                   vec_X,
                   bool_balance_only,
                   int_J = 3,
                   real_cutoff = 0,
                   bool_max_test = TRUE)
{

  #---
  # ERROR HANDLING
  #! CHECK IF int.dimZ > 0
  if (length(df_Z) == 0) {stop("Error: df_Z must not be empty.")}
  if (length(vec_X) == 0) {stop("Error: df_X must not be empty.")}
  if (!(length(vec_X) == length(t(vec_X)))) {stop("Error: vec_X should be a vector of a single variable.")}
  if (!is.data.frame(df_Z)) {stop("Error: df_Z must be a data.frame of pre-determined covariates.")}
  if (max(is.nan(vec_X))) {stop("Error: vec_X should not contain NaN.")}
  for (i in length(df_Z)) {if (max(is.nan(df_Z[,i]))) {stop("Error: df_Z should not contain NaN.")}}
  #! CHECK IF int_J is integer >= 1
  if (length(int_J) == 0) {stop("Error: int_J must not be empty")}
  if ((int_J %% 1)) {stop("Error: int_J need to be an integer.")}
  #! CHECK IF real_cutoff is scalar real value
  if (length(real_cutoff) == 0) {stop("Error: real_cutoff must not be empty")}
  if (length(real_cutoff) > 1) {stop("Error: real_cutoff must be a scalar")}
  #! CHECK IF df_Z is data.frame
  #! CHECK IF vec_X is a vector of a single variable
  if (!(length(vec_X) == length(df_Z[,1]))) {stop("Error: Number of observations should match for vec_X and df_Z.")}
  #! CHECK IF bool.BalanceOnly is boolean
  if (!(bool_balance_only == 0) & !(bool_balance_only == 1)) {stop("Error: bool.BalanceOnly need to be either 0 or 1")}

  #---
  # Retrieve the colnames of df_Z
  vec_col_names_Z <- colnames(df_Z)
  int_dim_Z <- length(vec_col_names_Z)
  colnames(df_Z) <- NULL

  #---
  # Normalize vec_X to have 0 at the cutoff
  vec_X <- vec_X - real_cutoff

  # Load the dataset
  # as the colnames removed, each column of df_Z is labeled as vec_Z.1, vec_Z.2 ,...
  df_data <- data.frame(vec_Z = df_Z,vec_X = vec_X)

  if (bool_balance_only)
  {
    list_result <- return_result_balance(df_data = df_data,int_dim_Z=int_dim_Z,int_J = int_J,bool_max_test=bool_max_test)
  } else {
    list_result <- return_result_joint(df_data = df_data,int_dim_Z = int_dim_Z,int_J = int_J,bool_max_test=bool_max_test)
  }

  return(list_result)
}
