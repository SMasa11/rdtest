#'
#' Main Test Function
#'
#' @param Z data frame of pre-determined covariates
#' @param X vector of the running variable X
#' @param bool_balance Boolean TRUE if joint balance test of z should be run
#' @param int_J integer of nearest neighbor for the variance estimation <= 3
#' @param real_cutoff scalar real of cutoff value, default = 0
#' @param bool_max_test Boolean,
#'
#' @export
#'

rdtest <- function(Z,
                   X,
                   bool_balance,
                   int_J = 3,
                   real_cutoff = 0,
                   bool_max_test = TRUE)
{

  #! CHECK IF int.dimZ > 0
  if (length(Z) == 0) {stop("Z must not be empty")}
  if (length(X) == 0) {stop("X must not be empty")}

  if (!(length(X) == length(t(X)))) {stop("X should be a vector of a single variable")}
  if (!is.data.frame(Z)) {stop("Z must be a data.frame of pre-determined covariates")}
  if (max(is.nan(X))) {stop("X should not contain NaN")}
  for (i in length(Z)) {if (max(is.nan(Z[,i]))) {stop("Z should not contain NaN")}}

  #! CHECK IF int.J is integer >= 1

  if (length(int_J) == 0) {stop("int.J must not be empty")}
  if ((int_J %% 1)) {stop("int.J need to be an integer.")}
  #! CHECK IF real.cutoff is scalar real value
  if (length(real_cutoff) == 0) {stop("real.cutoff must not be empty")}
  if (length(real_cutoff) > 1) {stop("real.cutoff must be a scalar")}
  #! CHECK IF Z is data.frame
  #! CHECK IF X is a vector of a single variable
  if (!(length(X) == length(Z[,1]))) {stop("Number of observations should match for X and Z.")}

  #! CHECK IF bool_balance is boolean
  if (!(bool_balance == 0) & !(bool_balance == 1)) {stop("bool_balance need to be either 0 or 1")}

  # Retrieve the colnames of Z
  vec_col_name_Z <- colnames(Z)
  int_dim_Z <- length(vec_col_name_Z)
  colnames(Z) <- NULL

  # Normalize X to have 0 at the cutoff
  X <- X - real_cutoff

  # Load the dataset
  # as the colnames removed, Z is labeled as vec.Z.1, vec.Z.2 ,...
  df_data <- data.frame(vec.Z = Z,vec.X = X)

  list_result <- return_results_joint(df_data = df_data,
                                      int_dim_Z = int_dim_Z,
                                      int_J = int_J,
                                      bool_max_test = bool_max_test,
                                      bool_balance = bool_balance)

  return(list_result)
}
