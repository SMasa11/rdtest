#'
#' Main Test Function
#'
#' @param df_Z Data.frame, pre-determined covariates
#' @param vec_X Vector, running variable X
#' @param bool_joint Boolean, TRUE if joint test instead of balance test only.
#' @param int_J integer of nearest neighbor for the variance estimation <= 3
#' @param real_cutoff scalar real of cutoff value, default = 0
#' @param bool_max_test_V_inv Boolean, option for using inverse V for max test.
#' @param bool_L2_std Boolean, for using standardized t stat for L2 test.
#' @param bool_max_test Boolean,
#'
#' @export
#'

rdtest <- function(df_Z,
                   vec_X,
                   bool_joint = TRUE,
                   int_J = 3,
                   real_cutoff = 0,
                   bool_max_test = FALSE,
                   bool_max_test_V_inv = FALSE,
                   bool_L2_std = TRUE)
{

  #! CHECK IF int.dimZ > 0
  if (length(df_Z) == 0) {stop("df_Z must not be empty")}
  if (length(vec_X) == 0) {stop("vec_X must not be empty")}

  if (!(length(vec_X) == length(t(vec_X)))) {
    stop("vec_X should be a vector of a single variable")}
  if (!is.data.frame(df_Z)) {
    stop("df_Z must be a data.frame of pre-determined covariates")}
  if (max(is.nan(vec_X))) {
    stop("vec_X should not contain NaN")}
  for (i in length(df_Z)) {
    if (max(is.nan(df_Z[,i]))) {stop("df_Z should not contain NaN")}}

  #! CHECK IF int.J is integer >= 1
  if (length(int_J) == 0) {stop("int.J must not be empty")}
  if ((int_J %% 1)) {stop("int.J need to be an integer.")}
  #! CHECK IF real.cutoff is scalar real value
  if (length(real_cutoff) == 0) {stop("real.cutoff must not be empty")}
  if (length(real_cutoff) > 1) {stop("real.cutoff must be a scalar")}
  #! CHECK IF Z is data.frame
  #! CHECK IF X is a vector of a single variable
  if (!(length(vec_X) == length(df_Z[,1]))) {
    stop("Number of observations should match for X and Z.")}

  #! CHECK IF bool_joint is boolean
  if (!(bool_joint == 0) & !(bool_joint == 1)) {
    stop("bool_joint need to be either 0 or 1")}

  # Retrieve the colnames of Z
  vec_col_name_Z <- colnames(df_Z)
  int_dim_Z <- length(vec_col_name_Z)
  colnames(df_Z) <- NULL

  # Normalize X to have 0 at the cutoff
  vec_X <- vec_X - real_cutoff

  # Load the dataset
  # as the colnames removed, Z is labeled as vec.Z.1, vec.Z.2 ,...
  df_data <- data.frame(vec_Z = df_Z,vec_X = vec_X)

  list_result <- return_result_joint(
    df_data = df_data,
    int_dim_Z = int_dim_Z,
    int_J = int_J,
    bool_max_test = bool_max_test,
    bool_max_test_V_inv = bool_max_test_V_inv,
    bool_L2_std = bool_L2_std,
    bool_joint = bool_joint)

  return(list_result)
}
