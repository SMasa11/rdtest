#'
#' Main Test Function
#'
#' @param Z data frame of pre-determined covariates
#' @param X vector of the running variable X
#' @param bool.BalanceOnly Boolean TRUE if joint balance test of z should be run
#' @param int.J integer of nearest neighbor for the variance estimation <= 3
#' @param real.cutoff scalar real of cutoff value, default = 0
#'
#' @export
#'

rdtest <- function(Z,X,bool.BalanceOnly,int.J=3,real.cutoff=0,bool.maxTest=TRUE)
{

  #! CHECK IF int.dimZ > 0
  if (length(Z) == 0) {stop("Z must not be empty")}
  if (length(X) == 0) {stop("X must not be empty")}

  if (!(length(X) == length(t(X)))) {stop("X should be a vector of a single variable")}
  if (!is.data.frame(Z)) {stop("Z must be a data.frame of pre-determined covariates")}
  if (max(is.nan(X))) {stop("X should not contain NaN")}
  for (i in length(Z)) {if (max(is.nan(Z[,i]))) {stop("Z should not contain NaN")}}

  #! CHECK IF int.J is integer >= 1

  if (length(int.J) == 0) {stop("int.J must not be empty")}
  if ((int.J %% 1)) {stop("int.J need to be an integer.")}
  #! CHECK IF real.cutoff is scalar real value
  if (length(real.cutoff) == 0) {stop("real.cutoff must not be empty")}
  if (length(real.cutoff) > 1) {stop("real.cutoff must be a scalar")}
  #! CHECK IF Z is data.frame
  #! CHECK IF X is a vector of a single variable
  if (!(length(X) == length(Z[,1]))) {stop("Number of observations should match for X and Z.")}

  #! CHECK IF bool.BalanceOnly is boolean
  if (!(bool.BalanceOnly == 0) & !(bool.BalanceOnly == 1)) {stop("bool.BalanceOnly need to be either 0 or 1")}

  # Retrieve the colnames of Z
  vec.colNamesZ <- colnames(Z)
  int.dimZ <- length(vec.colNamesZ)
  colnames(Z) <- NULL

  # Normalize X to have 0 at the cutoff
  X <- X - real.cutoff

  # Load the dataset
  # as the colnames removed, Z is labeled as vec.Z.1, vec.Z.2 ,...
  df.data <- data.frame(vec.Z = Z,vec.X = X)

  if (bool.BalanceOnly)
  {
    list.result <- return_result_balance(df.data = df.data,int.dimZ=int.dimZ,int.J = int.J,bool.maxTest=bool.maxTest)
  } else {
    list.result <- return_results_joint(df.data = df.data,int.dimZ = int.dimZ,int.J = int.J,bool.maxTest=bool.maxTest)
  }

  return(list.result)
}
