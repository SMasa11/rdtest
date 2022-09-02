####
#'
#' return_random_stat
#'
#' @param int_num_simul Integer, number of simulation draws.
#' @param int_dim Integer, dimension to generate the draws.
#' @param mat_Cor Matrix, correlation matrix of the target random variable.
#'
####

return_random_stat <- function(int_num_simul,
                               int_dim,
                               mat_Cor)
{
    vec_random_stat <- rep(0,int_num_simul)
    #Generates the same random normal vector with the seed of 373568, restoring the currently loading seed
    mat_random_normal <-
      matrix(R.utils::withSeed({
        rnorm((int_dim)*int_num_simul,0,1)},seed=373568)
        ,int_dim,int_num_simul)
    U <- svd(mat_Cor)$u
    V <- svd(mat_Cor)$v
    D <- diag(sqrt(svd(mat_Cor)$d))
    mat_CorSqrt <- U %*% D %*% t(V)
    mat_randomDrawStats <- mat_CorSqrt %*% mat_random_normal
    return(mat_randomDrawStats)
}
