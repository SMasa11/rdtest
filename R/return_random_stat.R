#'
#' return_random_stat
#'
#' an internal function for critical value computation
#'
#' @param int_num_simul Integer, number of simulation draws.
#' @param int_dim Integer, dimension to generate the draws.
#' @param mat_cor Matrix, correlation matrix of the target random variable.
#'

return_random_stat <- function(int_num_simul,
                               int_dim,
                               mat_cor)
{
    vec_random_stat <- rep(0,int_num_simul)
    #Generates the same random normal vector with the seed of 373568, restoring the currently loading seed
    mat_random_normal <-
      matrix(R.utils::withSeed({
        stats::rnorm((int_dim)*int_num_simul,0,1)},seed=373568)
        ,int_dim,int_num_simul)
    U <- svd(mat_cor)$u
    V <- svd(mat_cor)$v
    D <- diag(sqrt(svd(mat_cor)$d))
    mat_cor_sqrt <- U %*% D %*% t(V)
    mat_random_draw_stat <- mat_cor_sqrt %*% mat_random_normal
    return(mat_random_draw_stat)
}
