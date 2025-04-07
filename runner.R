# runner.R
devtools::document()
devtools::load_all()
devtools::build()
devtools::install()

set.seed(3863224)

Z <- data.frame(
    var1 = rnorm(1000),
    var2 = rnorm(1000),
    var3 = rnorm(1000),
    var4 = rnorm(1000),
    var5 = rnorm(1000))
vec_X <- rnorm(1000)

res <- rdtest::rdtest(
    Z = Z,
    vec_X = vec_X,
    bool_joint = FALSE,
    int_J = 3,
    real_cutoff = 0,
    bool_max_test = TRUE)
summary(res)

res <- rdtest::rdtest(
    Z = Z,
    vec_X = vec_X,
    bool_joint = TRUE,
    int_J = 3,
    real_cutoff = 0,
    bool_max_test = TRUE,
    bool_equivalence = TRUE,
    real_epsilon = 0.1)
summary(res)