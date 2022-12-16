library(rdtest)

vec_X = rnorm(100)
vec_Z = rnorm(100)
df_Z = data.frame(Z1 = rnorm(100), Z2 = rnorm(100))

testthat::test_that("empty Z is an error", {
  testthat::expect_error(rdtest::rdtest(Z = data.frame(), vec_X = vec_X),
                         regexp = "Z must not be empty")
})

testthat::test_that("empty X is an error", {
  testthat::expect_error(rdtest::rdtest(Z = vec_Z, vec_X = c()),
                         regexp = "vec_X must not be empty")
})

testthat::test_that("vec_X should not be multi-dimensional", {
  testthat::expect_error(rdtest::rdtest(Z = df_Z,
                                        vec_X = data.frame(X1=rnorm(100),X2=rnorm(100))),
                         regexp = "vec_X should be a vector of a single variable")
})

testthat::test_that("Z should be either df or a vector", {
  testthat::expect_error(rdtest::rdtest(Z = "test", vec_X = vec_X),
                         regexp = "Z must be either a data.frame
            or a vector of pre-determined covariates")
})

testthat::test_that("NaN in Z or X is detected", {
  testthat::expect_error(rdtest::rdtest(Z = c(vec_Z,NaN),vec_X=c(vec_X,1)),
                         regexp = "Z should not contain NaN or NA")
  testthat::expect_error(rdtest::rdtest(Z = c(vec_Z,0),vec_X=c(vec_X,NA)),
                         regexp = "vec_X should not contain NaN or NA")
  testthat::expect_error(rdtest::rdtest(Z = c(vec_Z,0),vec_X=c(vec_X,NaN)),
                         regexp = "vec_X should not contain NaN or NA")
})

testthat::test_that("real_cutoff is scalar real value", {
  testthat::expect_error(rdtest::rdtest(Z = df_Z, vec_X = vec_X,
                                        real_cutoff = c()),
                         regexp = "real_cutoff must not be empty")
  testthat::expect_error(rdtest::rdtest(Z = df_Z, vec_X = vec_X,
                                        real_cutoff = c(0,0)),
                         regexp = "real_cutoff must be a scalar")
})

testthat::test_that("Number of obs should match", {
  testthat::expect_error(rdtest::rdtest(Z = c(vec_Z,0), vec_X = vec_X),
                         regexp = "Number of observations should match for X and Z")
})

testthat::test_that("bool_joint is boolean", {
  testthat::expect_error(rdtest::rdtest(Z = df_Z, vec_X = vec_X,
                                        bool_joint = 2),
                         regexp = "bool_joint need to be either 0 or 1")
})


