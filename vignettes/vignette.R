## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(rdtest)
library(haven)
set.seed(1)
N <- 1000
Z <- data.frame(var1 = rnorm(N),
                var2 = rnorm(N),
                var3 = rnorm(N),
                var4 = rnorm(N),
                var5 = rnorm(N))
vec_X <- rnorm(N)

## ----basic usage--------------------------------------------------------------
res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X)
ans <- summary(res_rd)

## ----Not recommended----------------------------------------------------------
res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X, bool_L2_std = FALSE)
ans <- summary(res_rd)

## ----other statistics---------------------------------------------------------
res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X, bool_max_test = TRUE, bool_equivalence = TRUE)
ans <- summary(res_rd)

## ----1947 change--------------------------------------------------------------
data_path <- system.file("extdata", "018_births.dta", package = "rdtest")
df_018 <- haven::read_stata(data_path)
df_018_1 <- data.frame(df_018[df_018$r == 1,c("births","still","s1mob1")])
res_rd <- rdtest::rdtest(Z=df_018_1[,c("births","still")],vec_X = df_018_1$s1mob1)
ans <- summary(res_rd)

res_rd <- rdtest::rdtest(Z=df_018_1[,c("births","still")],vec_X = df_018_1$s1mob1,bool_max_test = TRUE,bool_equivalence = TRUE)
ans <- summary(res_rd)

## ----1972 change--------------------------------------------------------------
df_018 <- haven::read_stata(data_path)
df_018_2 <- data.frame(df_018[df_018$r == 2,c("births","still","s2mob1")])
vec_nan <- is.nan(df_018_2[,1]) + is.na(df_018_2[,1])
vec_nan <- vec_nan + is.nan(df_018_2[,2]) + is.na(df_018_2[,2])
vec_nan <- vec_nan + is.nan(df_018_2[,3]) + is.na(df_018_2[,3])
df_018_2 <- df_018_2[vec_nan == 0,]

res_rd <- rdtest::rdtest(Z=df_018_2[,c("births","still")],vec_X = df_018_2$s2mob1)
ans <- summary(res_rd)

res_rd <- rdtest::rdtest(Z=df_018_2[,c("births","still")],vec_X = df_018_2$s2mob1,bool_max_test = TRUE,bool_equivalence = TRUE)
ans <- summary(res_rd)

## -----------------------------------------------------------------------------
data_path <- system.file("extdata", "165_analysis.dta", package = "rdtest")
df_165 <- haven::read_stata(data_path)
list_var <- c("nb_registered_R1",
              "nb_candidates_R1",
              "prop_registered_turnout_R1",
              "distance_voteshare_cand12_R1",
              "predicted_assignment")
df_165 <- data.frame(df_165[,c(list_var,"running")])
vec_nan <- is.nan(df_165[,1]) + is.na(df_165[,1])
for (i in seq(2,6)) {
  vec_nan <- vec_nan + is.nan(df_165[,i]) + is.na(df_165[,i])
}
print(sum((vec_nan > 0)))
df_165 <- df_165[vec_nan == 0,]

res_rd <- rdtest::rdtest(Z=df_165[,list_var],vec_X = df_165$running)
ans <- summary(res_rd)

res_rd <- rdtest::rdtest(Z=df_165[,list_var],vec_X = df_165$running,bool_max_test = TRUE,bool_equivalence = TRUE)
ans <- summary(res_rd)

