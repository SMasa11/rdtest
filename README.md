
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rdtest

<!-- badges: start -->

[![R-CMD-check](https://github.com/SMasa11/rdtest/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SMasa11/rdtest/actions/workflows/R-CMD-check.yaml)
<!-- badges: end --> 

*rdtest* offers a joint diagnostic procedure for
the testable restrictions of regression discontinuity (RD) designs. This
package combines the existing packages of diagnostic testing (rddensity
and rdrobust) to offer a set of joint test statistics to evaluate the
validity of a given RD design. The theoretical background of this
package is provided in [Fusejima, Ishihara, and Sawada (2023,
arXiv)](https://arxiv.org/abs/2205.04345). This work is supported by
JSPS KAKENHI Grant Number JP20J20046 (Fusejima); JP20J00900 (Ishihara);
and 21K13269 (Sawada).

## Background

The sufficient condition for identification in RD design, local
randomization, has has two implications: continuous density function and
balancing of pre-treatment covariates. Violation of any of these
testable restriction raises aa concern on the sufficient condition.
(Nevertheless, these testable restrictions are not directly related to
identification itself. See [Ishihara and Sawada
(2023)](https://arxiv.org/abs/2009.07551) for further details.)

Both testable restrictions are implemented in R packages, for example,
[*rddensity*](https://rdpackages.github.io/rddensity/) for the former
and [*rdrobust*](https://rdpackages.github.io/rdrobust/) for the latter.
As pointed out in [Fusejima, et
al. (2023)](https://arxiv.org/abs/2205.04345), however, multiple tests
(median number of tests is 12) are run in practice without appropriate
size control for the most of studies. Consequently, these RD practices
suffer from massive size distortion due to *multiple testing problem*.
Among 60 papers studied, 35% of studies reject at least one of these
testable restrictions; however, rejection of each test separately is
approximately 5%. Hence, the RD designs in these studies appear to
maintain the null hypotheses; however, not a small fraction of studies
contain false rejection of the null hypotheses more frequently than it
should be. Our *rdtest* package is the resolution to this multiple
testing problem.

## Installation

You can install the development version of rdtest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SMasa11/rdtest")
```

## Usage

*rdtest* package minimally expect a pair of variables:

- Z : a vector or a data.frame of pre-determined covariates whose
  conditional expectation should be continous at the cutoff under local
  randomization.

- vec_X: a vector of running variable.

For example,

``` r
library(rdtest)
set.seed(1)
N <- 1000
Z <- data.frame(var1 = rnorm(N),
                var2 = rnorm(N),
                var3 = rnorm(N),
                var4 = rnorm(N),
                var5 = rnorm(N),
                var6 = rnorm(N),
                var7 = rnorm(N),
                var8 = rnorm(N),
                var9 = rnorm(N),
                var10 = rnorm(N))
vec_X <- rnorm(N)
```

In default, a version of Wald statistics (`sWald`, standardized Wald
statistics. See [Fusejima, et
al. (2023)](https://arxiv.org/abs/2205.04345) for its detail).

``` r
res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X)
ans <- summary(res_rd)
#> Call: rdtest::rdtest(Z = Z, vec_X = vec_X)
#> 
#>  Naive test results:
#> 
#> ========  ======  ===============
#> Variable  t-stat  (Naive) p-value
#> ========  ======  ===============
#> density   -0.699  0.242          
#> var1      0.469   0.32           
#> var2      -0.34   0.367          
#> var3      0.868   0.193          
#> var4      1.171   0.121          
#> var5      1.935   0.026          
#> var6      1.47    0.071          
#> var7      -1.175  0.12           
#> var8      0.627   0.265          
#> var9      -0.871  0.192          
#> var10     -0.017  0.493          
#> ========  ======  ===============
#> 
#> 
#>  Joint test result: covariates and density.
#> 
#> ===  ======================  ==============  ==============  =============
#> \    Test Type               Test Statistic  Critical Value  Joint p-value
#> ===  ======================  ==============  ==============  =============
#> 95%  standardized Wald test  11.434          20.639          0.398        
#> ===  ======================  ==============  ==============  =============
```

With `bool_max_test` option TRUE, `rdtest` returns max test result
instead of `sWald` test result.

``` r
res_rd <- rdtest::rdtest(Z = Z, vec_X = vec_X, bool_max_test = TRUE)
ans <- summary(res_rd)
#> Call: rdtest::rdtest(Z = Z, vec_X = vec_X, bool_max_test = TRUE)
#> 
#>  Naive test results:
#> 
#> ========  ======  ===============
#> Variable  t-stat  (Naive) p-value
#> ========  ======  ===============
#> density   -0.699  0.242          
#> var1      0.469   0.32           
#> var2      -0.34   0.367          
#> var3      0.868   0.193          
#> var4      1.171   0.121          
#> var5      1.935   0.026          
#> var6      1.47    0.071          
#> var7      -1.175  0.12           
#> var8      0.627   0.265          
#> var9      -0.871  0.192          
#> var10     -0.017  0.493          
#> ========  ======  ===============
#> 
#> 
#>  Joint test result: covariates and density.
#> 
#> ===  =========  ==============  ==============  =============
#> \    Test Type  Test Statistic  Critical Value  Joint p-value
#> ===  =========  ==============  ==============  =============
#> 95%  Max test   3.782           8.08            0.417        
#> ===  =========  ==============  ==============  =============
```

In the above example, var5 has naive p-value of 0.026, but neither
*sWald* test nor *max* test reject the joint null hypothesis of the
whole 10 covariates and the density satisfy the continuity condition as
designed.

## Which tests should be used?

Based on our simulation exercises, we recommend the following
rule-of-thumb:

- *max* test is recommended when a researcher has a particularly
  concerning covariate among up to 5 covariates. For example, a proxy
  variable for the outcome or a placebo test may be of interest.

- *sWald* test is recommended when a researcher has no prior concern on
  a specific covariate in violation of the continuity. In simulation,
  *sWald* test demonstrates a better performance than *max* test for 10
  and 25 covariates in its size control.

- Nevertheless, we recommend limiting one’s focus on a fewer set of
  covariates (at most 10) for either tests.
