Documentation For hdjmtTest
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

The `hdjmtTest` package performs high-dimensional mediator screening and
testing in the presence of high correlation among mediators using
Ridge-HOLP and Approximate Orthogonalization methods. In addition to
accounting for correlation, the method is designed to detect mediators
that exhibit joint effects. This vignette demonstrates how to use the
primary functions with a simulated dataset, `ExampleData`.

# Structure of ExampleData

The example data provided with the package contains:

- information on `n = 200` samples and `p = 2000` potential active
  mediators.
- **Variables**:
  - `M`: A 200x2000 matrix of mediators generated using the compound
    symmetry covariance structure to introduce high correlation among
    the mediators, with $\rho= 0.8$.
  - `x`: A vector of length 200 representing exposures.
  - `y`: A vector of length 200 representing outcomes.
  - `alp_vec`: A parameter vector of length 200 that relates the
    exposure variable to the mediators. The first 8 values are
    `-0.5127127, 0.6597036, 0.6756640, 0.5235137, 0.9305369, -0.9827865, -0.8941141, -0.9230220`
    with the rest being zeros.
  - `beta_vec`: A parameter vector of length 200 that relates the
    mediators to the outcome variable. The first 12 values are
    `-0.8033093, 0.9360975, 0.8185305, -0.7951502, -0.8783739, 0.8940459, -0.6911509, -0.8524771, -0.6812097, -0.8285034, -0.5986530, -0.9639383`
    with the rest being zeros.
- **Active Mediators**: given the non-zero values in `alp_vec` and
  `beta_vec`, only the first 8 mediators are truly active.

# How to install package from github

``` r
#library(devtools)

#devtools::install_github('samuelOsarfo/hdjmtTest',force=TRUE)
```

If package “qvalue” is not found, please first install “qvalue” package
through Bioconductor:
<https://www.bioconductor.org/packages/release/bioc/html/qvalue.html>

``` r
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("qvalue")

#devtools::install_github('samuelOsarfo/hdjmtTest',force=TRUE)
```

# Summary of Functions

## `medsc_holp`

- **Description**: Performs mediator screening using the Ridge-HOLP
  method.
- **Usage**: `medsc_holp(y, x, M, d = NULL, r = 1)`
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `M`: Matrix of mediators.
  - `COV.S`: a `data.frame` or `matrix` of covariates (optional).
  - `d`: Desired number of mediators to select (optional). Default value
    is $d= n/\log(n)$.
  - `r`: Ridge-HOLP penalty parameter (default is 1).

## `app_orth`

- **Description**: Fits the Approximate Orthogonalization model to test
  the significance of mediators.
- **Usage**:  
  `app_orth(y, x, chosen_M, COV.S=NULL, k = 1)`  
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `chosen_M`: Matrix of mediators selected during screening.
  - `COV.S`: a `data.frame` or `matrix` of covariates (optional).
  - `k`: Scalar used to compute projection directions (default is 1).

## `get_active_med.hima`, `get_active_med.mod`, `get_active_med.hdmt`

- **Description**: These are wrapper functions that combine screening
  and testing to identify active mediators using different
  joint-significance test variants.
- **Usage**:  
  `get_active_med.<variant>(y, x, M, COV.S=NULL, pval.adjust='HDMT', d=NULL, r=1, k=1)`  
- **Variants**:
  - `get_active_med.hima`: Uses the `null_estimation` function from the
    original HIMA repository.
  - `get_active_med.mod`: Uses a modified version of the
    `null_estimation` function.
  - `get_active_med.hdmt`: Uses the implementation in the HDMT package.
- **Arguments**:
  - `y`: Outcome vector.
  - `x`: Exposure vector.
  - `M`: Matrix of mediators.
  - `COV.S`: A `data.frame` or `matrix` of covariates (optional).
  - `pval.adjust`: Specifies which method to use for controlling FWER or
    FDR in the joint significance testing. Either `HDMT` (default) or
    `Bonferroni`.
  - `d`: Desired number of mediators to select (optional). Default is
    $d = n / \log(n)$.
  - `r`: Ridge penalty for HOLP (default = 1).
  - `k`: Scalar for computing projection directions for AO (default =
    1).

``` r
# Load the ExampleData
library(hdjmtTest)
data(ExampleData)

# Extract the components
y <- ExampleData$y
x <- ExampleData$x
M <- ExampleData$M
```

# Mediator Screening using Ridge-HOLP

``` r
# Perform Ridge-HOLP screening
chosen_ind <- medsc_holp(y, x, M)

print("Indexes of selected mediators:")
#> [1] "Indexes of selected mediators:"
print(chosen_ind)
#>  [1]    1    2    3    4    5    6    7    8    9   13   14   15   16  122  143
#> [16]  158  225  231  335  508  564  617  621  735  770  841  874  913 1121 1132
#> [31] 1141 1198 1199 1261 1283 1312 1341 1713
```

# Fitting Approximate Orthogonalization (AO)

``` r
# Apply AO on the chosen mediators
chosen_med <- M[, medsc_holp(y, x, M)]
ao_result <- app_orth(y, x, chosen_med)

print("Test statistics for selected mediators:")
#> [1] "Test statistics for selected mediators:"
print(ao_result$ts)
#>  [1] -3.72024578  2.67612211  3.00562363 -3.72810486 -3.97726571  3.09896630
#>  [7] -5.99513882 -5.11769518 -4.31088274 -1.37258620 -0.09246053  1.56924975
#> [13] -0.93316491 -0.62605707 -0.29834471 -2.71818331 -2.94093210 -2.52667928
#> [19]  0.65210500 -0.48697785  1.05131616  0.95514011 -2.39445202  0.08746192
#> [25] -0.29276861 -0.51804528  0.53994704  0.85002145 -1.85963831  1.08960119
#> [31]  0.61536533  0.65494676  1.47138359  0.60772835 -1.08963456  1.56615734
#> [37]  0.16283545  0.65304735

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 1.990290e-04 7.447950e-03 2.650368e-03 1.929251e-04 6.971225e-05
#>  [6] 1.941971e-03 2.033117e-09 3.092919e-07 1.626041e-05 1.698810e-01
#> [11] 9.263321e-01 1.165898e-01 3.507348e-01 5.312775e-01 7.654401e-01
#> [16] 6.564146e-03 3.272263e-03 1.151466e-02 5.143334e-01 6.262740e-01
#> [21] 2.931134e-01 3.395068e-01 1.664522e-02 9.303044e-01 7.696990e-01
#> [26] 6.044267e-01 5.892336e-01 3.953132e-01 6.293672e-02 2.758889e-01
#> [31] 5.383135e-01 5.125020e-01 1.411874e-01 5.433677e-01 2.758742e-01
#> [36] 1.173118e-01 8.706480e-01 5.137258e-01
```

# Identifying Active Mediators

``` r
# Using HDMT (HIMA null_estimation)
active_mediators.hima <- get_active_med.hima(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  07:34:08 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  07:34:08 PM
#> Step 3: Joint Significance Testing   -----  07:34:10 PM
#> Complete!!   07:34:13 PM
print(active_mediators.hima)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT (modified null_estimation)
active_mediators.mod <- get_active_med.mod(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  07:34:13 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  07:34:13 PM
#> Step 3: Joint Significance Testing   -----  07:34:14 PM
#> Complete!!   07:34:14 PM
print(active_mediators.mod)
#> [1] 1 2 3 4 5 6 7 8

# Using HDMT package implementation
active_mediators.hdmt <- get_active_med.hdmt(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  07:34:14 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  07:34:14 PM
#> Step 3: Joint Significance Testing   -----  07:34:15 PM
#> Complete!!   07:34:15 PM
print(active_mediators.hdmt)
#> [1] 1 2 3 4 5 6 7 8

# Using Bonferroni correction
active_mediators_Bonferroni <- get_active_med.hima(y, x, M, pval.adjust='bonferroni')
#> Step 1: Ridge-HOLP Screening   -----  07:34:15 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  07:34:15 PM
#> Step 3: Joint Significance Testing   -----  07:34:16 PM
#> Complete!!   07:34:16 PM
print(active_mediators_Bonferroni)
#> [1] 1 4 5 7 8
```

# Competing packages

``` r
## Install the HIMA package
# install.packages('HIMA')

suppressMessages(library(HIMA))
#> Warning: package 'HIMA' was built under R version 4.3.3
#> Warning: package 'ncvreg' was built under R version 4.3.3
#> Warning: package 'glmnet' was built under R version 4.3.3
#> Warning in .recacheSubclasses(def@className, def, env): undefined subclass
#> "ndiMatrix" of class "replValueSp"; definition not updated

HIMA::hima_dblasso(x, M, y)
#> Step 1: Sure Independent Screening ...  (7:34:18 PM)
#> Step 2: De-biased Lasso Estimates ...   (7:34:18 PM)
#> Step 3: Joint significance test ...     (7:34:25 PM)
#> Done!     (7:34:25 PM)
#>   Index  alpha_hat   alpha_se   beta_hat   beta_se        IDE      rimp
#> 1     1 -0.4989049 0.06159061 -0.7496780 0.2385597  0.3740180  7.329343
#> 2     3  0.6351005 0.05489418  1.2367930 0.2811790  0.7854878 15.392600
#> 3     4  0.5432349 0.05966641 -0.8091499 0.2513469 -0.4395585  8.613689
#> 4     5  0.7310222 0.04849277 -1.0239189 0.2910554 -0.7485075 14.667924
#> 5     6 -0.7606745 0.04613191  0.9481523 0.2956501 -0.7212353 14.133492
#> 6     7 -0.7180654 0.04946084 -1.4278684 0.2765547  1.0253029 20.092072
#> 7     8 -0.7294041 0.04861566 -1.3832008 0.2858711  1.0089124 19.770879
#>           pmax
#> 1 1.675018e-03
#> 2 1.089531e-05
#> 3 1.285240e-03
#> 4 4.348918e-04
#> 5 1.341231e-03
#> 6 2.429447e-07
#> 7 1.307913e-06
```

Out of the 8 active mediators,
`M[, 1], M[,3], M[,4], M[,5], M[,6], M[,7], M[,8]` were identified as
active by HIMA (specifically HIMA2).

# Reference

- Wang, X., & Leng, C. (2016). High dimensional ordinary least squares
  projection for screening variables. *Journal of the Royal Statistical
  Society Series B: Statistical Methodology*, 78(3), 589–611.
- Battey, H. S., & Reid, N. (2023). On inference in high-dimensional
  regression. *Journal of the Royal Statistical Society Series B:
  Statistical Methodology*, 85(1), 149–175.
