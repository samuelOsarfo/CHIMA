Documentation For CHIMA
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

The `CHIMA` package performs high-dimensional mediator screening and
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

#devtools::install_github('samuelOsarfo/CHIMA',force=TRUE)
```

If package “qvalue” is not found, please first install “qvalue” package
through Bioconductor:
<https://www.bioconductor.org/packages/release/bioc/html/qvalue.html>

``` r
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("qvalue")

#devtools::install_github('samuelOsarfo/CHIMA',force=TRUE)
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

## `get_active_med`

- **Description**: This is wrapper function that combine screening and
  testing to identify active mediators.

- **Usage**:  
  `get_active_med(y, x, M, COV.S=NULL, pval.adjust='HDMT', d=NULL, r=1, k=1, alpha=0.05)`

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
  - `alpha`: Target FDR level (default =0.05)

``` r
# Load the ExampleData
library(CHIMA)
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

print("beta_j's estimates for selected mediators:")
#> [1] "beta_j's estimates for selected mediators:"
print(ao_result$bhat)
#>  [1] -0.85454891  0.72798593  0.91841399 -0.88883970 -1.23070731  0.99816884
#>  [7] -1.68936358 -1.52036730 -1.01303848 -0.41023057 -0.02754499  0.36738506
#> [13] -0.27628426 -0.12108618 -0.06281622 -0.57080517 -0.61709430 -0.54748920
#> [19]  0.14286997 -0.10754716  0.21732091  0.19888873 -0.48871179  0.01647673
#> [25] -0.05994742 -0.10851191  0.11491766  0.17757867 -0.37963078  0.24259622
#> [31]  0.12993300  0.13889385  0.30378477  0.12870413 -0.24002365  0.34967104
#> [37]  0.03512818  0.13361500

print("Test statistics for selected mediators:")
#> [1] "Test statistics for selected mediators:"
print(ao_result$ts)
#>  [1] -3.69847699  2.66046294  2.98803640 -3.70629008 -3.95399298  3.08083289
#>  [7] -5.96005863 -5.08774929 -4.28565787 -1.36455459 -0.09191950  1.56006738
#> [13] -0.92770456 -0.62239374 -0.29659896 -2.70227802 -2.92372341 -2.51189456
#> [19]  0.64828924 -0.48412833  1.04516445  0.94955116 -2.38044103  0.08695014
#> [25] -0.29105549 -0.51501397  0.53678757  0.84504760 -1.84875675  1.08322546
#> [31]  0.61176456  0.65111438  1.46277388  0.60417227 -1.08325863  1.55699307
#> [37]  0.16188263  0.64922608

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 2.168970e-04 7.803331e-03 2.807761e-03 2.103174e-04 7.685766e-05
#>  [6] 2.064225e-03 2.521475e-09 3.623379e-07 1.821992e-05 1.723931e-01
#> [11] 9.267620e-01 1.187440e-01 3.535608e-01 5.336830e-01 7.667727e-01
#> [16] 6.886615e-03 3.458719e-03 1.200850e-02 5.167979e-01 6.282948e-01
#> [21] 2.959470e-01 3.423404e-01 1.729193e-02 9.307111e-01 7.710089e-01
#> [26] 6.065433e-01 5.914144e-01 3.980843e-01 6.449295e-02 2.787084e-01
#> [31] 5.406935e-01 5.149727e-01 1.435293e-01 5.457291e-01 2.786936e-01
#> [36] 1.194721e-01 8.713983e-01 5.161923e-01
```

# Identifying Active Mediators

``` r
# Using HDMT For FDR control
active_mediators<- get_active_med(y, x, M)
#> Step 1: Ridge–HOLP Screening   -----  10:34:02  am
#> Step 2: Approximate Orthogonalization Estimates   -----  10:34:02 AM
#> Step 3: Joint Significance Testing   -----  10:34:02 AM
#> Complete!!   10:34:02 AM
print(active_mediators)
#>   Index  alpha_hat   beta_hat      P_value
#> 1     1 -0.6113066 -0.8545489 1.838289e-04
#> 2     2  0.8772111  0.7279859 7.134070e-03
#> 3     3  0.8735573  0.9184140 2.512810e-03
#> 4     4  0.6831304 -0.8888397 1.781345e-04
#> 5     5  1.1462583 -1.2307073 6.370095e-05
#> 6     6 -1.1805587  0.9981688 1.835465e-03
#> 7     7 -1.0269161 -1.6893636 1.666244e-09
#> 8     8 -1.1017281 -1.5203673 2.671939e-07
```

# Reference

- Wang, X., & Leng, C. (2016). High dimensional ordinary least squares
  projection for screening variables. *Journal of the Royal Statistical
  Society Series B: Statistical Methodology*, 78(3), 589–611.
- Battey, H. S., & Reid, N. (2023). On inference in high-dimensional
  regression. *Journal of the Royal Statistical Society Series B:
  Statistical Methodology*, 85(1), 149–175.
