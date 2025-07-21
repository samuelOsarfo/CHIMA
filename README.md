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
  `get_active_med(y, x, M, COV.S=NULL, pval.adjust='HDMT', d=NULL, r=1, k=1)`

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

print("Test statistics for selected mediators:")
#> [1] "Test statistics for selected mediators:"
print(ao_result$ts)
#>  [1] -3.70761795  2.66703840  2.99542148 -3.71545035 -3.96376546  3.08844731
#>  [7] -5.97478920 -5.10032391 -4.29625008 -1.36792715 -0.09214669  1.56392316
#> [13] -0.92999742 -0.62393201 -0.29733202 -2.70895683 -2.93094953 -2.51810283
#> [19]  0.64989152 -0.48532488  1.04774762  0.95189802 -2.38632440  0.08716505
#> [25] -0.29177485 -0.51628685  0.53811427  0.84713618 -1.85332604  1.08590270
#> [31]  0.61327656  0.65272364  1.46638919  0.60566551 -1.08593595  1.56084125
#> [37]  0.16228273  0.65083067

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 2.092180e-04 7.652293e-03 2.740659e-03 2.028421e-04 7.377674e-05
#>  [6] 2.012054e-03 2.303879e-09 3.390727e-07 1.737116e-05 1.713349e-01
#> [11] 9.265815e-01 1.178356e-01 3.523724e-01 5.326722e-01 7.662130e-01
#> [16] 6.749512e-03 3.379277e-03 1.179889e-02 5.157623e-01 6.274459e-01
#> [21] 2.947549e-01 3.411487e-01 1.701773e-02 9.305403e-01 7.704588e-01
#> [26] 6.056541e-01 5.904982e-01 3.969192e-01 6.383564e-02 2.775220e-01
#> [31] 5.396935e-01 5.139345e-01 1.425423e-01 5.447369e-01 2.775073e-01
#> [36] 1.185612e-01 8.710832e-01 5.151558e-01
```

# Identifying Active Mediators

``` r
# Using HDMT For FDR control
active_mediators<- get_active_med(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  03:10:07 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:10:08 PM
#> Step 3: Joint Significance Testing   -----  03:10:09 PM
#> Complete!!   03:10:09 PM
print(active_mediators)
#> [1] 1 2 3 4 5 6 7 8


# Using Bonferroni correction
active_mediators_Bonferroni <- get_active_med(y, x, M, pval.adjust='bonferroni')
#> Step 1: Ridge-HOLP Screening   -----  03:10:09 PM
#> Step 2: Approximate Orthogonalization Estimates   -----  03:10:09 PM
#> Step 3: Joint Significance Testing   -----  03:10:10 PM
#> Complete!!   03:10:10 PM
print(active_mediators_Bonferroni)
#> [1] 1 4 5 7 8
```

# Reference

- Wang, X., & Leng, C. (2016). High dimensional ordinary least squares
  projection for screening variables. *Journal of the Royal Statistical
  Society Series B: Statistical Methodology*, 78(3), 589–611.
- Battey, H. S., & Reid, N. (2023). On inference in high-dimensional
  regression. *Journal of the Royal Statistical Society Series B:
  Statistical Methodology*, 85(1), 149–175.
