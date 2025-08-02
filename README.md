Documentation For CHIMA
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

The `CHIMA` package performs high-dimensional mediator screening and
testing in the presence of high correlation among mediators using
Ridge-HOLP and Approximate Orthogonalization methods. This vignette demonstrates how to use the
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
# Load ExampleData
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
#>  [1] -3.77233689  2.71359333  3.04770856 -3.78030601 -4.03295563  3.14235822
#>  [7] -6.07908313 -5.18935347 -4.37124399 -1.39180523 -0.09375517  1.59122248
#> [13] -0.94623115 -0.63482316 -0.30252215 -2.75624348 -2.98211121 -2.56205800
#> [19]  0.66123581 -0.49379655  1.06603676  0.96851404 -2.42797929  0.08868657
#> [25] -0.29686797 -0.52529898  0.54750741  0.86192350 -1.88567708  1.10485786
#> [31]  0.62398172  0.66411737  1.49198599  0.61623780 -1.10489169  1.58808678
#> [37]  0.16511548  0.66219136

print("P-values for selected mediators:")
#> [1] "P-values for selected mediators:"
print(ao_result$pval)
#>  [1] 1.617257e-04 6.655782e-03 2.305934e-03 1.566357e-04 5.507968e-05
#>  [6] 1.675929e-03 1.208717e-09 2.110255e-07 1.235406e-05 1.639814e-01
#> [11] 9.253036e-01 1.115595e-01 3.440307e-01 5.255438e-01 7.622541e-01
#> [16] 5.846946e-03 2.862680e-03 1.040539e-02 5.084611e-01 6.214499e-01
#> [21] 2.864070e-01 3.327877e-01 1.518321e-02 9.293310e-01 7.665673e-01
#> [26] 5.993754e-01 5.840302e-01 3.887296e-01 5.933848e-02 2.692212e-01
#> [31] 5.326396e-01 5.066152e-01 1.357028e-01 5.377376e-01 2.692065e-01
#> [36] 1.122667e-01 8.688531e-01 5.078486e-01
```

# Identifying Active Mediators

``` r
# Using HDMT For FDR control
active_mediators<- get_active_med(y, x, M)
#> Step 1: Ridge-HOLP Screening   -----  20:27:19  PM
#> Step 2: Approximate Orthogonalization Estimates   -----  20:27:19 PM
#> Step 3: Joint Significance Testing   -----  20:27:20 PM
#> Complete!!   20:27:20 PM
print(active_mediators)
#>   Index  alpha_hat   beta_hat      P_value
#> 1     1 -0.6113066 -0.8545489 1.923956e-04
#> 2     2  0.8772111  0.7279859 7.312360e-03
#> 3     3  0.8735573  0.9184140 2.590775e-03
#> 4     4  0.6831304 -0.8888397 1.864699e-04
#> 5     5  1.1462583 -1.2307073 6.708103e-05
#> 6     6 -1.1805587  0.9981688 1.895790e-03
#> 7     7 -1.0269161 -1.6893636 1.867633e-09
#> 8     8 -1.1017281 -1.5203673 2.905756e-07
```

# Competing Package: HIMA

``` r
library(HIMA)
#> HIMA version 2.3.0
#> To access full functionality of HIMA, please make sure this version is current.
#> 
#> Citation:
#>   1. Zhang H, Zheng Y, Zhang Z, Gao T, Joyce B, Yoon G, Zhang W, Schwartz J,
#>      Just A, Colicino E, Vokonas P, Zhao L, Lv J, Baccarelli A, Hou L, Liu L.
#>      Estimating and Testing High-dimensional Mediation Effects in Epigenetic Studies.
#>      Bioinformatics. 2016.
#>      PMID: 27357171; PMCID: PMC5048064.
#> 
#>   2. Zhang H, Zheng Y, Hou L, Zheng C, Liu L.
#>      Mediation Analysis for Survival Data with High-Dimensional Mediators.
#>      Bioinformatics. 2021.
#>      PMID: 34343267; PMCID: PMC8570823.
#> 
#>   3. Zhang H, Chen J, Feng Y, Wang C, Li H, Liu L.
#>      Mediation effect selection in high-dimensional and compositional microbiome data.
#>      Stat Med. 2021.
#>      PMID: 33205470; PMCID: PMC7855955.
#> 
#>   4. Perera C, Zhang H, Zheng Y, Hou L, Qu A, Zheng C, Xie K, Liu L.
#>      HIMA2: high-dimensional mediation analysis and its application in epigenome-wide DNA methylation data.
#>      BMC Bioinformatics. 2022.
#>      PMID: 35879655; PMCID: PMC9310002.
#> 
#>   5. Zhang H, Hong X, Zheng Y, Hou L, Zheng C, Wang X, Liu L.
#>      High-Dimensional Quantile Mediation Analysis with Application to a Birth Cohort Study of Mother-Newborn Pairs.
#>      Bioinformatics. 2024.
#>      PMID: 38290773; PMCID: PMC10873903.
#> 
#>   6. Bai X, Zheng Y, Hou L, Zheng C, Liu L, Zhang H.
#>      An Efficient Testing Procedure for High-dimensional Mediators with FDR Control.
#>      Statistics in Biosciences. 2024.
#> ************************************************************************************************************************

 HIMA::dblassoHIMA(x, M, y)
#> Step 1: Sure Independent Screening ...  (8:27:20 PM)
#> Step 2: De-biased Lasso Estimates ...   (8:27:21 PM)
#> Step 3: Joint significance test ...     (8:27:29 PM)
#> Done!     (8:27:29 PM)
#>   Index  alpha_hat   alpha_se   beta_hat   beta_se        IDE      rimp
#> 1     1 -0.4989049 0.06159061 -0.7522887 0.2407388  0.3753205  7.361963
#> 2     3  0.6351005 0.05489418  1.2125514 0.2837474  0.7700919 15.105458
#> 3     4  0.5432349 0.05966641 -0.8081071 0.2536428 -0.4389920  8.610887
#> 4     5  0.7310222 0.04849277 -1.0270749 0.2937140 -0.7508145 14.727329
#> 5     6 -0.7606745 0.04613191  0.9458118 0.2983506 -0.7194550 14.112207
#> 6     7 -0.7180654 0.04946084 -1.4210193 0.2790808  1.0203848 20.014986
#> 7     8 -0.7294041 0.04861566 -1.4025766 0.2884823  1.0230451 20.067169
#>           pmax
#> 1 1.778553e-03
#> 2 1.925594e-05
#> 3 1.442524e-03
#> 4 4.707808e-04
#> 5 1.523680e-03
#> 6 3.547092e-07
#> 7 1.162553e-06
```
HIMA identifies 7 out of the 8 active mediators


# Reference

- Wang, X., & Leng, C. (2016). High dimensional ordinary least squares
  projection for screening variables. *Journal of the Royal Statistical
  Society Series B: Statistical Methodology*, 78(3), 589–611.
- Battey, H. S., & Reid, N. (2023). On inference in high-dimensional
  regression. *Journal of the Royal Statistical Society Series B:
  Statistical Methodology*, 85(1), 149–175.
