###########################################################################################
#### Function to compute the p-values to test H_{0}: \alpha_{j}=0 based on the OLS ####
###########################################################################################

# Inputs
# x : n-dimensional vectors of exposure
# chosen_M : n by q matrix of mediators chosen by some screening method
# COV.S: n by z matrix or data.frame of covariates


# Output
# ts: a vector of test statistics each of which asymptotically follow std. normal under the null
# pval: a vector of p-values to test \alpha_{j}=0


comp_alpha <- function(x, chosen_M, COV.S = NULL) {


  XC <- scale(x)

  p <- ncol(chosen_M)
  ts <- alpha_est <- pval <- rep(NA, p)


  if(!is.null(COV.S)) {
    COV.S <- scale(COV.S)
    XC <- cbind(XC, COV.S)
  }

  for(j in 1:p) {
    res <- stats::coef(summary(stats::lm(chosen_M[, j] ~ XC)))

    ts[j] <- res[2, 3] ; pval[j] <- res[2, 4]; alpha_est[j] <- res[2, 1]
  }

  return(list(ts = ts, pval = pval, alpha_est = alpha_est))
}


###########################################################################################
#### 17. function to conduct the joint significance test as the existing approaches do ####
###########################################################################################

# Inputs
# chosen_ind: a vector of indices for chosen mediators
# pval_alp: a vector of RAW p-values to test H_{0}: \alpha_{j} = 0 where \alpha_{j} is the effect of the exposure on the j-th mediator
# pval_beta: a vector of RAW p-values to test H_{0}: \beta_{j} = 0 where \beta_{j} is the effect of the j-th mediator on the response
# method: a method used to adjust each p-value that will be plugged as an argument of the p.adjust function, see ?p.adjust for more details. HDMT is the default
# alpha: a significance level
# output
# which_sig: an index set for active mediators that turn out to be significant

js_test <- function(chosen_ind, pval_alp, pval_beta, method = "HDMT", alpha = 0.05) {

  # Check for valid method
  method <- tolower(method)
  if (!method %in% c("hdmt", "bonferroni")) {
    stop("Invalid method. Choose either 'HDMT' or 'bonferroni'.")
  }

  if (method == "hdmt") {
    # Add small noise to avoid ties
    PA <- cbind(pval_alp, pval_beta)
    N0 <- nrow(PA) * ncol(PA)
    input_pvalues <- PA + matrix(stats::runif(N0, 0, 1e-10), nrow(PA), 2)

    # Estimate null proportions
    nullprop <- HDMT::null_estimation(input_pvalues)

    # Estimate FDR
    fdrcut <- HDMT::fdr_est(
      alpha00 = nullprop$alpha00,
      alpha01 = nullprop$alpha01,
      alpha10 = nullprop$alpha10,
      alpha1 = nullprop$alpha1,
      alpha2 = nullprop$alpha2,
      input_pvalues,
      exact = 0
    )

    which_sig <- chosen_ind[fdrcut <= alpha]
  }

  if (method == "bonferroni") {
    adjp_alp <- stats::p.adjust(pval_alp, method = "bonferroni")
    adjp_beta <- stats::p.adjust(pval_beta, method = "bonferroni")
    max_pval <- pmax(adjp_alp, adjp_beta)
    which_sig <- chosen_ind[max_pval <= alpha]
  }

  return(which_sig)
}
