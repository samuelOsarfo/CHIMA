#' Identify Active Mediators
#'
#' @param y  vector of outcomes
#' @param x vector of exposures
#' @param M a \code{data.frame} or \code{matrix} of mediators
#' @param COV.S a \code{data.frame} or \code{matrix} of covariates
#' @param pval.adjust specifies which method to use for controlling FWER/FDR in the joint significance testing. Either \code{'HDMT'} (default) or \code{'Bonferroni'}
#' @param d the number of screened mediators. Default value is \eqn{d =  n/\log(n)}.
#' @param r  a penalty parameter for the Ridge-HOLP. Default value is `1`
#' @param k  a scalar for computing projection directions for AO. Default value is `1`.
#' @param alpha  Target FDR level. Default value is `0.05`.
#'
#' @return A data frame with columns
#'   \describe{
#'     \item{Index}{Column index of active mediators in \code{M}.}
#'     \item{alpha_hat}{Estimated \eqn{\alpha_j}.}
#'     \item{beta_hat}{Estimated \eqn{\beta_j}.}
#'     \item{P_value}{Raw joint p‑value \eqn{\max(p_{\alpha j}, p_{\beta j})}.}
#'   }
#'   The data frame has zero rows if no active mediator is detected.#'
#' @export
#'
#' @examples
#' data(ExampleData) # Load ExampleData from the package
#' y <- ExampleData$y # 'y' is a vector of outcomes in ExampleData
#' x <- ExampleData$x #  'x' is a vector of exposures in ExampleData
#' M <- ExampleData$M #  'M' is a matrix of mediators in ExampleData
#'
#' # Get active mediators
#' results <- get_active_med(y, x, M)
#'
#' # Print the results of active mediators
#' print(results)
#'
#'
get_active_med <- function(y, x, M, COV.S=NULL, pval.adjust='HDMT', d=NULL, r=1,  k=1, alpha=0.05){

    #screen mediators
    msg("Step 1: Ridge–HOLP Screening   -----  ", format(Sys.time(), "%H:%M:%S  %p"))
    chosen_ind <- medsc_holp(y, x, M, COV.S, d, r)


    #get bhats, p-values, and test statistics for \beta_j's in outcome model
    msg("Step 2: Approximate Orthogonalization Estimates   -----  ", format(Sys.time(), "%H:%M:%S %p"))
    ao_obj <- app_orth(y, x, M[, chosen_ind], COV.S, k)


    #get p-values and test statistics for \alpha_j's in mediator model
    alp_all <- comp_alpha(x, M[,chosen_ind], COV.S)

    #get index for active mediators in chosen mediators
    msg("Step 3: Joint Significance Testing   -----  ", format(Sys.time(), "%H:%M:%S %p"))
    active_index <- tryCatch(
      js_test(chosen_ind,
              alp_all$pval,
              ao_obj$pval,
              method = pval.adjust,
              alpha),
      error = function(e) {
        msg("js_test() failed: ", conditionMessage(e),
             ". Returning empty result.")
        integer(0)
      }
    )

    if (length(active_index) == 0L || all(is.na(active_index))) {
      msg("No active mediators detected. Returning empty result.")
      return(empty_active_df())
    }

     #get raw pvalues for active mediators
     Pvalues <-cbind(ao_obj$pval[active_index], alp_all$pval[active_index])
     Pvalues<- apply(Pvalues, 1, max)

    #results
    msg("Complete!!   ", format(Sys.time(), "%H:%M:%S %p"))


    data.frame(Index=active_index,
                          alpha_hat=alp_all$alpha_est[active_index],
                          beta_hat=ao_obj$bhat[active_index],
                          P_value=Pvalues,
                          row.names=colnames(M)[active_index]
                          )

}
