% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_active_med.R
\name{get_active_med}
\alias{get_active_med}
\title{Identify Active Mediators}
\usage{
get_active_med(
  y,
  x,
  M,
  COV.S = NULL,
  pval.adjust = "HDMT",
  d = NULL,
  r = 1,
  k = 1,
  alpha = 0.05
)
}
\arguments{
\item{y}{vector of outcomes}

\item{x}{vector of exposures}

\item{M}{a \code{data.frame} or \code{matrix} of mediators}

\item{COV.S}{a \code{data.frame} or \code{matrix} of covariates}

\item{pval.adjust}{specifies which method to use for controlling FWER/FDR in the joint significance testing. Either \code{'HDMT'} (default) or \code{'Bonferroni'}}

\item{d}{the number of screened mediators. Default value is \eqn{d =  n/\log(n)}.}

\item{r}{a penalty parameter for the Ridge-HOLP. Default value is \code{1}}

\item{k}{a scalar for computing projection directions for AO. Default value is \code{1}.}

\item{alpha}{Target FDR level. Default value is \code{0.05}.}
}
\value{
A data frame with columns
\describe{
\item{Index}{Column index of active mediators in \code{M}.}
\item{alpha_hat}{Estimated \eqn{\alpha_j}.}
\item{beta_hat}{Estimated \eqn{\beta_j}.}
\item{P_value}{Raw joint p‑value \eqn{\max(p_{\alpha j}, p_{\beta j})}.}
}
The data frame has zero rows if no active mediator is detected.
}
\description{
Identify Active Mediators
}
\examples{
data(ExampleData) # Load ExampleData from the package
y <- ExampleData$y # 'y' is a vector of outcomes in ExampleData
x <- ExampleData$x #  'x' is a vector of exposures in ExampleData
M <- ExampleData$M #  'M' is a matrix of mediators in ExampleData

# Get active mediators
results <- get_active_med(y, x, M)

# Print the results of active mediators
print(results)


}
