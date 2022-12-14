% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grpreg_functions.R
\name{testMethods}
\alias{testMethods}
\title{testMethods}
\usage{
testMethods(
  a = 0.9,
  b = 0.5,
  n = 500,
  p = 500,
  ng = 10,
  truegroup = sort(c(rep(1:(p/ng), ng))),
  withinGroupNoise = 0,
  chooseBetas = FALSE,
  methods = "lasso"
)
}
\arguments{
\item{a}{Within group correlation (value from -1 to 1) for data generation}

\item{b}{Between group correlation (value from -1 to 1) for data generation}

\item{n}{Sample size for data generation}

\item{p}{Number of input variables for data generation}

\item{ng}{Number of groups for data generation}

\item{truegroup}{True grouping structure of input variables (p-dim vector with numerical true group labels)}

\item{withinGroupNoise}{Proportion of coefficients that switch from category 0 to category 1}

\item{chooseBetas}{True or false variable that indicates to take only last ng coefficients from each group to run methods on}

\item{methods}{Specified method ("lasso","ridge","group lasso","threshold group lasso","sparse group lasso","elastic net","group bridge")}
}
\description{
Generates design matrix with user-specified grouping structure and runs user-specified group regularization methods
}
\examples{
testMethods(a=0.9,b = 0.3,n=600,p=500,ng=10,methods=c("lasso","ridge","group lasso","threshold group lasso","elastic net","sparse group lasso","group bridge"))

}
