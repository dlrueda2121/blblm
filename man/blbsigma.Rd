% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{blbsigma}
\alias{blbsigma}
\title{Sigma}
\usage{
blbsigma(fit)
}
\arguments{
\item{fit}{blblm object}
}
\value{
double
}
\description{
Computes sigma from fit
}
\examples{
blbsigma(blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100))
}
