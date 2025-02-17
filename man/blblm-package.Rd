% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\docType{package}
\name{blblm-package}
\alias{blblm-package}
\alias{_PACKAGE}
\title{blblm: Linear Regression with Bag of Little Bootstraps}
\description{
Implements linear regression through Bag of Little Bootstraps procedure. Combines features from bootstrap and subsampling with the option to use multiple cores for faster processing.
    Samples n from b with replacement. Samples s times into sizes of b. Then for each subsample it resamples each r times until sample size is n. Computes the bootstrap statistic for each bootstrap sample. Computes the statistic from the bootstrap statistics, and finally takes the average of the statistics.
}
\details{
Linear Regression with Little Bag of Bootstraps
}
