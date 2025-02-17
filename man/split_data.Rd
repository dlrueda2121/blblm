% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{split_data}
\alias{split_data}
\title{Splits data}
\usage{
split_data(data, m)
}
\arguments{
\item{data}{dataframe}

\item{m}{integer}
}
\value{
list
}
\description{
Splits data into m parts of approximated equal sizes
}
\examples{
split_data(data = mtcars, m = 3)
}
