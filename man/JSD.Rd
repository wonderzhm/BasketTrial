% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JSD.R
\name{JSD}
\alias{JSD}
\title{JSD (Fujikawa et al., 2020) Method for Bayesian Basket Trial}
\usage{
JSD(nDat, yDat, be.a0 = NULL, be.b0 = NULL, epsilon = 2, tau = 0.3)
}
\arguments{
\item{nDat}{a vector of length B for the sample size in each basket.}

\item{yDat}{a vector of length B for the number of responses in each basket.}

\item{be.a0}{a vector of length B for beta prior parameter a0 in each basket.}

\item{be.b0}{a vector of length B for beta prior parameter b0 in each basket.}

\item{epsilon}{the global control parameter in the JSD model.}

\item{tau}{the threshold parameter in the JSD model.}
}
\value{
It returns a list including the posterior beta parameters and similarity matrix.
}
\description{
JSD (Fujikawa et al., 2020) Method for Bayesian Basket Trial
}
\examples{
JSD(nDat = c(25, 25, 25, 25, 25), yDat = c(2,9,11,13,20))
}
\references{
Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
A Bayesian basket trial design that borrows information across strata based on
the similarity between the posterior distributions of the response probability.
Biometrical Journal, 62(2), 330-338.
}
