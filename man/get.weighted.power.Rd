% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.weighted.power.R
\name{get.weighted.power}
\alias{get.weighted.power}
\title{Get Weighted Type I Error and Power.}
\usage{
get.weighted.power(object, Q, s0 = 100, s1 = 0)
}
\arguments{
\item{object}{returned by \link{post.infer}.}

\item{Q}{a vector of length B for the efficacy cutoff in each basket.}

\item{s0}{Setting s0=100 the weighted power reduces to type I error under global null.
Please use this default.}

\item{s1}{Setting s1=0 gives equal weight for calculating weighted power across scenarios.
Please use this default.}
}
\value{
It returns a list with \code{error.tw} for average basket-wise type I error rate (BWER) under global null,
\code{bwer} for BWERs for all null baskets,
\code{power.cdr} for average true positive rate (TPR) across scenarios except global null,
\code{power.ccr} for average correct classification rate (CCR) across scenarios except global null.
}
\description{
Get weighted type I error (WE) and power(WP) cross all scenarios.
including family wise (fwer) or trial wise (twer) or false discovery rate (fdr).
}
\examples{
N <- rbind(
c(10, 25),
c(10, 25),
c(10, 25),
c(10, 25),
c(10, 25)) # interim sample size and total sample size for each indication
scenarios <- rbind( c(0.15, 0.15, 0.15, 0.15, 0.15), c(0.3, 0.3, 0.3, 0.3, 0.3) )
res <- generate.data(N = N, ORRs = scenarios, ntrial = 1000, seed = 343809)
post <- post.infer(res, pnull = rep(0.15,5), stopbounds = cbind(c(1,1,1,1,1)),
ModelFit = "localPP", method = "PEB", a = 2, delta = 0.3)
(Q <- get.Q.bwer(post, alpha = 0.1, digits = 3, Qclust = rep(1, 5)))
(powers <- get.weighted.power(object = post, Q = Q))
}
