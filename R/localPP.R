#' Local Power Prior for Bayesian Basket Trial
#'
#' @param nDat a vector of length B for the sample size in each basket.
#' @param yDat a vector of length B for the number of responses in each basket.
#' @param be.a0 a vector of length B for beta prior parameter a0 in each basket.
#' @param be.b0 a vector of length B for beta prior parameter b0 in each basket.
#' @param a the global control parameter in the local PP 3-component framework.
#' @param delta the threshold parameter in the local PP 3-component framework.
#' @param method either \code{PEB} for the pairwise empirical Bayes or \code{GEB}
#' for the global empirical Bayes.
#' @param symmetry logical variable to indicate whether the similarity matrix
#' will be set to be symmetric; default is \code{FALSE}.
#'
#' @return It returns a list including the posterior beta parameters and similarity matrix.
#' @importFrom stats optimize optim
#' @export
#'
#' @examples
#' localPP(nDat = c(25, 25, 25, 25, 25), yDat = c(2,9,11,13,20),
#' be.a0 = rep(0.5, 5), be.b0 = rep(0.5, 5), a = 4, delta = 1, method = "PEB")
#' @references
#' Zhou, H., Shen, R., Wu, S., & He, P. (2023). A Bayesian Basket Trial Design Using Local Power Prior.
#' arXiv preprint arXiv:2312.15352.
localPP <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL, a = 1,
                    delta = 0.4, method = "PEB", symmetry = FALSE){
  n <- nDat
  y <- yDat
  B <- length(n)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  if(B==1){
    return(list(a.post = be.a0 + y, b.post = be.b0 + n - y, sm = matrix(1)))
  }else{
    # maximum borrowing for each basket
    amax <- rep(NA, B)
    for(i in 1:B) amax[i] <- min(a*n[i]/sum(n[-i]), 1)
    # similarity matrix
    sm <- matrix(1, B, B)
    if(method=="PEB"){
      fn <- function(theta, yi, ni, y0, n0, a0, b0) {
        (lbeta(a0+yi+theta*y0, b0+(ni-yi)+theta*(n0-y0))
         -lbeta(a0+theta*y0, b0+theta*(n0-y0)))
      }
      for(i in 1:B){
        for(j in (1:B)[-i]){
          res <- optimize(f = fn, interval = c(0, 1),
                          yi = y[i], ni = n[i], y0 = y[j], n0 = n[j],
                          a0 = be.a0[i], b0 = be.b0[i], maximum = TRUE)
          sm[i,j] <- res$maximum*(amax[i])*(abs(y[i]/n[i]-y[j]/n[j])<delta)
        }
      }
    }else{
      fn <- function(theta, yi, ni, y0, n0, a0, b0) {
        -2*(lbeta(a0+yi+sum(theta*y0), b0+(ni-yi)+sum(theta*(n0-y0)))
            -lbeta( a0+sum(theta*y0), b0+sum(theta*(n0-y0)) ))
      }
      for(i in 1:B){
        res <- optim(par=rep(0.5,B-1), fn = fn, method = "L-BFGS-B",
                     lower = rep(0,B-1), upper = rep(1,B-1),
                     yi = y[i], ni = n[i], y0 = y[-i], n0 = n[-i],
                     a0 = be.a0[i], b0 = be.b0[i])
        sm[i,-i] <- res$par*(amax[i])*(abs(y[i]/n[i]-y[-i]/n[-i])<delta)
      }
    }
    if (symmetry) sm <- (t(sm) + sm)/2
    ab.post <- matrix(NA, B, 2)
    for(i in 1:B){
      ab.post[i,1] <- be.a0[i]+sum(sm[i,]*y)
      ab.post[i,2] <- be.b0[i]+sum(sm[i,]*(n-y))
    }
    return(list(a.post = ab.post[,1], b.post = ab.post[,2],
                sm = sm))
  }
}

