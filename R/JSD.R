#' JSD (Fujikawa et al., 2020) Method for Bayesian Basket Trial
#'
#' @param nDat a vector of length B for the sample size in each basket.
#' @param yDat a vector of length B for the number of responses in each basket.
#' @param be.a0 a vector of length B for beta prior parameter a0 in each basket.
#' @param be.b0 a vector of length B for beta prior parameter b0 in each basket.
#' @param epsilon the global control parameter in the JSD model.
#' @param tau the threshold parameter in the JSD model.
#'
#' @return It returns a list including the posterior beta parameters and similarity matrix.
#' @importFrom stats integrate dbeta
#' @export
#'
#' @examples
#' JSD(nDat = c(25, 25, 25, 25, 25), yDat = c(2,9,11,13,20))
#' @references
#' Fujikawa, K., Teramukai, S., Yokota, I., & Daimon, T. (2020).
#' A Bayesian basket trial design that borrows information across strata based on
#' the similarity between the posterior distributions of the response probability.
#' Biometrical Journal, 62(2), 330-338.
## JSD (Fujikawa et al., 2020)
JSD <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL, epsilon = 2, tau = 0.3){
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
    # similarity matrix
    sm <- matrix(1, B, B)
    for(i in 1:(B-1)){
      fi <- function(x) dbeta(x, be.a0[i]+y[i], be.b0[i] + n[i] - y[i])
      for(j in (i+1):B){
        fj <- function(x) dbeta(x, be.a0[j]+y[j], be.b0[j] + n[j] - y[j])
        tmp1 <- integrate(function(x) fi(x)*log(2*fi(x)/(fi(x)+fj(x))),
                          lower=0.0001, upper=0.9999 )$value
        tmp2 <- integrate(function(x) fj(x)*log(2*fj(x)/(fi(x)+fj(x))),
                          lower=0.0001, upper=0.9999 )$value
        smij <- 1-(tmp1+tmp2)/2
        #sm[i,j] <- sm[j,i] <- smij*a*(abs(yo[i]/no[i]-yo[j]/no[j])<=delta)
        #pij <- prop.test(yo[c(i,j)], no[c(i,j)])$p.value
        #sm[i,j] <- sm[j,i] <- smij*a*(pij>=delta)
        sm[i,j] <- sm[j,i] <- (smij^epsilon)*((smij^epsilon)>tau)
      }
    }
    ab.post <- matrix(NA, B, 2)
    for(i in 1:B){
      ab.post[i,1] <- be.a0[i]+sum(sm[i,]*y)
      ab.post[i,2] <- be.b0[i]+sum(sm[i,]*(n-y))
    }
    return(list(a.post = ab.post[,1], b.post = ab.post[,2], sm = sm))
  }
}
