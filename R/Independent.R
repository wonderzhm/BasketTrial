#' Independent Beta Prior for Bayesian Basket Trial
#'
#' @param nDat a vector of length B for the sample size in each basket.
#' @param yDat a vector of length B for the number of responses in each basket.
#' @param be.a0 a vector of length B for beta prior parameter a0 in each basket.
#' @param be.b0 a vector of length B for beta prior parameter b0 in each basket.
#'
#' @return It returns a list including the posterior beta parameters.
#' @export
#'
#' @examples
#' Independent(nDat = c(25, 25, 25, 25, 25), yDat = c(2,9,11,13,20),
#' be.a0 = rep(0.5, 5), be.b0 = rep(0.5, 5))
Independent <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL){
  B <- length(nDat)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  return(list(a.post = be.a0 + yDat, b.post = be.b0 + nDat - yDat))
}

