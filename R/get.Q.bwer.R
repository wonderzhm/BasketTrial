#' Get Efficacy Cutoff based on Basket-Wise Error Rate (BWER) Control
#'
#' @param object returned by \link{post.infer}.
#' @param alpha basket-wise type I error control.
#' @param digits number of digits in the cutoffs.
#' @param Qclust \code{NULL} means all cutoffs are different; If there are B=5 baskets and
#' \code{Qclust=(1,1,2,2,2)}, it means cutoffs for the first two baskets will be the same
#' and another cutoff will be tuned separately for baskets 3-5.
#'
#' @return It returns the efficacy cutoffs.
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' N <- rbind(
#' c(10, 25),
#' c(10, 25),
#' c(10, 25),
#' c(10, 25),
#' c(10, 25)) # interim sample size and total sample size for each indication
#' scenarios <- rbind( c(0.15, 0.15, 0.15, 0.15, 0.15), c(0.3, 0.3, 0.3, 0.3, 0.3) )
#' res <- generate.data(N = N, ORRs = scenarios, ntrial = 1000, seed = 343809)
#' post <- post.infer(res, pnull = rep(0.15,5), stopbounds = cbind(c(1,1,1,1,1)),
#' ModelFit = "localPP", method = "PEB", a = 2, delta = 0.3)
#' (Q <- get.Q.bwer(post, alpha = 0.1, digits = 3, Qclust = rep(1, 5)))
#' Qmat <- array(NA, dim = dim(post$postprob))
#' for(i in 1:5) Qmat[,,i] <- Q[i]
#' apply(post$postprob>Qmat, c(1,3), mean)
get.Q.bwer <- function(object, alpha = 0.1, digits = 3, Qclust = NULL){
  #object: returned by post.infer()
  #Qclust: =NULL means all Qs are different;
  #        If there are B=5 baskets and Qclust=(1,1,2,2,2), it means Q for
  #        the first two baskets will be the same and another Q will be used for
  #        baskets 3-5.
  res.post <- object$postprob ### dim(nS, ntrial, B)
  pnull <- object$pnull
  Nmax <- object$N[,ncol(object$N)]
  nS <- dim(res.post)[1]
  ntrial <- dim(res.post)[2]
  B <- dim(res.post)[3]
  ORRs <- object$ORRs
  Q.final <- rep(NA, B)
  if(is.null(Qclust)) Qclust <- 1:B
  Q.unique <- unique(Qclust)
  nQ <- length(Q.unique)
  for(i in 1:nQ){
    i.ind <- which(Qclust==Q.unique[i])
    xxx <- as.vector(res.post[1, , i.ind])
    Q <- quantile(xxx, probs = 1-alpha, names = FALSE)
    Q0 <- floor(Q*10^digits)/(10^digits)
    Q1 <- ceiling(Q*10^digits)/(10^digits)
    error0 <- mean(xxx>Q0)
    error1 <- mean(xxx>Q1)
    Qfinal <- ifelse(abs(error0-alpha)<abs(error1-alpha), Q0, Q1)
    error.final <- mean(xxx>Qfinal)
    # for independent models, the error may be much higher than alpha
    Q.final[i.ind] <- ifelse(error.final<(alpha*1.05), Qfinal, Qfinal+1/(10^digits))
  }
  Q.final
}
