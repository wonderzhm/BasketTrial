#' Posterior Inference for Simulated Basket Trial Data
#'
#' It generates posterior probabilities P(p_j > pnull) after all interim analysis
#' and calculates rates for early stopping, number of patients and estimated ORR.
#'
#' @param object returned from \link{generate.data}.
#' @param pnull B by 1 vector of null response rates, where B is the number of baskets.
#' @param stopbounds B by (stage-1) matrix: stopping boundaries for each basket at each interim.
#' @param clusterk only needed for parallel computing.
#' @param nperclust only needed for parallel computing.
#' @param beta.a0 a vector of length B for beta prior parameter a0 in each basket.
#' @param beta.b0 a vector of length B for beta prior parameter b0 in each basket.
#' @param seed random seed for reproducibility.
#' @param ModelFit the method function, e.g., \code{localPP}, \code{JSD}, and other user defined methods.
#' @param ... additional arguments passed to the method function defined by \code{ModelFit}.
#'
#' @return It returns a list including \code{data}, \code{N}, and \code{ORRs}, where \code{data} is an
#' array with \code{dim=c(nS, ntrial, B, stage)}.
#' @importFrom stats rbinom pbeta
#' @export
#'
#' @examples
#' N <- rbind(
#' c(10, 25),
#' c(10, 25),
#' c(10, 25)) # interim sample size and total sample size for each indication
#' scenarios <- rbind( c(0.15, 0.15, 0.15), c(0.3, 0.3, 0.3) )
#' res <- generate.data(N = N, ORRs = scenarios, ntrial = 20, seed = 2024)
#' post <- post.infer(res, pnull = rep(0.15,3), stopbounds = cbind(c(1,1,1)),
#' ModelFit = "localPP", method = "PEB")
#' apply(post$earlystop, c(1,3), mean) # early stopping for each basket in each scenario
#' apply(post$npts, c(1,3), mean) # average number of pts for each basket in each scenario
#' apply(post$est, c(1,3), mean) # average ORR estimate for each basket in each scenario
post.infer <- function(object, pnull, stopbounds = NULL, clusterk = NULL, nperclust = NULL,
                       beta.a0 = pnull, beta.b0 = 1-pnull, seed = 987897, ModelFit, ...) {
  # object: returned from generate.data
  # pnull: B-vector of null response rates
  # stopbounds: B by (stage-1) matrix: stopping boundaries for each indication at each interim
  # beta.a0, beta.b0: B-vector of prior values for the beta prior for each indication
  ORRs <- object$ORRs
  N <- object$N
  if(is.null(clusterk) | is.null(nperclust)){
    dat <- object$data
  }else{
    dat <- object$data[, (clusterk-1)*nperclust+(1:nperclust), , ,drop=FALSE]
  }
  stage <- dim(N)[2] # number of analyses (interim+final)
  Nmax <- N[,stage]
  B <- ncol(ORRs) # Number of Indications
  nS <- nrow(ORRs) # number of scenarios
  ntrial <- dim(dat)[2]
  Fit <- get(ModelFit)
  is.mc <- !(ModelFit %in% c("Independent", "localPP", "JSD", "localMEM"))
  is.MEM <- (ModelFit == "MEM")

  # Simulate Trials
  res.post <- array(NA, dim = c(nS, ntrial, B)) # posterior probabilities
  res.estop <- array(NA, dim = c(nS, ntrial, B)) # early stop status
  res.pts <- array(NA, dim = c(nS, ntrial, B)) # final enrolled number of patients
  res.est <- array(NA, dim = c(nS, ntrial, B)) # estimate of ORR

  set.seed(seed)

  for (s in 1:nS){
    for (trial in 1:ntrial) {
      # read observations for each Indication
      yobs <- array( dat[s, trial, , ], dim = dim(dat)[3:4] )
      if(stage==1){
        y <- yobs[, stage]
        n <- N[,stage]
        if(is.mc){
          if(is.MEM){
            fit0 <- Fit(nDat = n, yDat = y, p0 = pnull[1], ...)
            pp <- fit0$post_prob
            phat <- fit0$phat
          }else{
            fit0 <- Fit(nDat = n, yDat = y, ...)
            pp <- colMeans(fit0>matrix(pnull, nrow(fit0), ncol(fit0), byrow = TRUE))
            phat <- colMeans(fit0)
          }
        }else{
          fit0 <- Fit(nDat = n, yDat = y, be.a0 = beta.a0, be.b0 = beta.b0, ...)
          pp <- pbeta(pnull, fit0$a.post, fit0$b.post, lower.tail = FALSE)
          phat <- fit0$a.post/(fit0$a.post+fit0$b.post)
        }
        res.post[s, trial, ] <- pp
        res.estop[s, trial, ] <- rep(0, B)
        res.pts[s, trial, ] <- n
        res.est[s, trial, ] <- phat
      }else{
        Last_Stage = rep(1, B) # Keep Track of Last Stage
        y <- yobs[, 1]
        n <- N[, 1]
        stop.flag <- rep(FALSE, B)
        # interim analysis
        for (j in 1:(stage-1)){
          ind <- which(stop.flag==FALSE)
          if(length(ind)==0) break
          stop.flag[ind] <- (y[ind] <= stopbounds[ind,j])
          Last_Stage <- Last_Stage + !stop.flag
          y <- yobs[cbind(1:B, Last_Stage)]
          n <- N[cbind(1:B, Last_Stage)]
        }
        res.post[s, trial, ] <- pbeta(pnull, beta.a0+y, beta.b0+n-y, lower.tail = FALSE)
        res.estop[s, trial, ] <-  stop.flag+0
        res.est[s, trial, ] <- ((beta.a0+y)/(beta.b0+n-y))
        res.pts[s, trial, ] <- n

        # Final Analysis, where n contains the Final Sample Size
        ind.left <- which(stop.flag==FALSE)
        if(length(ind.left)>1){
          if(is.mc){
            if(is.MEM){
              fit0 <- Fit(nDat = n[ind.left], yDat = y[ind.left], p0 = pnull[1], ...)
              pp <- fit0$post_prob
              phat <- fit0$phat
            }else{
              fit0 <- Fit(nDat = n[ind.left], yDat = y[ind.left], ...)
              pp <- colMeans(fit0>matrix(pnull[ind.left], nrow(fit0), ncol(fit0), byrow = TRUE))
              phat <- colMeans(fit0)
            }
          }else{
            fit0 <- Fit(nDat = n[ind.left], yDat = y[ind.left],
                        be.a0 = beta.a0[ind.left], be.b0 = beta.b0[ind.left], ...)
            pp <- pbeta(pnull[ind.left], fit0$a.post, fit0$b.post, lower.tail = FALSE)
            phat <- (fit0$a.post/(fit0$a.post+fit0$b.post))
          }
          res.post[s, trial, ind.left] <- pp
          res.est[s, trial, ind.left] <- phat
        }
      }
    }
  }
  return (list(earlystop = res.estop, postprob = res.post, npts = res.pts, est = res.est,
               N = N, ORRs = ORRs, pnull = pnull, stopboundary = stopbounds))
}
