#' Generate Data for A Basket Trial Design
#'
#' @param N a matrix with \code{dim=c(B, stage)}, where \code{B} is the number of baskets and
#' \code{stage} is the # of analyses (interim+final)
#' @param ORRs a matrix with \code{dim=c(nS, B)}, where \code{nS} is the number of trial scenarios for objective
#' response rates.
#' @param ntrial the total number of trials simulated.
#' @param seed random seed for reproducibility.
#'
#' @return It returns a list including \code{data}, \code{N}, and \code{ORRs}, where \code{data} is an
#' array with \code{dim=c(nS, ntrial, B, stage)}.
#' @importFrom stats rbinom
#' @export
#'
#' @examples
#' N <- rbind(
#' c(10, 25),
#' c(10, 25),
#' c(10, 25)) # interim sample size and total sample size for each indication
#' scenarios <- rbind( c(0.15, 0.15, 0.15), c(0.3, 0.3, 0.3) )
#' res <- generate.data(N = N, ORRs = scenarios, ntrial = 20, seed = 2024)
generate.data <- function(N, ORRs, ntrial = 10000, seed = 987897){
  # N: matrix with dim=(B, stage), where stage is the # of analyses (interim+final)
  # ORRs: a matrix with dim = (nS, B)
  set.seed(seed)
  if(is.vector(ORRs)){
    ORRs <- matrix(ORRs, nrow = 1)
  }
  if(is.vector(N)){
    N <- matrix(N)
  }
  stage <- dim(N)[2] # number of analyses (interim+final)
  B <- ncol(ORRs) # Number of Indications
  nS <- nrow(ORRs) # number of scenarios
  res <- array(NA, dim = c(nS, ntrial, B, stage))
  for (s in 1:nS){
    for (trial in 1:ntrial) {
      # Generate Observations for each Indication
      res[s, trial, , 1] <- rbinom(B, N[, 1], ORRs[s,])
      if(stage>1){
        for(j in 2:stage){
          res[s, trial, , j] <- res[s, trial, , j-1] + rbinom(B, N[, j]-N[, j-1], ORRs[s,])
        }
      }
    }
  }
  return(list(data = res, N = N, ORRs = ORRs))
}
