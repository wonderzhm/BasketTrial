library(shiny)
library(DT)
library(rhandsontable)
library(xtable)
library(partitions)
library(parallel)
library(doParallel)
library(foreach)
Independent <- function(nDat, yDat, be.a0 = NULL, be.b0 = NULL){
  B <- length(nDat)
  if(is.null(be.a0) | is.null(be.b0)) {
    be.a0 = rep(.15, B)
    be.b0 = rep(.85, B)
  }
  return(list(a.post = be.a0 + yDat, b.post = be.b0 + nDat - yDat))
}

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

get.weighted.power <- function(object, Q, s0 = 100, s1 = 0){
  # Setting s0=100 the weighted power reduces to type I error under global null
  # Setting s1=0 gives equal weight for calculating weighted power across scenarios.
  #object: returned by post.infer()
  res.post <- object$postprob ### dim(nS, ntrial, B)
  pnull <- object$pnull
  nS <- dim(res.post)[1]
  ntrial <- dim(res.post)[2]
  B <- dim(res.post)[3]
  ORRs <- object$ORRs
  pnull.mat <- matrix(pnull, nS, B, byrow = TRUE)
  ## weights
  H0status <- (ORRs<=pnull.mat)
  nH0 <- rowSums(H0status)
  ind0 <- which(nH0!=0)
  wei0 <- rep(0, nS)
  wei0[ind0] <-  nH0[ind0]^s0/sum(nH0[ind0]^s0)
  H1status <- (ORRs>pnull.mat)
  nH1 <- rowSums(H1status)
  ind1 <- which(nH1!=0)
  wei1 <- rep(0, nS)
  wei1[ind1] <-  nH1[ind1]^s1/sum(nH1[ind1]^s1)
  
  # Errors and powers
  Qarray <- array(NA, dim = c(nS, ntrial, B))
  for(i in 1:B) Qarray[,,i] <- Q[i]
  res.rej <- (res.post>Qarray)
  fdrs <- rep(NA, nS)
  twerrors <- rep(NA, nS)
  fwerrors <- rep(NA, nS)
  cdrs <- rep(NA, nS)
  adrs <- rep(NA, nS)
  ccrs <- rep(NA, nS)
  bwer <- matrix(NA, nS, B)
  
  for(i in 1:nS){
    rate.rej <- matrix(res.rej[i,,], ntrial, B)
    x0 <- rate.rej[,ORRs[i,]<=pnull,drop=FALSE]
    x1 <- rate.rej[,ORRs[i,]>pnull,drop=FALSE]
    n.rej <- rowSums(rate.rej)
    n.fd <- ifelse(n.rej==0, 0, rowSums(x0)/n.rej)
    fdrs[i] <- ifelse(length(x0)==0, NA, mean(n.fd))
    fwerrors[i] <- ifelse(length(x0)==0, NA, mean(apply(x0, 1, any)))
    twerrors[i] <- ifelse(length(x0)==0, NA, mean(colMeans(x0)))
    #fwps[i] <- ifelse(length(x1)==0, NA, mean(apply(x1, 1, any)))
    #twps[i] <- ifelse(length(x1)==0, NA, mean(colMeans(x1)))
    cdrs[i] <- ifelse(length(x1)==0, NA, mean(rowMeans(x1)))
    adrs[i] <- ifelse(length(x1)==0, NA, mean((rowSums(x1)-rowSums(x0))/ncol(x1)))
    ccrs[i] <- ifelse(length(x1)==0, NA, mean((rowSums(x1)+(ncol(x0)-rowSums(x0)))/B))
    bwer[i,] <- colMeans(rate.rej)
    bwer[i,ORRs[i,]>pnull] <- NA
  }
  w.fdr <- weighted.mean(fdrs, w = wei0)
  w.twe <- weighted.mean(twerrors, w = wei0)
  w.fwe <- weighted.mean(fwerrors, w = wei0)
  w.cdr <- weighted.mean(cdrs, w = wei1)
  w.adr <- weighted.mean(adrs, w = wei1)
  w.ccr <- weighted.mean(ccrs, w = wei1)
  
  return(list(Q=Q, error.fdr = w.fdr, error.tw=w.twe, error.fw=w.fwe,
              power.cdr=w.cdr, power.adr = w.adr, power.ccr = w.ccr,
              ind.error.fdr = fdrs, ind.error.tw = twerrors, ind.error.fw = fwerrors,
              ind.power.cdr=cdrs, ind.power.adr=adrs, ind.power.ccr=ccrs, bwer = bwer,
              w0 = wei0, w1 = wei1))
}

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
parseNumericInput <- function(x) {
  nums <- as.numeric(unlist(strsplit(gsub(" ", "", x), ",")))
  nums[!is.na(nums)]
}
post.infer.parallel <- function(nclust, nperclust, data.object, pnull, stopbounds = NULL, 
                                beta.a0 = pnull, beta.b0 = 1-pnull, seed = 38098, ModelFit, ...) {
  if(nclust > detectCores()) stop("nclust needs to be smaller than detectCores()")
  
  ## Create Clusters
  cl <- makeCluster(nclust, type = "PSOCK")
  registerDoParallel(cl)
  
  ## Export Functions to the Cluster
  clusterCall(cl, function() {source(file.path("..", "functions", "functions.R"))})
  clusterExport(cl, list("nperclust", "data.object", "pnull", "stopbounds", "beta.a0", "beta.b0", "seed"))
  
  ## Parallel Computing
  results <- parLapply(cl, (1:nclust), fun = post.infer, object = data.object, 
                       pnull = pnull, stopbounds = stopbounds, nperclust = nperclust, 
                       beta.a0 = beta.a0, beta.b0 = beta.b0, seed = seed,
                       ModelFit = ModelFit, ...)
  
  ## Close Cluster
  stopCluster(cl)
  
  ## Combine the Results
  resultsAll <- results[[1]]
  for(i in 2:nclust){
    resultsAll$earlystop <- abind::abind(resultsAll$earlystop, results[[i]]$earlystop, along=2)
    resultsAll$postprob <- abind::abind(resultsAll$postprob, results[[i]]$postprob, along=2)
    resultsAll$npts <- abind::abind(resultsAll$npts, results[[i]]$npts, along=2)
    resultsAll$est <- abind::abind(resultsAll$est, results[[i]]$est, along=2)
  }
  return(resultsAll)
}
server <- function(input, output, session) {
  
  # ------------------------
  # Trial Design
  # ------------------------
  
  output$basket_sample_sizes <- renderUI({
    lapply(1:input$num_baskets, function(i) {
      numericInput(paste0("n_basket_", i),
                   paste("Sample size for Basket", i),
                   value = 20, min = 1)
    })
  })
  
  output$interim_inputs <- renderUI({
    if (input$num_interim > 0) {
      lapply(1:input$num_baskets, function(b) {
        lapply(1:input$num_interim, function(k) {
          fluidRow(
            column(6, numericInput(paste0("ia_", b, "_", k),
                                   paste("Basket", b, "pts at IA", k),
                                   value = 5, min = 1)),
            column(6, numericInput(paste0("stop_", b, "_", k),
                                   paste("Stop if responses <="),
                                   value = -1))
          )
        })
      })
    }
  })
  
  observeEvent(input$p0, {
    updateNumericInput(session, "b1", value = input$p0)
    updateNumericInput(session, "b2", value = 1 - input$p0)
  })
  
  observe({
    updateSelectInput(session, "design_method", selected = "local-PP-PEB")
  })
  
  cutoff_values <- reactiveVal(NULL)
  
  observeEvent(input$calculate_cutoffs, {
    withProgress(message = 'Calculating efficacy cutoffs...', value = 0, {
      
      nB <- input$num_baskets
      K <- input$num_interim
      
      incProgress(0.1, detail = "Setting up trial design...")
      
      if (K == 0) {
        nDat <- matrix(sapply(1:nB, function(i) input[[paste0("n_basket_", i)]]), 
                       nrow = nB, ncol = 1)
      } else {
        nDat <- matrix(NA, nrow = nB, ncol = K + 1)
        for (b in 1:nB) {
          for (k in 1:K) {
            nDat[b, k] <- input[[paste0("ia_", b, "_", k)]]
          }
          nDat[b, K + 1] <- input[[paste0("n_basket_", b)]]
        }
      }
      
      ORRs <- matrix(rep(input$p0, nB), nrow = 1)
      
      incProgress(0.2, detail = "Generating simulation data...")
      
      ntrial_cutoffs <- ifelse(is.null(input$ntrial_design), 5000, input$ntrial_design)
      data.object <- generate.data(N = nDat, ORRs = ORRs, ntrial = ntrial_cutoffs)
      pnull <- rep(input$p0, nB)
      
      stopbounds <- NULL
      if (K > 0) {
        stopbounds <- matrix(NA, nrow = nB, ncol = K)
        for (b in 1:nB) {
          for (k in 1:K) {
            stopbounds[b, k] <- input[[paste0("stop_", b, "_", k)]]
          }
        }
      }
      
      method <- input$design_method
      
      incProgress(0.4, detail = paste("Running", method, "analysis..."))
      
      postinf <- NULL
      if(method == "Independent"){
        postinf <- tryCatch({
          post.infer(object = data.object, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = method)
        }, error = function(e) {
          showNotification(paste("Error in", method, ":", e$message), type = "error")
          return(NULL)
        })
      } else if (method == "local-PP-GEB"){
        if (is.null(input$a) || is.null(input$delta)) {
          showNotification("Please set parameters 'a' and 'delta' for local-PP-GEB method", type = "error")
          return(NULL)
        }
        postinf <- tryCatch({
          post.infer(object = data.object, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = "localPP", method = "GEB", a = input$a, delta = input$delta)
        }, error = function(e) {
          showNotification(paste("Error in", method, ":", e$message), type = "error")
          return(NULL)
        })
      } else if (method == "local-PP-PEB"){
        if (is.null(input$a) || is.null(input$delta)) {
          showNotification("Please set parameters 'a' and 'delta' for local-PP-PEB method", type = "error")
          return(NULL)
        }
        postinf <- tryCatch({
          post.infer(object = data.object, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = "localPP", method = "PEB", a = input$a, delta = input$delta)
        }, error = function(e) {
          showNotification(paste("Error in", method, ":", e$message), type = "error")
          return(NULL)
        })
      } else if (method == "JSD"){
        if (is.null(input$epsilon) || is.null(input$tau)) {
          showNotification("Please set parameters 'epsilon' and 'tau' for JSD method", type = "error")
          return(NULL)
        }
        postinf <- tryCatch({
          post.infer(object = data.object, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = "JSD", epsilon = input$epsilon, tau = input$tau)
        }, error = function(e) {
          showNotification(paste("Error in", method, ":", e$message), type = "error")
          return(NULL)
        })
      }
      
      if(is.null(postinf)) return(NULL)
      
      incProgress(0.8, detail = "Computing efficacy cutoffs...")
      
      if (K == 0) {
        sample_sizes <- nDat[, 1]
      } else {
        sample_sizes <- nDat[, ncol(nDat)]
      }
      
      unique_sizes <- unique(sample_sizes)
      Qclust <- match(sample_sizes, unique_sizes)
      
      Q <- get.Q.bwer(postinf, alpha = input$alpha, Qclust = Qclust)
      cutoffs <- data.frame(
        Basket = paste("Basket", 1:nB),
        Efficacy_Cutoff = Q
      )
      cutoff_values(cutoffs)
      
      incProgress(1.0, detail = "Complete!")
      
      showNotification("Efficacy cutoffs calculated successfully!", type = "message", duration = 3)
    })
  })
  
  output$efficacy_cutoff_table <- renderTable({
    cutoff_data <- cutoff_values()
    if(!is.null(cutoff_data)) {
      cutoff_data$Efficacy_Cutoff <- round(cutoff_data$Efficacy_Cutoff, 3)
    }
    cutoff_data
  }, digits = 3)
  
  # ------------------------
  # Simulation
  # ------------------------
  
  scenarios_data <- reactiveVal()
  sim_results <- reactiveVal(NULL)
  
  observeEvent(input$num_baskets, {
    req(input$num_baskets)
    nB <- input$num_baskets
    if(nB > 0 && nB <= 20) {
      current_data <- scenarios_data()
      
      if(is.null(current_data) || ncol(current_data) != nB) {
        # Create data frame without Scenario column
        initial_data <- data.frame(matrix(NA, nrow = 3, ncol = nB))
        colnames(initial_data) <- paste0("Basket_", 1:nB)
        rownames(initial_data) <- paste0("S", 1:3)
        
        for(i in 1:nB) {
          if(i <= 3) {
            initial_data[, i] <- c(0.15, 0.15, 0.30)
          } else if(i == 4 || i == 5) {
            initial_data[, i] <- c(0.15, 0.30, 0.30)
          } else {
            initial_data[, i] <- c(0.15, 0.15, 0.30)
          }
        }
        
        scenarios_data(initial_data)
        showNotification(paste("Scenario table updated for", nB, "baskets"), type = "message", duration = 2)
      }
    }
  })
  
  output$simulation_scenarios_ui <- renderUI({
    req(input$num_baskets)
    
    nB <- input$num_baskets
    if(nB <= 0 || nB > 20) return(div("Please configure Trial Design first"))
    
    trial_design_configured <- TRUE
    for(i in 1:nB) {
      if(is.null(input[[paste0("n_basket_", i)]])) {
        trial_design_configured <- FALSE
        break
      }
    }
    
    tagList(
      if(!trial_design_configured) {
        div(class = "alert alert-warning",
            style = "margin: 10px 0;",
            "Trial Design not configured. Using K=", 
            ifelse(is.null(input$num_interim) || input$num_interim == 0, "0", input$num_interim),
            ", ", nB, " baskets, and a sample size of 20 per basket."
        )
      },
      
      h4("Scenario Configuration"),
      div(style = "margin: 10px 0;",
          actionButton("add_scenario", "Add Row", 
                       class = "btn btn-primary btn-sm", 
                       style = "margin-right: 5px;"),
          actionButton("delete_scenario", "Delete Row", 
                       class = "btn btn-warning btn-sm")
      ),
      
      p("Enter response rates (0-1) for each basket under different scenarios:"),
      rHandsontableOutput("scenario_input_table"),
      
      div(class = "alert alert-info",
          style = "margin: 10px 0;",
          strong("Performance Tip:"), 
          " For faster results, use ≤1000 simulations for testing, ≤100 for quick checks."
      )
    )
  })
  
  output$scenario_input_table <- renderRHandsontable({
    req(scenarios_data())
    
    data <- scenarios_data()
    if(is.null(data) || nrow(data) == 0) return(NULL)
    
    # Use row names directly as row headers, no Scenario column in the table
    rhandsontable(data, 
                  width = "100%", 
                  height = "200px",
                  rowHeaders = rownames(data)) %>%
      hot_cols(colWidths = rep(80, ncol(data))) %>%
      hot_col(1:ncol(data), type = "numeric", format = "0.000")
  })
  
  observeEvent(input$add_scenario, {
    current_data <- scenarios_data()
    nB <- input$num_baskets
    
    if(is.null(current_data)) {
      current_data <- data.frame(matrix(NA, nrow = 0, ncol = nB))
      colnames(current_data) <- paste0("Basket_", 1:nB)
    }
    
    if(nrow(current_data) >= 20) {
      showNotification("Maximum 20 scenarios allowed", type = "warning")
      return()
    }
    
    # Create new row with default values
    new_row <- data.frame(matrix(0.15, nrow = 1, ncol = nB))
    colnames(new_row) <- colnames(current_data)
    
    # Add new row
    updated_data <- rbind(current_data, new_row)
    
    # Renumber all row names as S1, S2, S3...
    rownames(updated_data) <- paste0("S", 1:nrow(updated_data))
    
    scenarios_data(updated_data)
    
    showNotification("Scenario added successfully!", type = "message", duration = 2)
  })
  
  observeEvent(input$delete_scenario, {
    current_data <- scenarios_data()
    
    if(is.null(current_data) || nrow(current_data) <= 1) {
      showNotification("Cannot delete - need at least one scenario!", type = "warning")
      return()
    }
    
    # Remove last row
    updated_data <- current_data[-nrow(current_data), , drop = FALSE]
    
    # Renumber all row names as S1, S2, S3...
    rownames(updated_data) <- paste0("S", 1:nrow(updated_data))
    
    scenarios_data(updated_data)
    showNotification("Scenario deleted successfully!", type = "message", duration = 2)
  })
  
  observeEvent(input$scenario_input_table, {
    updated_data <- hot_to_r(input$scenario_input_table)
    current_data <- isolate(scenarios_data())
    
    # Handle row count changes from manual table editing
    if(nrow(updated_data) != nrow(current_data)) {
      # Renumber all rows as S1, S2, S3... IMMEDIATELY
      rownames(updated_data) <- paste0("S", 1:nrow(updated_data))
      
      # For new rows with NA, set default value
      if(any(is.na(updated_data))) {
        updated_data[is.na(updated_data)] <- 0.15
      }
      
      # Force immediate update
      isolate(scenarios_data(updated_data))
    } else {
      # Preserve existing row names
      rownames(updated_data) <- rownames(current_data)
      
      # Validate numeric data - skip validation if NA to speed up
      if(any(is.na(updated_data))) {
        return()
      }
      
      if(any(updated_data < 0, na.rm = TRUE) || any(updated_data > 1, na.rm = TRUE)) {
        return()
      }
      
      isolate(scenarios_data(updated_data))
    }
  }, priority = 1000)
  
  observeEvent(input$run_simulation, {
    sim_results(NULL)
    
    # Check if efficacy cutoffs have been calculated
    if(is.null(cutoff_values())) {
      showNotification("Please calculate efficacy cutoffs in the Trial Design tab first!", type = "error")
      return()
    }
    
    if(is.null(input$num_baskets) || is.null(input$num_simulations) || is.null(scenarios_data())) {
      showNotification("Please configure Trial Design first and fill in all required fields", type = "error")
      return()
    }
    
    if(input$num_baskets <= 0 || input$num_simulations <= 0) {
      showNotification("Number of baskets and simulations must be positive!", type = "error")
      return()
    }
    
    if(input$num_simulations > 1000) {
      showNotification(paste("Running", input$num_simulations, "simulations may take several minutes."), 
                       type = "warning", duration = 5)
    }
    
    nB <- input$num_baskets
    K <- ifelse(is.null(input$num_interim), 0, input$num_interim)
    
    trial_design_configured <- TRUE
    for(i in 1:nB) {
      if(is.null(input[[paste0("n_basket_", i)]])) {
        trial_design_configured <- FALSE
        break
      }
    }
    
    if(!trial_design_configured) {
      showNotification("Using default sample sizes (20 per basket) since Trial Design not configured", 
                       type = "message", duration = 3)
    }
    
    if (K == 0) {
      nDat <- numeric(nB)
      for(i in 1:nB) {
        if(trial_design_configured) {
          sample_size <- input[[paste0("n_basket_", i)]]
          if(is.null(sample_size) || is.na(sample_size) || sample_size <= 0) {
            showNotification(paste("Invalid sample size for basket", i, "in Trial Design"), type = "error")
            return()
          }
          nDat[i] <- sample_size
        } else {
          nDat[i] <- 20
        }
      }
      Nmat <- matrix(nDat, nrow = nB, ncol = 1)
    } else {
      Nmat <- matrix(NA, nrow = nB, ncol = K + 1)
      for (b in 1:nB) {
        for (k in 1:K) {
          if(trial_design_configured) {
            interim_size <- input[[paste0("ia_", b, "_", k)]]
            if(is.null(interim_size) || is.na(interim_size) || interim_size <= 0) {
              showNotification(paste("Invalid interim size for basket", b, "stage", k, "in Trial Design"), type = "error")
              return()
            }
            Nmat[b, k] <- interim_size
          } else {
            Nmat[b, k] <- 5
          }
        }
        if(trial_design_configured) {
          final_size <- input[[paste0("n_basket_", b)]]
          if(is.null(final_size) || is.na(final_size) || final_size <= 0) {
            showNotification(paste("Invalid final size for basket", b, "in Trial Design"), type = "error")
            return()
          }
          Nmat[b, K + 1] <- final_size
        } else {
          Nmat[b, K + 1] <- 20
        }
      }
    }
    
    scenario_df <- scenarios_data()
    if(is.null(scenario_df) || nrow(scenario_df) == 0) {
      showNotification("Please configure scenarios!", type = "error")
      return()
    }
    
    # Extract ORRs directly - no Scenario column to remove
    ORRs <- as.matrix(scenario_df)
    mode(ORRs) <- "numeric"
    
    if(any(is.na(ORRs)) || any(ORRs < 0) || any(ORRs > 1)) {
      showNotification("Invalid response rates detected!", type = "error")
      return()
    }
    
    if(!is.null(input$set_seed) && !is.na(input$set_seed)) {
      set.seed(input$set_seed)
    }
    
    # Pass row names as scenario names
    result <- run_simulation_core(Nmat, ORRs, input$num_simulations, input$set_seed, nB, K, rownames(scenario_df))
    if(!is.null(result)) {
      sim_results(result)
      showNotification("Simulation completed successfully!", type = "message")
    }
  })
  
  run_simulation_core <- function(Nmat, ORRs, num_simulations, seed_val, nB, K, scenario_names) {
    withProgress(message = 'Running simulation...', value = 0, {
      
      incProgress(0.1, detail = "Generating simulated data...")
      
      sim_data <- tryCatch({
        generate.data(N = Nmat, ORRs = ORRs,
                      ntrial = num_simulations, 
                      seed = seed_val)
      }, error = function(e) {
        showNotification(paste("Error generating data:", e$message), type = "error")
        return(NULL)
      })
      
      if(is.null(sim_data)) return(NULL)
      
      # Use p0 from Trial Design tab
      pnull <- rep(input$p0, nB)
      
      # Get the efficacy cutoffs from Trial Design
      Q <- cutoff_values()$Efficacy_Cutoff
      
      stopbounds <- NULL
      if (K > 0) {
        stopbounds <- matrix(-1, nrow = nB, ncol = K)
      }
      
      selected_method <- input$design_method
      if(is.null(selected_method)) selected_method <- "Independent"
      
      incProgress(0.5, detail = paste("Running", selected_method, "method..."))
      
      postinf <- NULL
      if(selected_method == "Independent"){
        postinf <- tryCatch({
          post.infer(object = sim_data, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = "Independent")
        }, error = function(e) {
          showNotification(paste("Error in", selected_method, ":", e$message), type = "error")
          return(NULL)
        })
      } else if(selected_method == "JSD"){
        epsilon <- ifelse(is.null(input$epsilon), 3, input$epsilon)
        tau <- ifelse(is.null(input$tau), 0.5, input$tau)
        postinf <- tryCatch({
          post.infer(object = sim_data, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = "JSD", tau = tau, epsilon = epsilon)
        }, error = function(e) {
          showNotification(paste("Error in", selected_method, ":", e$message), type = "error")
          return(NULL)
        })
      } else if(selected_method == "local-PP-GEB"){
        if (is.null(input$a) || is.null(input$delta)) {
          showNotification("Please set parameters 'a' and 'delta' for local-PP-GEB method", type = "error")
          return(NULL)
        }
        postinf <- tryCatch({
          post.infer(object = sim_data, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = "localPP", method = "GEB", a = input$a, delta = input$delta)
        }, error = function(e) {
          showNotification(paste("Error in", selected_method, ":", e$message), type = "error")
          return(NULL)
        })
      } else if(selected_method == "local-PP-PEB"){
        if (is.null(input$a) || is.null(input$delta)) {
          showNotification("Please set parameters 'a' and 'delta' for local-PP-PEB method", type = "error")
          return(NULL)
        }
        postinf <- tryCatch({
          post.infer(object = sim_data, pnull = pnull, stopbounds = stopbounds, 
                     beta.a0 = rep(input$b1, nB), beta.b0 = rep(input$b2, nB), 
                     ModelFit = "localPP", method = "PEB", a = input$a, delta = input$delta)
        }, error = function(e) {
          showNotification(paste("Error in", selected_method, ":", e$message), type = "error")
          return(NULL)
        })
      }
      
      if(is.null(postinf)) return(NULL)
      
      incProgress(0.8, detail = "Computing operating characteristics...")
      
      # Use the efficacy cutoffs from Trial Design (not recalculating)
      res <- get.weighted.power(postinf, Q = Q)
      
      # Calculate operating characteristics
      power_matrix <- matrix(NA, nrow = nrow(ORRs), ncol = nB)
      for(s in 1:nrow(ORRs)) {
        for(b in 1:nB) {
          power_matrix[s, b] <- mean(postinf$postprob[s, , b] > Q[b])
        }
      }
      
      incProgress(0.9, detail = "Formatting results...")
      
      # Use scenario_names passed from rownames
      result_data <- data.frame(
        Scenario = scenario_names,
        stringsAsFactors = FALSE
      )
      
      for(b in 1:nB) {
        result_data[[paste0("Basket ", b)]] <- round(power_matrix[, b], 3)
      }
      
      result_data$FDR <- round(res$ind.error.fdr, 3)
      result_data$FWER <- round(res$ind.error.fw, 3)
      result_data$Power <- round(res$ind.power.cdr, 3)
      result_data$CCR <- round(res$ind.power.ccr, 3)
      
      incProgress(1.0, detail = "Complete!")
      
      return(result_data)
    })
  }
  
  output$simulation_results_table1 <- renderTable({
    req(sim_results())
    
    results_data <- sim_results()
    
    if(is.null(results_data) || nrow(results_data) == 0) {
      return(data.frame(Message = "No results available"))
    }
    
    results_display <- results_data
    
    basket_cols <- grep("^Basket", names(results_display))
    summary_cols <- c("FDR", "FWER", "Power", "CCR")
    
    for(col in c(basket_cols, which(names(results_display) %in% summary_cols))) {
      if(col <= ncol(results_display)) {
        results_display[, col] <- ifelse(is.na(results_display[, col]), 
                                         "", 
                                         sprintf("%.3f", as.numeric(results_display[, col])))
      }
    }
    
    return(results_display)
  }, caption = "Simulation Results", 
  caption.placement = "top", 
  striped = TRUE, 
  hover = TRUE)
  
  output$simulation_results_table2 <- renderUI({
    req(sim_results())
    
    div(style = "margin: 20px 0;",
        downloadButton("download_sim_results", "Download Results CSV", 
                       class = "btn btn-success")
    )
  })
  
  output$download_sim_results <- downloadHandler(
    filename = function() {
      paste0("basket_trial_simulation_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(sim_results())
      write.csv(sim_results(), file, row.names = FALSE)
    }
  )
  
  output$simulation_progress <- renderUI({
    if(!is.null(sim_results())) {
      div(class = "alert alert-success", 
          style = "margin-top: 20px;",
          "Simulation completed successfully! Results are displayed above.")
    }
  })
  
  # ------------------------
  # Trial Analysis
  # ------------------------
  
  output$trial_analysis_table_inputs <- renderUI({
    nB <- input$num_baskets_analysis
    lapply(1:nB, function(b) {
      fluidRow(
        column(6, numericInput(paste0("n_subj_", b),
                               paste("Subjects in Basket", b),
                               value = 20, min = 1)),
        column(6, numericInput(paste0("n_resp_", b),
                               paste("Responses in Basket", b),
                               value = 5, min = 0))
      )
    })
  })
  
  trial_analysis_results <- reactiveVal(NULL)
  
  observeEvent(input$analyze_trial, {
    B <- input$num_baskets_analysis
    req(B)
    
    nDat <- sapply(1:B, function(i) input[[paste0("n_subj_", i)]])
    yDat <- sapply(1:B, function(i) input[[paste0("n_resp_", i)]])
    
    if(any(yDat > nDat)) {
      showNotification("Number of responses cannot exceed number of subjects!", type = "error")
      return()
    }
    
    method <- input$design_method
    fitFun <- switch(method,
                     "Independent" = Independent,
                     "JSD" = JSD,
                     "local-PP-PEB" = localPP,
                     "local-PP-GEB" = localPP)
    
    fit <- if(method %in% c("Independent", "JSD")) {
      if (method == "JSD") {
        fitFun(nDat, yDat, be.a0 = rep(input$b1, B), be.b0 = rep(input$b2, B),
               epsilon = ifelse(is.null(input$epsilon), 2, input$epsilon),
               tau = ifelse(is.null(input$tau), 0.3, input$tau))
      } else {
        fitFun(nDat, yDat, be.a0 = rep(input$b1, B), be.b0 = rep(input$b2, B))
      }
    } else {
      fitFun(nDat, yDat,
             be.a0 = rep(input$b1, B),
             be.b0 = rep(input$b2, B),
             a = ifelse(is.null(input$a), 1, input$a),
             delta = ifelse(is.null(input$delta), 0.4, input$delta),
             method = ifelse(method=="local-PP-PEB","PEB","GEB"))
    }
    
    # Calculate posterior probability that p > p0
    p0 <- input$p0
    post_prob_greater_p0 <- pbeta(p0, fit$a.post, fit$b.post, lower.tail = FALSE)
    
    efficacy_cutoffs <- if(!is.null(cutoff_values())) {
      cutoff_values()$Efficacy_Cutoff[1:B]
    } else {
      rep(0.5, B)
    }
    
    trial_analysis_results(data.frame(
      Basket = paste("Basket", 1:B),
      `Posterior Prob of p > p0` = round(post_prob_greater_p0, 3),
      Efficacy_Cutoff = round(efficacy_cutoffs, 3),
      Promising = ifelse(post_prob_greater_p0 > efficacy_cutoffs, "Yes", "No"),
      check.names = FALSE
    ))
  })
  
  output$trial_analysis_results <- renderTable({
    trial_analysis_results()
  }, digits = 3)
}