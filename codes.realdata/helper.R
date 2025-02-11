
# load package

loadPackages <- function(package_list) {
  for (package in package_list) {
    if (require(package, character.only = TRUE, quietly = TRUE)) {
      # Package is installed, load it
      message(paste("Package", package, "loaded successfully."))
    } else {
      # Package is NOT installed, install and then load
      message(paste("Package", package, "not found. Installing..."))
      tryCatch({
        install.packages(package, dependencies = TRUE)  # Install with dependencies
        if (require(package, character.only = TRUE)) {
          message(paste("Package", package, "installed and loaded successfully."))
        } else {
          stop(paste("Failed to load package", package, "after installation."))
        }
      }, error = function(e) {
        message(paste("Error installing package", package, ":", e$message))
        # You could choose to stop the script here with stop() if a critical package fails.
      })
    }
  }
  invisible(NULL)
}


#-------------------------------------------------------------------------------
#- prediction performance with deviance difference
#-------------------------------------------------------------------------------

predMods <- function(X.v, obj, type = "link", ifCVMods = FALSE){
  
  predFun <- function(m, newdata, type){
    beta <- m$beta
    # beta name is null for PCA
    if (!is.null(names(beta)))
      s <- names(beta)[-1]
    else
      s <- seq_along(beta[-1])
    
    pred <- cbind(1, X.v[,s, drop = FALSE]) %*% beta
    
    if (type == "response")
      pred <- 1 / (1 + exp(-pred))    
    
    return(pred)
  }
  
  indStk <- grepl("stk", names(obj))
  
  re.nonstk <- mapply(predFun, m = obj[!indStk], 
                      MoreArgs = list(newdata = X.v, type = type), 
                      SIMPLIFY = FALSE)
  
  # prediction of stacking
  
  re.stk <- mapply(predict.stk.bin, 
                   obj = obj[indStk],
                   MoreArgs = list(newdata = X.v, type = type, ifCVMods = ifCVMods),
                   SIMPLIFY = FALSE)
  
  # flatten stk prediction
  re.stk <- unlist(re.stk, recursive = FALSE)
  
  re <- c(re.nonstk, re.stk)
  return(re)
}

getAUC <- function(haty, y, weights = NULL) {
  
  # remove zero weights
  if (is.null(weights))
    weights <- rep(1, length(y))
  
  indZero <- weights == 0
  roc <- WeightedROC::WeightedROC(as.numeric(haty)[!indZero], 
                                  as.numeric(y)[!indZero], 
                                  weight = weights[!indZero])
  auc <- WeightedROC::WeightedAUC(roc)
  
  # roc <- pROC::roc(as.numeric(y)[!indZero], as.numeric(haty)[!indZero])
  return(as.numeric(auc))
}

getBrier <- function(haty, y, weights = NULL) {
  
  if (is.null(weights))
    weights <- rep(1, length(y))
  
  haty <- drop(haty)
  brier <- sum(weights * (haty - y)^2) / sum(weights)
  return(brier)
}

getDev <- function(haty, y, weights) {
  
  if (is.null(weights))
    weights <- rep(1, length(y))
  
  haty <- drop(haty)
  dev <- -2 * (sum(y * log(haty) * weights) + sum(((1- y)) * log((1 - haty)) * weights))
  return(dev)
}


# 
# #- AUC
# auc <- function(haty, y){
#   p <- 1 / (1 + exp(-haty))
#   roc <- pROC::roc(as.numeric(y), as.numeric(p))
#   as.numeric(roc$auc)
# }
# 
# #- Misclassification
# mis <- function(haty, y){
#   y <- as.numeric(y)
#   haty <- as.numeric(haty > 0.5)
#   sum(abs(y - haty)) / (length(y))
# }
# 
# #- deviance
# dev <- function(haty, y){
#   p <- 1 / (1 + exp(-haty))
#   p <- drop(p)
#   dev <- -2 * (sum(y * log(p)) + sum(((1- y)) * log((1 - p))))
#   dev
# }
# 
# 
# brier <- function(haty, y) {
#   p <- 1 / (1 + exp(-haty))
#   p <- drop(p)
#   brier <- mean((p - y)^2)
#   brier
# }
# 
# auc <- unlist(lapply(re, auc, y = y.v))
# mis <- unlist(lapply(re, mis, y = y.v))
# dev <- unlist(lapply(re, dev, y = y.v))
# brier <- unlist(lapply(re, brier, y = y.v))
# 
# 
# list(auc = auc, mis = mis, dev = dev, brier = brier)



# function to import data
# update, also return surv object for IPCW weights

getData <- function(dir, datName, log2 = T, toBinary = T,
                    cutTimes = NULL) {
  
  # datName, data name
  # log2, whether log2 transformation on TPM (molecular information)
  # cutTimes, a vector of years. If null, all times will be used
  
  dir <- paste0(dirname(dir), "/", basename(dir), "/")
  datDir  <- paste0(dir, datName, ".RData")
  load(datDir)
  
  y <- survival::Surv(X[ ,"time"], event = X[ ,"status"])
  X <- X[ ,-c(1:3)]
  
  b <- as.numeric(!grepl("clin", colnames(X))) + 1
  
  if (log2) {
    X[ ,b == 2] <- log2(X[ ,b == 2] + 0.001)
  }
  
  if (is.null(cutTimes))
    cutTimes <- sort(y[ ,1], decreasing = F)
  else
    cutTimes <- cutTimes * 365
  
  #- trans to binary response
  if (toBinary){
    
    # cutTime, from years to days
    obj <- biY(y, cutTimes = cutTimes)
    X <- X[obj$sel,]
    
    bestCut <- obj$bestCut
    if (!is.null(cutTimes))
      bestCut = bestCut / 365
    
    surv <- y[obj$sel]
    y <- obj$biY
  }
  
  #- Rename molecules since the original name is too long
  colnames(X) <- c(colnames(X)[b == 1], 
                   paste0("m", seq_len(sum(b == 2))))
  
  #- return
  .GlobalEnv$X <- as.matrix(X)
  .GlobalEnv$y <- y
  .GlobalEnv$blocks <- b
  if (toBinary) {
    .GlobalEnv$bestCut <- bestCut
    .GlobalEnv$surv <- surv 
  }
}


getSplits <- function(n, nsim, nfolds, seed, foldid.internal = T){
  
  # n, sample size
  # nsim, number of duplicated runs
  
  set.seed(seed)
  
  sampleSplits <- list()
  foldids <- list()
  
  for(i in seq_len(nsim)){
    samsplit <- sample(rep(1:10, length = n), n)
    sampleSplits[[i]] <- samsplit
    foldids[[i]] <- sapply(1:nfolds, 
                           function(id){ 
                             nlength = sum(samsplit != id)
                             sample(rep(1:10, length = nlength), nlength)}, 
                           simplify = F)
  }
  
  if (nsim == 1) {
    sampleSplits <- unlist(sampleSplits)
    foldids <- foldids[[1]]
  }
  
  if (foldid.internal) {
    
    foldid.internals <- list()
    for(i in seq_len(nsim)){
      tmp <- list()
      for(j in seq_along(foldids[[i]])){
        sams <- foldids[[i]][[j]]
        tmp[[j]] <- sapply(1:nfolds, 
                           function(x){ 
                             nlength = sum(sams != x)
                             sample(rep(1:10, length = nlength), nlength)
                           }, simplify = F)
      }
      foldid.internals[[i]] <- tmp
    }
    .GlobalEnv$foldid.internals <- foldid.internals
  }
  
  #- return
  .GlobalEnv$sampleSplits <- sampleSplits
  .GlobalEnv$foldids <- foldids
  
}


# Use IPCW weights

getIPCW <- function(X.t, X.v, surv.t, surv.v, bestCut) {
  
  censT <- Surv(surv.t[ ,1], 1 - surv.t[ ,2])
  censV <- Surv(surv.v[ ,1], 1 - surv.v[ ,2])
  # fit cox model
  mod  <- coxph(censT ~ ., data = as.data.frame(X.t))
  lp.t <- predict(mod, as.data.frame(X.t))
  lp.v <- predict(mod, as.data.frame(X.v))
  
  # high, high risk groups. 1, censoring
  s <- surv.v[,1] >= bestCut | surv.v[,2] == 1
  # if (any(!s)) {
  #   warning("Several early censored cases will have zero weights")
  # }
  # 
  # high risk group
  high <- surv.v[, 1] < bestCut
  
  # get times: high time + bestCut
  # time should be in increasing order
  times <- sort(c(censV[, 1], bestCut), decreasing = F)
  
  # get cumulative hazard
  chz <- cumhz(lp.t, censT, times = times)
  probs <- mapply(function(x) {exp(-x * exp(lp.v))}, x = chz, SIMPLIFY = T)
  
  # get probabilities
  timeUse <- ifelse(high, surv.v[, 1], bestCut)
  ind     <- match(timeUse, times)
  
  p <- c()
  for(i in seq_along(timeUse)) {
    p <- c(p, probs[i, ind[i]])
  }
  
  weights <- ifelse(s, 1 / p, 0)
  setNames(weights, rownames(X.v))
}

#-------------------------------------------------------------------------------

