
fitfun.basic <- function(X, y, 
                         blocks = rep(1, ncol(X)),
                         weights = rep(1, nrow(X)),
                         foldid, 
                         ...) {
  
  require(glmnet)
  require(dplyr)
  
  # clinical models
  t1 <- Sys.time()
  s  <- blocks == 1
  cvs  <- cv.glmnet(X[,s], y, family = "binomial",
                    foldid = foldid,
                    weights = weights)
  clin <- glmnet(X[,s], y, family = "binomial",
                 lambda = cvs$lambda.min,
                 weights = weights)
  t2 <- Sys.time()
  # 
  beta <- coef(clin)[,1]
  beta <- beta[beta != 0]
  clin <- list(beta = beta, time = t2 - t1)
  
  # mol
  t1 <- Sys.time()
  s  <- blocks == 2
  cvs <- cv.glmnet(X[,s], y, family = "binomial",
                   foldid = foldid,
                   weights = weights)

  mol <- glmnet(X[,s], y, family = "binomial",
                lambda = cvs$lambda.min,
                weights = weights)

  t2 <- Sys.time()
  beta <- coef(mol)[,1]
  beta <- beta[beta != 0]

  mol <- list(beta = beta,
              time = t2 - t1)
   
  # naive
  t1 <- Sys.time()
  cvs <- cv.glmnet(X, y, family = "binomial",
                   foldid = foldid,
                   weights = weights)
  naive <- glmnet(X, y, family = "binomial",
                  lambda = cvs$lambda.min,
                  weights = weights)
  t2 <- Sys.time()
  beta <- coef(naive)[,1]
  beta <- beta[beta != 0]
  naive <- list(beta = beta,
                time = t2 - t1)
  
  # 
  t1 <- Sys.time()
  cv.ipf <- cv.ipflas(X, y, family = "binomial",
                      foldid = foldid,
                      blocks = blocks,
                      weights = weights,
                      alpha = 1,
                      penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)))

  ipf.las <- ipflas(X, y, family = "binomial",
                    weights = weights,
                    alpha = 1,
                    blocks = blocks,
                    pf = cv.ipf$pf.min,
                    lambda = cv.ipf$lambda.min)
  t2 <- Sys.time()

  beta <- coef(ipf.las)[,1]
  beta <- beta[beta != 0]

  ipf.las <- list(beta = beta,
                  pf   = cv.ipf$pf.min,
                  time = t2 - t1)
  
  re <- list(clin = clin,
             mol = mol,
             naive = naive,
             ipf.las = ipf.las)
  return(re)
  
}
  
fitfun.stk <- function(X, y, 
                       blocks = rep(1, ncol(X)),
                       weights = rep(1, nrow(X)),
                       foldid,
                       foldid.internal = NULL,
                       ...) {
  
  # clinical models used to test the equality of fitfun.stk
  t1 <- Sys.time()
  s  <- blocks == 1
  cvs  <- cv.glmnet(X[,s], y, family = "binomial",
                    foldid = foldid,
                    weights = weights)
  clin <- glmnet(X[,s], y, family = "binomial",
                 lambda = cvs$lambda.min,
                 weights = weights)
  t2 <- Sys.time()
  beta <- coef(clin)[,1]
  beta <- beta[beta != 0]
  clin <- list(beta = beta, time = t2 - t1)
  
  # stk
  t1 <- Sys.time()
  stk.pf <- stk.bin(X, y,
                    foldid = foldid,
                    foldid.internal = foldid.internal,
                    weights = weights,
                    seed = NULL,
                    nfold = 10,
                    method = "bypf",
                    blocks = blocks,
                    track = TRUE,
                    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                    optimFun = c("lasso", "ridge", "BS"),
                    ifWeights = TRUE,
                    ifshrinkage = TRUE,
                    truncValue = 0,
                    transPFun = function(x) x)
  t2 <- Sys.time()
  stk.pf$time = t2 - t1
  
  #- 2
  t1 <- Sys.time()
  stk.comp <- stk.bin(X, y,
                      foldid = foldid,
                      foldid.internal = foldid.internal,
                      weights = weights,
                      seed = NULL,
                      nfold = 10,
                      method = "bycomp",
                      blocks = blocks,
                      track = TRUE,
                      penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                      optimFun = c("lasso", "ridge", "BS"),
                      ifWeights = TRUE,
                      ifshrinkage = TRUE,
                      truncValue = 0,
                      transPFun = function(x) x)
  t2 <- Sys.time()
  stk.comp$time = t2 - t1
  
  re <- list(clin = clin,
             stk.pf = stk.pf,
             stk.comp = stk.comp)
  return(re)
}
  
fitfun.random.stk <- function(X, y, 
                              blocks = rep(1, ncol(X)),
                              weights = rep(1, nrow(X)),
                              foldid,
                              ...){
  
  # clinical models used to test the equality of fitfun.stk
  t1 <- Sys.time()
  s  <- blocks == 1
  cvs  <- cv.glmnet(X[,s], y, family = "binomial",
                    foldid = foldid,
                    weights = weights)
  clin <- glmnet(X[,s], y, family = "binomial",
                 lambda = cvs$lambda.min,
                 weights = weights)
  t2 <- Sys.time()
  beta <- coef(clin)[,1]
  beta <- beta[beta != 0]
  clin <- list(beta = beta, time = t2 - t1)
  
  # 
  sampleWeights <- ifelse(weights != 0, 1, 0)
  id <- which(weights != 0)
  nboot <- 100
  bootID <- list()
  for(j in seq_len(nboot)) {
    boot <- c()
    # two groups
    for(k in unique(y)) {
      ind <- which(y == k)
      tmp <- sample(ind, sum(sampleWeights[ind]), replace = TRUE, prob = sampleWeights[ind] / sum(sampleWeights[ind]))
      boot <- c(boot, tmp)
    }
    
    bootID[[j]] <- boot
  }
  

  t1 <- Sys.time()
  random.stk.pf <- random.stack.bin(X, y,
                                    bootID = bootID,
                                    randomP = table(blocks),
                                    nboot = 100,
                                    fastTune = 20,
                                    foldid = NULL,
                                    nfold = 10,
                                    blocks = blocks,
                                    weights = weights,
                                    method = "bypf",
                                    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                                    optimFun = c("lasso", "ridge", "BS"),
                                    limits = 0,
                                    ifWeights = TRUE,
                                    ifshrinkage = TRUE,
                                    truncValue = 0,
                                    transPFun = function(x) x,
                                    ifRandomSL = FALSE)
  t2 <- Sys.time()
  random.stk.pf$time <- t2 - t1
  
  t1 <- Sys.time()
  random.stk.comp <- random.stack.bin(X, y,
                                      bootID = bootID, 
                                      randomP = table(blocks), 
                                      nboot = 100,
                                      fastTune = 20,
                                      foldid = NULL,
                                      nfold = 10,
                                      blocks = blocks,
                                      weights = weights,
                                      method = "bycomp",
                                      penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                                      optimFun = c("lasso", "ridge", "BS"),
                                      limits = 0,
                                      ifWeights = TRUE,
                                      ifshrinkage = TRUE, 
                                      truncValue = 0, 
                                      transPFun = function(x) x,
                                      ifRandomSL = FALSE)
  t2 <- Sys.time()
  random.stk.comp$time <- t2 - t1
  
  re <- list(clin = clin,
             random.stk.pf = random.stk.pf,
             random.stk.comp = random.stk.comp)
  return(re)
}



