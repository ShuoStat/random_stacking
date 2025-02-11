
source("./helper.R")
source("./stk.glm_v5.R")
source("./stk.cox_v6.R")
source("./ipflasso.R")
source(".random.stk_v4.R")

# load package
pkgs <- c("survival", "dplyr", "glmnet", "doParallel", "doSNOW", "ggplot2", "tidyr", "tidyverse", "ggpubr", "WeightedROC")
loadPackages(pkgs)

#-------------------------------------------------------------------------------
#- run it
#- Load data
#-------------------------------------------------------------------------------

nsim   <- 10
nfolds <- 10
nCores <- 25
datNames <- c("BRCA", "HNSC", "KIRC", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "PAAD", "SKCM")

#-------------------------------------------------------------------------------
#- Get best cutoff for each data
#-------------------------------------------------------------------------------

bestCuts <- c()

for(nam in datNames) {
  
  cat("\n\n", nam, ":", as.character(Sys.time()), "\n\n")
  #- get data
  getData(dir = "../data", nam, log2 = F, toBinary  = T, cutTimes = c(1, 2, 3, 5))
  bestCuts <- c(bestCuts, bestCut)
  
}

names(bestCuts) <- datNames

#-------------------------------------------------------------------------------

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
                                    randomP = c(sum(blocks == 1), sum(blocks == 2)  / 3),
                                    nboot = 200,
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
                                      randomP = c(sum(blocks == 1), sum(blocks == 2)  / 3), 
                                      nboot = 200,
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

#-------------------------------------------------------------------------------
#- Prediction
#-------------------------------------------------------------------------------

IPCW = TRUE

for(nam in datNames) {
  
  cat("\n\n", nam, ":", as.character(Sys.time()), "\n\n")

  #- get data
  getData(dir = "../data/" ,nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  n <- nrow(X)
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = T)
  
  # 
  # X <- X[, 1:500]
  # blocks <- blocks[1:500]
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  foreach(i = 1 : (nfolds * nsim), 
          .packages = c("glmnet", "survival"), 
          .options.snow = opts) %dopar% {
            
            samID <- ceiling(i / nfolds)
            #- j, which fold of the i.nsim split being the test
            j <- i - (samID - 1) * nfolds
            #- foldid for i.nsim splits and j fold
            
            sam <- sampleSplits[[samID]]
            foldid <- foldids[[samID]][[j]]
            foldid.internal <- foldid.internals[[samID]][[j]]
            X.t <- X[sam != j, ]
            y.t <- y[sam != j]
            X.v <- X[sam == j, ]
            
            # add IPCW weights
            if (IPCW) {
              weights <- getIPCW(matrix(1, nrow = nrow(X.t), ncol = 1), 
                                 matrix(1, nrow = nrow(X.t), ncol = 1),
                                 y.t, y.t, 
                                 bestCut = bestCuts[nam] * 365)
            } else {
              weights <- rep(1, nrow(X.t))
              zero_weight <- y.t[,1] <=  bestCuts[nam] * 365 & y.t[,2] == 0
              weights[zero_weight] <- 0
            }
            
            y.t <- y.t[ ,1] < bestCuts[nam] * 365
            
            args <- list(X = X.t, y = y.t, 
                         blocks = blocks,
                         weights = weights,
                         foldid = foldid,
                         foldid.internal = foldid.internal)
            
            mod <- do.call(fitfun.random.stk, args)
            save(list = "mod", file = paste0("../output/random.stk_", nam, "_", samID, "_", j, "_200B.RData"))
          }
  stopCluster(cl)
}

#-------------------------------------------------------------------------------


