source("./helper.R")
source("./fit.fun.R")
source("./stk.glm_v5.R")
source("./stk.cox_v6.R")
source("./ipflasso.R")
source("./random.stk_v4.R")

# load package
# List of required packages
pkgs <- c("survival", "dplyr", "glmnet", "doParallel", "doSNOW", "ggplot2", "tidyr", "tidyverse", "ggpubr", "WeightedROC")

loadPackages(pkgs)

#-------------------------------------------------------------------------------
#- Parameter setting
#- 
#-------------------------------------------------------------------------------

nsim   <- 10 # number of replicaes
nfolds <- 10 # number of folds for cross-validation
nCores <- 25 # number of threads/cores to be used for parallel processing 
# used data
datNames <- c("BRCA", "HNSC", "KIRC", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "PAAD", "SKCM")

# loop over each dataset in datNames
for(nam in datNames) {
  
  # print dataset name and timestamp
  cat("\n\n", nam, ":", as.character(Sys.time()), "\n\n")
  
  #- get data, output X, y, and blocks to global environment
  getData(dir = "../data/" ,nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  n <- nrow(X) # sample size
  # splits
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = T) 
  
  # X <- X[, 1:500]
  # blocks <- blocks[1:500]
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  # initialize progress bar for parallel execution
  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # parallel running
  foreach(i = 1 : (nfolds * nsim), 
          .packages = c("glmnet", "survival"), 
          .options.snow = opts) %dopar% {
            
            samID <- ceiling(i / nfolds)
            #- j, which fold of the i.nsim split being the test
            j <- i - (samID - 1) * nfolds
            #- foldid for i.nsim splits and j fold
            
            # sample ID for training
            sam <- sampleSplits[[samID]]
            load(paste0("../../../output/", nam, "-", samID, "-", j, "-bi.basic_v6.RData"))
            save(list = "mod", file = paste0("../output/basic_", nam, "_", samID, "_", j, ".RData"))
            
            load(paste0("../../../output/", nam, "-", samID, "-", j, "-bi.stk_v6.RData"))
            save(list = "mod", file = paste0("../output/stk_", nam, "_", samID, "_", j, ".RData"))
            
            load(paste0("../../../output/", nam, "-", samID, "-", j, "-bi.random.stk_v6.RData"))
            save(list = "mod", file = paste0("../output/random.stk_", nam, "_", samID, "_", j, ".RData"))
            
          }
  
  stopCluster(cl)
}



# loop over each dataset in datNames
for(nam in datNames) {
  
  # print dataset name and timestamp
  cat("\n\n", nam, ":", as.character(Sys.time()), "\n\n")
  
  #- get data, output X, y, and blocks to global environment
  getData(dir = "../data/" ,nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  n <- nrow(X) # sample size
  # splits
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = T) 
  
  # X <- X[, 1:500]
  # blocks <- blocks[1:500]
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  # initialize progress bar for parallel execution
  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # parallel running
  foreach(i = 1 : (nfolds * nsim), 
          .packages = c("glmnet", "survival"), 
          .options.snow = opts) %dopar% {
            
            samID <- ceiling(i / nfolds)
            #- j, which fold of the i.nsim split being the test
            j <- i - (samID - 1) * nfolds
            #- foldid for i.nsim splits and j fold
            
            # sample ID for training
            sam <- sampleSplits[[samID]]

            load(paste0("../../../output/", nam, "-", samID, "-", j, "-bi.random.stk_v7.RData"))
            save(list = "mod", file = paste0("../output/random.stk_", nam, "_", samID, "_", j, "_P1_fast100.RData"))
          }
  
  stopCluster(cl)
}















