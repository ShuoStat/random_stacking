
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

#-------------------------------------------------------------------------------
#- Clinical information, Table 1
#-------------------------------------------------------------------------------

clinTab  <- c()
bestCuts <- c()

for(nam in datNames) {
  
  cat("\n\n", nam, ":", as.character(Sys.time()), "\n\n")
  #- get data
  getData(dir = "../data", nam, log2 = F, toBinary  = T, cutTimes = c(1, 2, 3, 5))
  ncase   <- nrow(X)
  nevents <- sum(y == T)
  nclin   <- length(colnames(X)[blocks == 1])
  nmol    <- length(colnames(X)[blocks == 2])
  bestCuts <- c(bestCuts, bestCut)
  
  clinTab <- rbind(clinTab, c(nam, ncase, nevents, nclin, nmol, bestCut))
}

names(bestCuts) <- datNames
colnames(clinTab) <- c("Data", "Sample size", "N(high risks)", "N(clinical)", "N(mol)", "cutTime")
print(xtable::xtable(clinTab), include.rownames = FALSE)
write.csv(clinTab, "../results/clinTab.csv")

#-------------------------------------------------------------------------------
#- Training models
#-------------------------------------------------------------------------------

# indicator for whether to use Inverse Probability of Censoring Weights (IPCW)
IPCW = TRUE
# list of modeling methods to be used
method = c("basic", "stk", "random.stk")

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
            # foldid 
            foldid <- foldids[[samID]][[j]]
            # foldid for interal CV
            foldid.internal <- foldid.internals[[samID]][[j]]
            
            X.t <- X[sam != j, ] # training data X
            y.t <- y[sam != j]   # training data y
            X.v <- X[sam == j, ] # validation data X
            
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
            
            # get categorical  y
            y.t <- y.t[ ,1] < bestCuts[nam] * 365
            
            # argument list for model training function 
            args <- list(X = X.t, y = y.t, 
                         blocks = blocks,
                         weights = weights,
                         foldid = foldid,
                         foldid.internal = foldid.internal)
            
            if ("basic" %in% method) {
              mod <- do.call(fitfun.basic, args)
              save(list = "mod", file = paste0("../output/basic_", nam, "_", samID, "_", j, ".RData"))
            }
            
            if ("stk" %in% method) {
              mod <- do.call(fitfun.stk, args)
              save(list = "mod", file = paste0("../output/stk_", nam, "_", samID, "_", j, ".RData"))
            }
            
            if ("random.stk" %in% method) {
              mod <- do.call(fitfun.random.stk, args)
              save(list = "mod", file = paste0("../output/random.stk_", nam, "_", samID, "_", j, ".RData"))
            }
          }
  stopCluster(cl)
}

#-----------------------------------------------------------------------------
#- Evaluate their predictions 
#-------------------------------------------------------------------------------

# get prediction for validation data
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
  
  # index for stk, basic, and random.stk
  indStk <- grepl("^stk", names(obj))
  indRanStk <- grepl("^random.stk", names(obj))
  indBasic <- !(indStk | indRanStk)
  
  # generate results
  re <- list()
  if (any(indBasic)) {
    reBasic<- mapply(predFun, m = obj[indBasic], 
                     MoreArgs = list(newdata = X.v, type = type), 
                     SIMPLIFY = FALSE)
    re <- c(re, reBasic)
  }
  
  # prediction of stacking
  if (any(indStk)) {
    reStk <- mapply(predict.stk.bin, 
                    obj = obj[indStk],
                    MoreArgs = list(newdata = X.v, type = type, ifCVMods = ifCVMods),
                    SIMPLIFY = FALSE)
    re <- c(re, unlist(reStk, recursive = FALSE))
  }
  
  if (any(indRanStk)) {
    reRanStk <- mapply(predict.random.stk.bin, 
                       obj = obj[indRanStk],
                       MoreArgs = list(newdata = X.v, type = type),
                       SIMPLIFY = FALSE)
    re <- c(re, unlist(reRanStk, recursive = FALSE))
  }
  return(re)
}


aucs <- list() # store AUC values for each dataset
devs <- list() # store deviance values for each dataset
briers <- list() # store Brier score values for each dataset
IPCW = FALSE # disable IPCW, not necessary for validation

for(nam in datNames){
  
  print(nam)
  print(Sys.time())
  getData(dir = "../data/", nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  n <- nrow(X) # sample size 
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = TRUE) # train-validation splits for cross-validation
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  re <- foreach(i = 1 : (nfolds * nsim),
                .packages = c("glmnet", "survival", "dplyr")) %dopar% {
                  
                  samID <- ceiling(i / nfolds)
                  #- j, which fold of the i.nsim split being the test
                  j <- i - (samID - 1) * nfolds
                  #- foldid for i.nsim splits and j fold
                  
                  sam <- sampleSplits[[samID]]
                  X.v <- X[sam == j,]
                  y.v <- y[sam == j]
                  
                  X.t <- X[sam != j,]
                  y.t <- y[sam != j]
                  
                  # load selected trained models
                  if ("basic" %in% method) {
                    load(paste0("../output/basic_", nam, "_", samID, "_", j, ".RData"))
                    basic <- mod
                  }
                  
                  if ("stk" %in% method) {
                    load(paste0("../output/stk_", nam, "_", samID, "_", j, ".RData"))
                    stk <- mod
                  }
                  
                  if ("random.stk" %in% method) {
                    load(paste0("../output/random.stk_", nam, "_", samID, "_", j, ".RData"))
                    random.stk <- mod
                  }
                  
                  # text equality of clin models
                  if (i == 1) {
                    identical(basic$clin$beta, stk$clin$beta)
                    identical(basic$clin$beta, random.stk$clin$beta)
                  }
                  
                  # Weights
                  if (IPCW) {
                    
                    weights <- getIPCW(matrix(1, nrow = nrow(X.t), ncol = 1), 
                                       matrix(1, nrow = nrow(X.v), ncol = 1),
                                       y.t, y.v, 
                                       bestCut = bestCuts[nam] * 365)
                    
                    y.v <- y.v[ ,1] < bestCuts[nam] * 365
                    
                  } else {
                    grp <- biY(y.v, cutTimes = bestCuts[nam] * 365)
                    X.v <- X.v[grp$sel,]
                    weights <- rep(1, nrow(X.v))
                    y.v <- grp$biY
                  }
                  
                  # Occasionally, y.v can be all T or F
                  if (length(table(y.v)) == 1) {
                    return(NULL)
                  } else {
                    # probabilities
                    probs.basic <- predMods(X.v, basic, type = "response")
                    probs.stk   <- predMods(X.v, stk[-1], type = "response", ifCVMods = FALSE)
                    probs.randon.stk <- predMods(X.v, random.stk[-1], type = "response")
                    probs <- c(probs.basic, probs.stk, probs.randon.stk)
                    
                    # results
                    auc   <- lapply(probs, getAUC, y = y.v, weights = weights)
                    dev   <- lapply(probs, getDev, y = y.v, weights = weights)
                    brier <- lapply(probs, getBrier, y = y.v, weights = weights)
                    
                    return(list(auc = unlist(auc), 
                                dev = unlist(dev),
                                brier = unlist(brier)))
                  }
                }
  
  stopCluster(cl)
  
  aucs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "auc"))
  devs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "dev"))
  briers[[nam]] <- Reduce(rbind, lapply(re, `[[`, "brier"))
}
          
exportFun <- function(objList, filename, ifCI = FALSE) {
  
  # objList, one of the aucs, devs, and briers 
  # filename, name of the output file
  # ifCI, indicate if export confidence error
  
  out <- cbind.data.frame(lapply(objList, colMeans)) %>%
    t() %>%
    as.data.frame() %>%
    select(!contains(c("lasso", "BS"))) %>%
    mutate_all( ~ formatC(., format = "f", digits = 3))

  n <- nrow(objList[[1]]) # number of duplicated runs
  
  if (ifCI){
    se <- cbind.data.frame(lapply(objList, function(x) apply(x, 2, function(xx) sd(xx) / sqrt(n)))) %>%
      t() %>%
      as.data.frame() %>%
      select(!contains(c("lasso", "BS"))) %>%
      mutate_all( ~ formatC(., format = "f", digits = 3))
    
    out <- mapply(function(x, y){
      paste0(x, "(", y, ")")
    }, x = out, y = se)
  }
  
  out <- cbind(Data = names(objList), out) %>% as.data.frame()
  openxlsx::write.xlsx(out, file = paste0("../results/random_", filename, ".xlsx"))
  return(out)
}

# export results
exportFun(aucs, "auc", ifCI = TRUE)
exportFun(devs, "dev", ifCI = TRUE)
exportFun(briers, "brier", ifCI = TRUE)

#-------------------------------------------------------------------------------
# Weight distribution
#-------------------------------------------------------------------------------

weights <- list()

for(nam in datNames){
  
  print(nam)
  print(Sys.time())
  getData(dir = "../data/", nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  
  n <- nrow(X)
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = TRUE)
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  re <- foreach(i = 1 : (nfolds * nsim),
                .packages = c("glmnet", "survival", "dplyr")) %dopar% {
                  
                  samID <- ceiling(i / nfolds)
                  #- j, which fold of the i.nsim split being the test
                  j <- i - (samID - 1) * nfolds
                  #- foldid for i.nsim splits and j fold
                  
                  # stk
                  load(paste0("../output/stk_", nam, "_", samID, "_", j, ".RData"))
                  stk <- mod[-1]
                  
                  # random stacking
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, ".RData"))
                  random.stk <- mod[-1]
                  
                  stk <- c(stk, random.stk)
                  lapply(stk, function(x) x$alphas$ridge[-1])
                }
  
  stopCluster(cl)
  
  for(i in seq_along(re)) {
    if (i == 1)
      re0 <- re[[1]]
    else
      re0 <- mapply(rbind, re0, re[[i]])
  }
  
  weights[[nam]] <-  re0
}

# make box plot
plot.box <- function(weights, sel = c("CP", "PF")){
  
  # sel, export the weights of Rstk(CP) or Rstk(PF)
  
  # obj = weights
  if (sel == "CP"){
    # output Rstk(CP)
    obj1 <-  lapply(weights, `[[`, "stk.comp")
    obj2 <-  lapply(weights, `[[`, "random.stk.comp")
  }
  
  if (sel == "PF"){
    # output Rstk
    obj1 <-  lapply(weights, `[[`, "stk.pf")
    obj2 <-  lapply(weights, `[[`, "random.stk.pf")
  }
  
  # obj1
  obj1 <- lapply(obj1, function(x){
    
    x %>% as.data.frame %>% 
      #- remove intercept
      tidyr::pivot_longer(everything(),
                          names_to = "mods",
                          values_to = "alpha") %>%
      mutate(method = "Basic Stacking")
  })
  
  obj1 <- bind_rows(obj1, .id = "Data")
  
  # obj2
  obj2 <- lapply(obj2, function(x){
    
    x %>% as.data.frame %>% 
      #- remove intercept
      tidyr::pivot_longer(everything(),
                          names_to = "mods",
                          values_to = "alpha") %>%
      mutate(method = "Random Stacking")
  })
  
  obj2 <- bind_rows(obj2, .id = "Data")
  
  obj <- rbind(obj1, obj2) %>% 
    mutate(x = paste0(mods, method),
           x = factor(x, levels = unique(x)))
  
  #- ggplot
  ggplot2::ggplot(obj, aes(x = x, y = alpha, fill = factor(method))) + 
    geom_boxplot(outliers = FALSE) + 
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2, color = "black") +
    theme_bw() + 
    scale_x_discrete(labels = rep(unique(obj$mods), times = 2)) +
    facet_wrap(vars(Data), nrow = 2) + 
    xlab("") +
    ylab("Weights") + 
    theme(legend.title = element_blank(), 
          legend.position = "bottom",
          legend.margin = margin(t = -15, b = 0, l = 0, r = 0, unit = "pt"))
}

# combination
p <- plot.box(weights, sel = "CP") +  scale_x_discrete(labels = c("Clin", "Mol", "Clin", "Mol"))
ggsave("../results/random.weights.comp.tiff", plot = p, width = 7, height = 5, units = "in", dpi = 300, compression = "lzw")

# pf

p <- plot.box(weights, sel = "PF") + 
  scale_x_discrete(labels = c("1-1", "1-2", "1-4", "1-8", "1-1", "1-2", "1-4", "1-8")) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave("../results/random.weights.pf.tiff", plot = p, width = 8, height = 5, units = "in", dpi = 300, compression = "lzw")


#-------------------------------------------------------------------------------
# evaluate the effects random feature size on prediction accuracy of random stacking
# 1/3 * P, results from addrun_P1_3.R
# 1/6 * P, results from addrun_P1_6.R
# sqrt(p), results from addrun_P_sqrt.R
# 2/3 * P, results from addrun_P2_3.R
#-------------------------------------------------------------------------------

aucs <- list()
devs <- list()
briers <- list()

IPCW = FALSE

for(nam in datNames){
  
  print(nam)
  print(Sys.time())
  getData(dir = "../data/", nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  
  n <- nrow(X)
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = TRUE)
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  re <- foreach(i = 1 : (nfolds * nsim),
                .packages = c("glmnet", "survival", "dplyr")) %dopar% {
                  
                  samID <- ceiling(i / nfolds)
                  #- j, which fold of the i.nsim split being the test
                  j <- i - (samID - 1) * nfolds
                  #- foldid for i.nsim splits and j fold
                  
                  sam <- sampleSplits[[samID]]
                  X.v <- X[sam == j,]
                  y.v <- y[sam == j]
                  
                  X.t <- X[sam != j,]
                  y.t <- y[sam != j]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, ".RData"))
                  names(mod) <- paste0(names(mod), "1_1")
                  rstk1_1 <- mod
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_1_3P.RData"))
                  names(mod) <- paste0(names(mod), "1_3")
                  rstk1_3 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_1_6P.RData"))
                  names(mod) <- paste0(names(mod), "1_6")
                  rstk1_6 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_sqrtP.RData"))
                  names(mod) <- paste0(names(mod), "_sqrt")
                  rstk_sqrt <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_2_3P.RData"))
                  names(mod) <- paste0(names(mod), "2_3")
                  rstk2_3 <- mod[-1]
                  
                  rstk <- c(rstk1_1, rstk2_3, rstk1_3, rstk1_6, rstk_sqrt)

                  # Weights
                  if (IPCW) {
                    weights <- getIPCW(matrix(1, nrow = nrow(X.t), ncol = 1), 
                                       matrix(1, nrow = nrow(X.v), ncol = 1),
                                       y.t, y.v, 
                                       bestCut = bestCuts[nam] * 365)
                    
                    y.v <- y.v[ ,1] < bestCuts[nam] * 365
                    
                  } else {
                    grp <- biY(y.v, cutTimes = bestCuts[nam] * 365)
                    X.v <- X.v[grp$sel,]
                    weights <- rep(1, nrow(X.v))
                    y.v <- grp$biY
                  }
                  
                  # Occasionally, y.v can be all T or F
                  if (length(table(y.v)) == 1) {
                    return(NULL)
                  } else {
                    probs <- predMods(X.v, rstk, type = "response")
                    
                    # results
                    auc   <- lapply(probs, getAUC, y = y.v, weights = weights)
                    dev   <- lapply(probs, getDev, y = y.v, weights = weights)
                    brier <- lapply(probs, getBrier, y = y.v, weights = weights)
                    
                    return(list(auc = unlist(auc), 
                                dev = unlist(dev),
                                brier = unlist(brier)))
                  }
                }
  
  stopCluster(cl)
  
  aucs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "auc"))
  devs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "dev"))
  briers[[nam]] <- Reduce(rbind, lapply(re, `[[`, "brier"))
}

# Visualization

plotData <- function(objList) {
  
  mean <- cbind.data.frame(lapply(objList, colMeans)) %>%
    t() %>%
    as.data.frame() %>%
    dplyr::select(!contains(c("lasso", "BS", "clin")))
  
  n <- nrow(objList[[1]])
  
  se <- cbind.data.frame(lapply(objList, function(x) apply(x, 2, function(xx) sd(xx) / sqrt(n)))) %>%
    t() %>%
    as.data.frame() %>%
    select(!contains(c("lasso", "BS", "clin")))
  
  mean <- mean %>%
    rownames_to_column("data") %>%
    pivot_longer(cols = -1, names_to = "method", values_to = "mean")
  
  se <- se %>%
    rownames_to_column("data") %>%
    pivot_longer(cols = -1, names_to = "method", values_to = "se")
  
  re <- dplyr::left_join(mean, se, by = c("data", "method")) %>%
    mutate(upper = mean + se,
           lower = mean - se)
  
  return(re)
}

pData <- function(DataPlot, select = c("Rstk(PF)", "Rstk(CP)"), xlab = "Random feature size", ylab = "") {
  
  DataPlot <- filter(DataPlot, method == select)
  pList <- list()
  
  for(i in unique(DataPlot$data)) {
    
    dat <- filter(DataPlot, data == i)
    
    pList[[i]] <- ggplot(dat, aes(x = size, y = mean, color = size))+ 
      geom_line(group = 1, color = 1) + 
      geom_pointrange(aes(ymin = lower, ymax = upper)) +
      # facet_wrap( ~ data, nrow = 2) +
      theme_bw() + 
      labs(title = i,
           x = xlab,
           y = ylab) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            legend.title = element_blank(),
            legend.position = "none")
  }
  
  ggpubr::ggarrange(plotlist = pList, nrow = 2, ncol = 5)
  
}

# auc, Rstk(CP)
DataPlot <- plotData(aucs)

DataPlot <- DataPlot %>% 
  mutate(size = case_when(str_detect(DataPlot$method, "1_1") ~ "1P",
                          str_detect(DataPlot$method, "2_3") ~ "2/3P",
                          str_detect(DataPlot$method, "1_3") ~ "1/3P",
                          str_detect(DataPlot$method, "1_6") ~ "1/6P",
                          str_detect(DataPlot$method, "sqrt") ~ "sqrt(P)"),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("1P", "2/3P", "1/3P", "1/6P", "sqrt(P)")))

p <- pData(DataPlot, select = "Rstk(CP)", ylab = "AUC")
ggsave(filename = paste0("../results/featureSizeAUC_comp.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)


# auc, Rstk(PF)
DataPlot <- plotData(aucs)

DataPlot <- DataPlot %>% 
  mutate(size = case_when(str_detect(DataPlot$method, "1_1") ~ "1P",
                          str_detect(DataPlot$method, "2_3") ~ "2/3P",
                          str_detect(DataPlot$method, "1_3") ~ "1/3P",
                          str_detect(DataPlot$method, "1_6") ~ "1/6P",
                          str_detect(DataPlot$method, "sqrt") ~ "sqrt(P)"),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("1P", "2/3P", "1/3P", "1/6P", "sqrt(P)")))

p <- pData(DataPlot, select = "Rstk(PF)", ylab = "AUC")
ggsave(filename = paste0("../results/featureSizeAUC_pf.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

# dev, Rstk(CP)
DataPlot <- plotData(devs)

DataPlot <- DataPlot %>% 
  mutate(size = case_when(str_detect(DataPlot$method, "1_1") ~ "1P",
                          str_detect(DataPlot$method, "2_3") ~ "2/3P",
                          str_detect(DataPlot$method, "1_3") ~ "1/3P",
                          str_detect(DataPlot$method, "1_6") ~ "1/6P",
                          str_detect(DataPlot$method, "sqrt") ~ "sqrt(P)"),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("1P", "2/3P", "1/3P", "1/6P", "sqrt(P)")))

p <- pData(DataPlot, select = "Rstk(CP)", ylab = "Dev")
ggsave(filename = paste0("../results/featureSizeDev_comp.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)


# dev, Rstk(PF)
DataPlot <- plotData(devs)

DataPlot <- DataPlot %>% 
  mutate(size = case_when(str_detect(DataPlot$method, "1_1") ~ "1P",
                          str_detect(DataPlot$method, "2_3") ~ "2/3P",
                          str_detect(DataPlot$method, "1_3") ~ "1/3P",
                          str_detect(DataPlot$method, "1_6") ~ "1/6P",
                          str_detect(DataPlot$method, "sqrt") ~ "sqrt(P)"),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("1P", "2/3P", "1/3P", "1/6P", "sqrt(P)")))

p <- pData(DataPlot, select = "Rstk(PF)", ylab = "AUC")
ggsave(filename = paste0("../results/featureSizeDev_pf.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# Evaluate the effects of number of bootstrap, using 1/3 P
# B = 200
# B = 150
# B = 100
# B = 50
#-------------------------------------------------------------------------------

aucs <- list()
devs <- list()
briers <- list()

IPCW = FALSE

for(nam in datNames){
  
  print(nam)
  print(Sys.time())
  getData(dir = "../data/", nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  
  n <- nrow(X)
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = TRUE)
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  re <- foreach(i = 1 : (nfolds * nsim),
                .packages = c("glmnet", "survival", "dplyr")) %dopar% {
                  
                  samID <- ceiling(i / nfolds)
                  #- j, which fold of the i.nsim split being the test
                  j <- i - (samID - 1) * nfolds
                  #- foldid for i.nsim splits and j fold
                  
                  sam <- sampleSplits[[samID]]
                  X.v <- X[sam == j,]
                  y.v <- y[sam == j]
                  
                  X.t <- X[sam != j,]
                  y.t <- y[sam != j]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_50B.RData"))
                  names(mod) <- paste0(names(mod), "50")
                  rstk50 <- mod
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_100B.RData"))
                  names(mod) <- paste0(names(mod), "100")
                  rstk100 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_150B.RData"))
                  names(mod) <- paste0(names(mod), "150")
                  rstk150 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_200B.RData"))
                  names(mod) <- paste0(names(mod), "200")
                  rstk200 <- mod[-1]
                  
                  rstk <- c(rstk50, rstk100, rstk150, rstk200)
                  
                  # Weights
                  if (IPCW) {
                    weights <- getIPCW(matrix(1, nrow = nrow(X.t), ncol = 1), 
                                       matrix(1, nrow = nrow(X.v), ncol = 1),
                                       y.t, y.v, 
                                       bestCut = bestCuts[nam] * 365)
                    
                    y.v <- y.v[ ,1] < bestCuts[nam] * 365
                    
                  } else {
                    grp <- biY(y.v, cutTimes = bestCuts[nam] * 365)
                    X.v <- X.v[grp$sel,]
                    weights <- rep(1, nrow(X.v))
                    y.v <- grp$biY
                  }
                  
                  # Occasionally, y.v can be all T or F
                  if (length(table(y.v)) == 1) {
                    return(NULL)
                  } else {
                    probs <- predMods(X.v, rstk, type = "response")
                    
                    # results
                    auc   <- lapply(probs, getAUC, y = y.v, weights = weights)
                    dev   <- lapply(probs, getDev, y = y.v, weights = weights)
                    brier <- lapply(probs, getBrier, y = y.v, weights = weights)
                    
                    return(list(auc = unlist(auc), 
                                dev = unlist(dev),
                                brier = unlist(brier)))
                  }
                }
  
  stopCluster(cl)
  
  aucs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "auc"))
  devs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "dev"))
  briers[[nam]] <- Reduce(rbind, lapply(re, `[[`, "brier"))
}

# auc
DataPlot <- plotData(aucs)
DataPlot <- DataPlot %>% 
  mutate(size = str_extract(DataPlot$method, "50|100|150|200"),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("50", "100", "150", "200")))

p <- pData(DataPlot, select = "Rstk(CP)", xlab = "Number of bootstraps", ylab = "")
ggsave(filename = paste0("../results/bootstrapSizeAUC_cp.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

p <- pData(DataPlot, select = "Rstk(PF)", xlab = "Number of bootstraps", ylab = "")
ggsave(filename = paste0("../results/bootstrapSizeAUC_pf.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

# dev
DataPlot <- plotData(devs)

DataPlot <- DataPlot %>% 
  mutate(size = str_extract(DataPlot$method, "50|100|150|200"),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("50", "100", "150", "200")))
p <- pData(DataPlot, select = "Rstk(CP)", xlab = "Number of bootstraps", ylab = "")
ggsave(filename = paste0("../results/bootstrapSizeDev_cp.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

p <- pData(DataPlot, select = "Rstk(PF)", xlab = "Number of bootstraps", ylab = "")
ggsave(filename = paste0("../results/bootstrapSizeDev_pf.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)


#-------------------------------------------------------------------------------
# The effects of lambda (fast random stacking) based on 1/3*P
# v tuning % = 100, v15
# v tuning % = 50, v16
# v tuning % = 20, v8
# v tuning % = 10, v17
#-------------------------------------------------------------------------------

tunings <- list()
IPCW = FALSE

aucs <- list()
devs <- list()
briers <- list()

for(nam in datNames){
  
  print(nam)
  print(Sys.time())
  getData(dir = "../data/", nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  
  n <- nrow(X)
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = TRUE)
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  re <- foreach(i = 1 : (nfolds * nsim),
                .packages = c("glmnet", "survival", "dplyr")) %dopar% {
                  
                  samID <- ceiling(i / nfolds)
                  #- j, which fold of the i.nsim split being the test
                  j <- i - (samID - 1) * nfolds
                  #- foldid for i.nsim splits and j fold
                  
                  sam <- sampleSplits[[samID]]
                  X.v <- X[sam == j,]
                  y.v <- y[sam == j]
                  
                  X.t <- X[sam != j,]
                  y.t <- y[sam != j]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_fast_100.RData"))
                  names(mod) <- paste0(names(mod), "_100")
                  rstk100 <- mod
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_fast_50.RData"))
                  names(mod) <- paste0(names(mod), "_50")
                  rstk50 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_fast_20.RData"))
                  names(mod) <- paste0(names(mod), "_20")
                  rstk20 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_fast_10.RData"))
                  names(mod) <- paste0(names(mod), "_10")
                  rstk10 <- mod[-1]
                  
                  rstk <- c(rstk100, rstk50, rstk20, rstk10)
                  
                  # Weights
                  if (IPCW) {
                    weights <- getIPCW(matrix(1, nrow = nrow(X.t), ncol = 1), 
                                       matrix(1, nrow = nrow(X.v), ncol = 1),
                                       y.t, y.v, 
                                       bestCut = bestCuts[nam] * 365)
                    
                    y.v <- y.v[ ,1] < bestCuts[nam] * 365
                    
                  } else {
                    grp <- biY(y.v, cutTimes = bestCuts[nam] * 365)
                    X.v <- X.v[grp$sel,]
                    weights <- rep(1, nrow(X.v))
                    y.v <- grp$biY
                  }
                  
                  # Occasionally, y.v can be all T or F
                  if (length(table(y.v)) == 1) {
                    return(NULL)
                  } else {
                    probs <- predMods(X.v, rstk, type = "response")
                    
                    # results
                    auc   <- lapply(probs, getAUC, y = y.v, weights = weights)
                    dev   <- lapply(probs, getDev, y = y.v, weights = weights)
                    brier <- lapply(probs, getBrier, y = y.v, weights = weights)
                    
                    return(list(auc = unlist(auc), 
                                dev = unlist(dev),
                                brier = unlist(brier)))
                  }
                }
  
  stopCluster(cl)
  
  aucs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "auc"))
  devs[[nam]] <-  Reduce(rbind, lapply(re, `[[`, "dev"))
  briers[[nam]] <- Reduce(rbind, lapply(re, `[[`, "brier"))
}

# auc
DataPlot <- plotData(aucs)
DataPlot <- DataPlot %>% 
  mutate(size = gsub("_", "", str_extract(DataPlot$method, "_100|_50|_20|_10")),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("100", "50", "20", "10")))

p <- pData(DataPlot, select = "Rstk(CP)", xlab = "% tunned parameters", ylab = "AUC")
ggsave(filename = paste0("../results/fastTuningAUC_cp.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

p <- pData(DataPlot, select = "Rstk(PF)", xlab = "% tunned parameters", ylab = "AUC")
ggsave(filename = paste0("../results/fastTuningAUC_pf.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

# dev
DataPlot <- plotData(devs)
DataPlot <- DataPlot %>% 
  mutate(size = gsub("_", "", str_extract(DataPlot$method, "_100|_50|_20|_10")),
         method = case_when(str_detect(DataPlot$method, "stk.comp") ~ "Rstk(CP)",
                            str_detect(DataPlot$method, "stk.pf") ~ "Rstk(PF)")) %>%
  mutate(size = factor(size, levels = c("100", "50", "20", "10")))

p <- pData(DataPlot, select = "Rstk(CP)", xlab = "% tunned parameters", ylab = "Deviance")
ggsave(filename = paste0("../results/fastTuningDev_cp.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

p <- pData(DataPlot, select = "Rstk(PF)", xlab = "% tunned parameters", ylab = "Deviance")
ggsave(filename = paste0("../results/fastTuningDev_pf.jpeg"), p, width = 12, height = 5, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# Computational complexity
# random stacking, P, 20% tuning, v6
# random stacking, P, 100% tuning, v7
# random stacking, 1/3P, 100% tuning v15
# random stacking, 1/3P, 20% tuning v8
# basic stacking
# IPFLASSO
#-------------------------------------------------------------------------------

times <- list()
for(nam in datNames){
  
  print(nam)
  print(Sys.time())
  getData(dir = "../data/", nam, log2 = T, toBinary  = F, cutTimes = c(1, 2, 3, 5))
  
  n <- nrow(X)
  getSplits(n, nsim, nfolds, seed = 6278392, foldid.internal = TRUE)
  
  #- parallel runs
  nCores <- pmin(detectCores() - 2, nCores)
  cl <- makeCluster(nCores)
  registerDoSNOW(cl)
  
  re <- foreach(i = 1 : (nfolds * nsim),
                .combine = "rbind",
                .packages = c("glmnet", "survival", "dplyr")) %dopar% {
                  
                  samID <- ceiling(i / nfolds)
                  #- j, which fold of the i.nsim split being the test
                  j <- i - (samID - 1) * nfolds
                  #- foldid for i.nsim splits and j fold
                  
                  sam <- sampleSplits[[samID]]
                  X.v <- X[sam == j,]
                  y.v <- y[sam == j]
                  
                  X.t <- X[sam != j,]
                  y.t <- y[sam != j]
                  
                  load(paste0("../output/basic_", nam, "_", samID, "_", j, ".RData"))
                  basic <- mod
  
                  load(paste0("../output/stk_", nam, "_", samID, "_", j, ".RData"))
                  stk <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, ".RData"))
                  names(mod) <- paste0(names(mod), "1P_0.2T") # 100% features, 2.% tuning
                  rstk6 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_P1_fast100.RData"))
                  names(mod) <- paste0(names(mod), "1P_1T") 
                  rstk7 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_fast_100.RData"))
                  names(mod) <- paste0(names(mod), "0.3P_1T")
                  rstk15 <- mod[-1]
                  
                  load(paste0("../output/random.stk_", nam, "_", samID, "_", j, "_fast_20.RData"))
                  names(mod) <- paste0(names(mod), "0.3P_0.2T")
                  rstk8 <- mod[-1]
                  
                  rstk <- c(basic, stk, rstk6, rstk7, rstk15, rstk8)
                  # extract time
                  time <- lapply(rstk, `[[`, "time")
                  # to minutes
                  # Function to convert time differences to minutes
                  convert_to_minutes <- function(time_diff) {
                    if (attr(time_diff, "units") == "secs") {
                      as.numeric(time_diff) / 60
                    } else if (attr(time_diff, "units") == "mins") {
                      as.numeric(time_diff)
                    } else if (attr(time_diff, "units") == "hours") {
                      as.numeric(time_diff) * 60
                    } else {
                      stop("Unknown time unit")
                    }
                  }
                  
                  # Apply the conversion function to each element in the list
                  time <- lapply(time, convert_to_minutes)
                  return(unlist(time))
                }
  stopCluster(cl)
  
  times[[nam]] <-  re
}

#- 

getTime <- function(times, datName) {
  
  time <- colMeans(times[[datName]])
  time <- data.frame(method = names(time), time = time) %>%
    # remove cli, mol, naive
    filter(!method %in% c("clin", "mol", "naive")) %>%
    arrange(time)
  
    # duplicate ipf.las
  
  time <- filter(time, str_detect(method, ".pf|ipf")) %>%
    mutate(group = "PF") %>%
    bind_rows(
      filter(time, str_detect(method, ".comp|ipf")) %>%
        mutate(group = "CP")) 
  
  # generate labels
  time <- mutate(time, 
                 x = case_when(str_detect(method, "ipf.las") ~ "IPF(las)",
                               str_detect(method, "^stk") ~ "Basic Stk",
                               str_detect(method, "1P_0.2T") ~ "Rstk(P, 20%)",
                               str_detect(method, "1P_1T") ~ "Rstk(P, 100%)",
                               str_detect(method, "0.3P_1T") ~ "Rstk(1/3P, 100%)",
                               str_detect(method, "0.3P_0.2T") ~ "Rstk(1/3P, 20%)"),
                 x = factor(x, levels = unique(x)))
  return(time)
}

timeggData <- getTime(times, "BRCA")
p1 <- ggplot(timeggData, aes(x = x, y = time)) +
  geom_bar(stat = "identity", width = 0.5) +
  # geom_text(aes(label = time), vjust = -0.5) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 140)) +
  labs(x = NULL, title = "BRCA", y = "Time/min") + 
  coord_flip() +
  facet_grid( ~ group)

ggsave("../results/Rstk_time.jpeg", plot = p1, width = 7, height = 3, units = "in", dpi = 600)

#-------------------------------------------------------------------------------













