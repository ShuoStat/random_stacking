################################################################################
# version 6
# gaussian method
# based on stk.glm_v6.R and above version
################################################################################

# history version
# versoin 1
# Uisng bootstrap, replace 10 CV re-sampling

# versoin 2
# stratified bootstrap for glm 

# version 3
# randomization in super learner

# version 4
# option for using randomization in super learner (only in glm)
# update random stacking in cox

# version 5
# allow customize the bootstrap aggregating function


################################################################################
# 1 cox section
################################################################################
#'
#' Base Learner for Stacked Cox 
#'
#' Stacking is a two-stage procedure. StkCoxBL fits the base learner of the survival stacking.
#'
#' @param X A matrix of predictor variables. Row, samples; Columns, variables
#' @param y Survival object.
#' @param times Quantiles of observed times to evaluate, default is NULL to calculate within function.
#' @param foldid Vector indicating fold IDs for cross-validation of stacking, auto-generated if NULL.
#' @param foldid.internal A list of internal fold IDs for nested cross-validation, auto-generated if NULL.
#' @param blocks Integer vector indicating block IDs for each predictor in X.
#' @param block.clin Index for clinical data block within blocks.
#' @param nfold Number of folds for cross-validation.
#' @param ifCox Logical indicating if partial likelihood optimization is used.
#' @param cox.time A specific time point for Cox model predictions, calculated if NULL.
#' @param method Methods to combine components, either "bycomp" or "bypf".
#' @param cens.method A method for handling censoring in the Brier score calculation, default "equal".
#' @param penalty.factors List of penalty factors if method = bypf.
#' @param aggFun bootstrap aggregation function, default is mean().
#' @param track Logical to track intermediate results.
#'
#' @return A list containing model fits, predictions, and possibly Cox predictions if `ifCox` is TRUE.
#' @examples
#'
#' @export random.StkCoxBL

random.StkCoxBL <- function(X, y, 
                            times = NULL,
                            bootID = NULL, 
                            randomP = c(ncol(X)), 
                            nboot = 100,
                            fastTune = 100, 
                            foldid = NULL,
                            nfold = 10,
                            blocks = rep(1, ncol(X)), 
                            block.clin = 1,
                            ifCox = TRUE, 
                            cox.time = NULL,
                            method = c("bycomp", "bypf"),
                            cens.method = "equal",
                            penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                            aggFun = mean,
                            ...) {
  
  require(dplyr)
  require(survival)
  
  # argument checking
  if (fastTune > nboot) {
    fastTune = nboot
    warning(paste0("fastTune is reduced to ", nboot))
  }
  
  # cens.method, a method for IBS censoring: "Breslow", "KM", "equal", "clinicalAdjusted", 
  # "medianSurv"
  # intercept, if intercept is required for lasso and ridge
  # method, either c("bycomp", "bypf") or function. When method is a function, X and y should be the mandatory arguments, and the output should be a list of glmnet object. 
  # if survival, add linear predictor `lp` in the output object
  
  n = nrow(X)
  p = ncol(X)
  b = sort(unique(blocks))
  
  if (is.null(bootID)) {
    bootID <- sapply(seq_len(nboot), function(x) sample(seq_len(n), n, replace = TRUE), simplify = FALSE)
  }
  
  brierDat <- list()
  coxDat   <- c()
  CVMods   <- list()
  
  lambdasMat <- c()
  maxLam  <- c()
  lambdas <- NULL
  lambdaSeqs <- NULL
    
  for(i in seq_along(bootID)) {  
    
    idSamples <- bootID[[i]]
    # select features by blocks
    idFeatures <- c()
    for(j in seq_along(b)) {
      idFeatures <- c(idFeatures, sample(which(blocks == b[j]), randomP[j]))
    }
    
    if (is.null(foldid)){
      fd = sample(rep(1:nfold, length = length(idSamples)), length(idSamples))
    } else {
      if (!is.list(foldid)) {
        fd = foldid
      } else{
        fd = foldid[[i]]
      }
    }
    
    if (is.function(method)) {
      # default arguments
      default.args <- list(X = data.matrix(X[idSamples, idFeatures]),
                           y = y[idSamples],
                           foldid = fd)
      # arguments in method function
      argsList <- formals(method)
      # arguments 
      argsList[names(default.args)] <- default.args
      argsList$blocks <- blocks[inFeatures]
      mods <- do.call(method, argsList)
      
    } else {
      
      if (method == "bypf") 
        mods <- pf.mods(X[idSamples, idFeatures], y[idSamples], 
                        family = "cox", 
                        foldid = fd, 
                        blocks = blocks[idFeatures], 
                        penalty.factors = penalty.factors,
                        lambdas = lambdas,
                        lambdaSeq = lambdaSeqs,
                        nlam = 50)
      
      if (method == "bycomp") 
        mods <- comp.mods(X[idSamples, idFeatures], y[idSamples], 
                          family = "cox", 
                          foldid = fd, 
                          blocks = blocks[idFeatures],
                          lambdas = lambdas,
                          lambdaSeqs = lambdaSeqs,
                          nlam = 50)
    }
    
    # get lambda
    if (i <= fastTune) {
      
      lambdasMat <- rbind(lambdasMat, unlist(lapply(mods, `[[`, "lambda")))
      
      #remove possible Inf and NA
      maxLam <- rbind(maxLam, unlist(lapply(mods, function(x) {
        tmp <- x[["lambdaSeq"]]
        max(tmp[is.finite(tmp)], na.rm = TRUE)
      })))
      
    } else if(i == (fastTune + 1)) {
      lambdas <- colMeans(lambdasMat)
      maxLam <- apply(maxLam, 2, max)
      
      lambdaSeqs <- mapply(function(maxL, minL) {
        exp(seq(log(maxL), log(minL), length.out = 50))
      }, maxL = maxLam, minL = lambdas, SIMPLIFY = FALSE)
    }
    
    # simplify the mods
    simMod <- lapply(mods, function(m){
      beta <- coef(m)[,1]
      beta <- beta[beta != 0]
      list(beta = beta, lp = as.vector(m$lp))
    })
    
    # generate y based on trainID and the saved y in the outcomes
    CVMods[[i]] <- list(mods = simMod, trainID = idSamples)
    
    #- should be more than one models
    if (length(mods) <= 1)
      stop("More than one models should be included")
    
    #- Get Surv Prob
    betas <- lapply(mods, function(x) {
      b <- coef(x)[,1]
      b[b != 0]
    })
    
    mats <- c()
    idSamplesRest <- setdiff(seq_len(n), idSamples)
    
    for(j in seq_along(betas)){
      
      beta  <- betas[[j]]
      probs <- survP(beta, X[idSamples, ], y[idSamples], X[idSamplesRest,], times)
      mat <- c()
      # Get data for Brier score
      # Two parameters for aggregating the formulas 
      samID <- c()
      timeID <- c()
      
      for(k in seq_along(times)) {
        
        bs <- brierScore(X[idSamples, blocks == block.clin],
                         X[idSamplesRest, blocks == block.clin],
                         y[idSamples,], y[idSamplesRest,], 
                         probs[ ,k], 
                         times[k],
                         cens.method = cens.method)
        
        mat <- rbind(mat, bs)
        samID <- c(samID, idSamplesRest[attr(bs, "sampleID")])
        timeID <- c(timeID, rep(times[k], nrow(bs)))
      }
      
      colnames(mat)[3] <- paste(j)
      if (j == 1) {
        mats <- mat
      } else {
        mats = cbind(mats, mat[, 3, drop = F])
      }
    }
    
    # brierDat will be created, it does not take much time anyway
    brierDat <- rbind(brierDat, cbind(mats, id = samID, time = timeID))
    
    # get survP at cox.time
    if (ifCox){
      
      if (is.null(cox.time))
        cox.time <- quantile(c(0, max(time[y[,2] == 1])), probs = 0.5)
      
      coxP <- c()
      for(j in seq_along(betas)){
        beta  <- betas[[j]]
        probs <- survP(beta, X[idSamples,], y[idSamples], X[idSamplesRest,], cox.time)
        coxP = cbind(coxP, probs)
      }
      coxDat <- rbind(coxDat, cbind(coxP, id = idSamplesRest))
    }
    #- get track information 
  }
  
  # bagging the prediction
  
  brierDat <- brierDat %>% 
    as.data.frame() %>%
    group_by(id, time) %>%
    summarise_at(vars(everything()), aggFun) %>%
    ungroup() %>%
    select(-c("id", "time"))
  
  # output brierDat, anyway
  out <- list(CVpred = as.matrix(brierDat),
              CVMods = CVMods,
              times = times,
              y = y,
              lambdas = lambdasMat,
              aggFun = aggFun,
              family = "cox")
  
  # output CVpredCox
  if(ifCox) {
    
    # bagging the prediction
    coxDat <- coxDat %>%
      as.data.frame() %>%
      group_by(id) %>%
      summarise_at(vars(everything()), aggFun) %>%
      ungroup() %>%
      arrange(id)
    
    out$CVpredCox <- as.matrix(coxDat[, !names(coxDat) %in% "id", drop = FALSE])
    out$cox.time  <- cox.time
    # y is only necessary in ifCox
    # update y, in case not all y selected
    out$y <- y[coxDat$id]
  }
  
  class(out) <- c("random.stk", "random.stk_cox", "stk", "stk_cox", "list") 
  # stk_cox used for merge
  return(out)
}

#-------------------------------------------------------------------------------
# super learner cox function
#-------------------------------------------------------------------------------

#' Stacked Cox Super Learner
#'
#' Applies a super learner approach to the outputs of the `StkCoxBL` function. 
#' It combines predictions from multiple models to optimize prediction accuracy.
#'
#' @param obj An object of class 'stk' typically output from `StkCoxBL`.
#' @param optimloss Optimization criteria, one of "ibsLoss", "logLik", or "PaLogLik".
#' @param optimFun Optimization function, either "lasso" or "ridge".
#' @param moreOptimFun Additional custom optimization functions, defaults to NULL. Using obj as the only input
#' @param intercept Logical indicating if an intercept should be included in the model.
#' @param foldid Vector indicating fold IDs for cross-validation, auto-generated if NULL.
#' @param foldid.cox Vector indicating fold IDs for Cox-specific cross-validation, auto-generated if NULL.
#' @param nfold Number of folds for cross-validation.
#' @param limits Numerical value setting limits for model coefficients.
#' @param ifWeights Logical indicating if weighting should be applied in model optimization.
#'
#' @return A list with optimized model weights and the input object enhanced with super learner results.
#' @examples
#' @export StkCoxSL

random.StkCoxSL <- function(obj, 
                            optimLoss = c("ibsLoss", "logLik", "PaLogLik"),
                            optimFun = c("lasso", "ridge"),
                            moreOptimFun = NULL, 
                            intercept = TRUE, 
                            foldid = NULL,
                            foldid.cox = NULL,
                            nfold = 10,
                            limits = 0,
                            ifWeights = TRUE,
                            truncValue = 0, 
                            transPFun = function(x) x,
                            ...) {
  
  # transform to the obj same format with StkCoxSL
  if ("CVpredCox" %in% names(obj)){
    obj$subMod$y <- obj$y
  }
  
  SL <- StkCoxSL(obj = obj, 
                 optimLoss = optimLoss,
                 optimFun = optimFun,
                 moreOptimFun = moreOptimFun,
                 intercept = intercept,
                 foldid = foldid,
                 foldid.cox = foldid.cox,
                 nfold = nfold,
                 limits = limits,
                 ifWeights = ifWeights,
                 truncValue = truncValue,
                 transPFun = transPFun)
  
  return(SL)
}

#-------------------------------------------------------------------------------

random.stack.cox <- function(X, y,
                             times = NULL,
                             bootID = NULL, 
                             randomP = c(ncol(X)), 
                             nboot = 100,
                             fastTune = nboot,
                             foldid = NULL,
                             nfold = 10,
                             blocks = rep(1, ncol(X)),
                             block.clin = 1,
                             ifCox = TRUE, 
                             cox.time = NULL,
                             method = c("bycomp", "bypf"),
                             cens.method = "equal",
                             penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                             aggFun = mean,
                             optimLoss = c("ibsLoss", "logLik", "PaLogLik"), 
                             optimFun = c("lasso", "ridge"),
                             moreOptimFun = NULL,
                             intercept = TRUE,
                             foldidSL = NULL,
                             foldidSL.cox = NULL,
                             nfoldSL = 10,
                             limits = 0,
                             ifWeights = TRUE,
                             truncValue = 0, 
                             transPFun = function(x) x,
                             ...){
  
  # stk base learner
  obj <- random.StkCoxBL(X, y,
                         times = times,
                         bootID = bootID, 
                         randomP = randomP, 
                         nboot = nboot,
                         fastTune = fastTune,
                         foldid = foldid,
                         nfold = nfold,
                         blocks = blocks,
                         block.clin = block.clin,
                         ifCox = ifCox, 
                         cox.time = cox.time,
                         method = method,
                         cens.method = cens.method,
                         penalty.factors = penalty.factors,
                         aggFun = aggFun,
                         ...)
  
  # stk super learner
  out <- random.StkCoxSL(obj = obj, 
                         optimLoss = optimLoss,
                         optimFun = optimFun,
                         moreOptimFun = moreOptimFun, 
                         intercept = intercept, 
                         foldid = foldidSL,
                         foldid.cox = foldidSL.cox,
                         nfold = nfoldSL,
                         limits = limits,
                         ifWeights = ifWeights,
                         truncValue = truncValue, 
                         transPFun = transPFun,
                         ...)
  
  return(out)
}

random.predBLSurvP <- function(obj, newdata, times = NULL) {
  obj$subMod$y <- obj$y
  predBLSurvP(obj, newdata, times = times, ifCVMods = TRUE)
}

random.predBLSurvPCox <- function(obj, newdata, times = NULL, ifCVMods = F) {
  obj$subMod$y <- obj$y
  obj$subMod$times <- obj$times
  predBLSurvPCox(obj, newdata = newdata, times = times, ifCVMods = TRUE)
}

random.predSLSurvP <- function(BLSurvP, alphas) {
  predSLSurvP(BLSurvP = BLSurvP, alphas = alphas)
}

random.predSLSurvPCox <- function(BLSurvPCox, alphas) {
  predSLSurvPCox(BLSurvPCox = BLSurvPCox, alphas = alphas) 
}

predict.random.stack.cox <- function(obj, newdata, times = NULL, ...) {
  
  if(is.null(times))
    times <- obj$times
  
  # get submodel prediction  
  alphas = obj$alpha
  
  # when length(alpha) == 1, alpha is not a lit
  # when one optim and one cen.method, then alpha is not a list
  if (!is.list(alphas))
    alphas <- list(alphas)
  
  # NOTE: alpha is a list in the new version
  # lp = obj$subMod$lp
  # y  = obj$subMod$y
  
  # define the out, should be a list
  alphaNames <- names(alphas)
  # alphas for predSLSurvP
  ind <- c(grep("^ibsLoss", alphaNames), grep("^logLik", alphaNames))
  indCox <- grep("PaLogLik", alphaNames)
  
  out <- vector("list", length(alphaNames))
  names(out) <- alphaNames
  # loop for the method 
  # for logLik and ibsLoss
  if (length(ind) > 0) {
    BLSurvP <- random.predBLSurvP(obj = obj, newdata = newdata, times = times)
    out[alphaNames[ind]] <- random.predSLSurvP(BLSurvP, alphas = alphas[ind])
  }
  
  if (length(indCox) > 0) {
    BLSurvPCox <- random.predBLSurvPCox(obj = obj, newdata = newdata, times = times)
    out[alphaNames[indCox]] <- random.predSLSurvPCox(BLSurvPCox, alphas[indCox])
  }
  
  notPredAlphas <- alphaNames[-c(ind, indCox)]
  if (length(notPredAlphas) > 0)
    warning(paste0("Prediction not support: ", 
                   paste0(notPredAlphas, collapse = ",")))
  
  if (length(out) == 1)
    out <- out[[1]]
  return(out)
}

################################################################################
# 2 glm section
################################################################################

random.StkBL <- function(X, y,
                         family,
                         bootID = NULL, 
                         randomP = table(blocks), 
                         nboot = 100,
                         fastTune = 100,
                         foldid = NULL,
                         nfold = 10,
                         method = c("bycomp", "bypf"),
                         blocks = rep(1, ncol(X)),
                         weights = rep(1, nrow(X)),
                         penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                         ...) {
  
  # length of Y
  require(dplyr)
  
  # argument checking
  if (fastTune > nboot) {
    fastTune = nboot
    warning(paste0("fastTune is reduced to ", nboot))
  }
  
  n = nrow(X)
  p = ncol(X)
  b = sort(unique(blocks))
  
  if (is.null(bootID)){
    sampleWeights <- ifelse(weights != 0, 1, 0)
    for(j in seq_len(nboot)) {
      # two groups
      if (family == "binomial") {
        boot <- c()
        for(k in unique(y)) {
          ind <- which(y == k)
          tmp <- sample(ind, sum(sampleWeights[ind]), replace = TRUE, prob = sampleWeights[ind] / sum(sampleWeights[ind]))
          boot <- c(boot, tmp)
        }
      } else {
        boot <- sample(seq_len(n), n, replace = TRUE, prob = sampleWeights / sum(sampleWeights))
      }
      bootID[[j]] <- boot
    }
  }
  
  # CV prediction
  CVMods <- list()
  pred   <- c() 
  
  lambdasMat <- c()
  maxLam  <- c()
  lambdas <- NULL
  lambdaSeqs <- NULL
  
  for(i in seq_along(bootID)) {  
    
    # print(i)
    idSamples <- bootID[[i]]
    # select features by blocks
    idFeatures <- c()
    for(j in seq_along(b)) {
      idFeatures <- c(idFeatures, sample(which(blocks == b[j]), randomP[j]))
    }
    
    # get foldid
    if (is.null(foldid)){
      fd = sample(rep(1:nfold, length = length(idSamples)), length(idSamples))
    } 
    
    # customized function as sub-models
    if (is.function(method)) {
      
      # default arguments
      default.args <- list(X = X[idSamples, idFeatures],
                           y = y[idSamples],
                           family = "binomial",
                           weights = weights[idSamples],
                           foldid = fd)
      # arguments in method function
      argsList <- formals(method)
      # arguments 
      argsList[names(default.args)] <- default.args
      mods <- do.call(method, argsList)
      
    } else {
      
      if (method == "bypf")
        mods <- pf.mods(X[idSamples, idFeatures], 
                        y[idSamples],
                        family = family,
                        foldid = fd,
                        blocks = blocks[idFeatures],
                        weights = weights[idSamples],
                        penalty.factors = penalty.factors,
                        lambdas = lambdas,
                        lambdaSeq = lambdaSeqs,
                        nlam = 50)
      
      if (method == "bycomp")
        mods <- comp.mods(X[idSamples, idFeatures], 
                          y[idSamples], 
                          family = family, 
                          foldid = fd, 
                          weights = weights[idSamples],
                          blocks = blocks[idFeatures],
                          lambdas = lambdas,
                          lambdaSeqs = lambdaSeqs,
                          nlam = 50)
    }
    
    # get lambda
    if (i <= fastTune) {
      lambdasMat <- rbind(lambdasMat, unlist(lapply(mods, `[[`, "lambda")))
      
      #remove possible Inf and NA
      maxLam <- rbind(maxLam, unlist(lapply(mods, function(x) {
        tmp <- x[["lambdaSeq"]]
        max(tmp[is.finite(tmp)], na.rm = TRUE)
        })))
    } else if(i == (fastTune + 1)) {
      lambdas <- colMeans(lambdasMat)
      maxLam <- apply(maxLam, 2, max)
      
      lambdaSeqs <- mapply(function(maxL, minL) {
        exp(seq(log(maxL), log(minL), length.out = 50))
      }, maxL = maxLam, minL = lambdas, SIMPLIFY = FALSE)
    }
    
    # simplify the mods
    simMod <- lapply(mods, function(m){
      beta <- coef(m)[,1]
      beta <- beta[beta != 0]
      beta
    })
    
    # generate y based on trainID and the saved y in the outcomes
    CVMods[[i]] <- list(mods = simMod, trainID = idSamples)
    idSamplesRest <- setdiff(seq_len(n), idSamples)
    tmp <- Reduce(cbind, predict.mods(mods, X[idSamplesRest,]))
    colnames(tmp) <- names(mods)
    pred <- rbind(pred, cbind(tmp, id = idSamplesRest))
  }
  
  # bagging the prediction
  pred <- pred %>% 
    as.data.frame() %>%
    group_by(id) %>%
    summarise_at(vars(everything()), mean) %>%
    arrange(-desc(id)) %>%
    ungroup() 
  
  id <- as.numeric(pull(pred, id))
  if (length(setdiff(seq_len(n), id)) > 0)
    cat(paste0("Cases were not included for super learner: ", paste0(setdiff(seq_len(n), id), collapse = ",")))
  pred <- pred %>% dplyr::select(-id)
  attr(pred, "id") <- id
  
  # colnames(pred) <- names(mods)
  out <- list(CVpred = pred,
              y = y,
              weights = weights[id],
              method = method,
              CVMods = CVMods,
              family = family)
  
  class(out) <- c("random.stk", "random.stk_glm", "stk", "stk_glm", "list")
  return(out)
}

random.StkSL <- function(obj, 
                         optimFun = c("lasso", "ridge", "solnp"),
                         moreOptimFun = NULL, 
                         foldid = NULL,
                         nfold = 10,
                         limits = 0,
                         ifWeights = TRUE,
                         ifshrinkage = TRUE,
                         truncValue = 0, 
                         transPFun = function(x) x,
                         randomSL = FALSE, 
                         bootID = NULL,
                         nboot = 100,
                         ifRandomSL = FALSE,
                         ...) {
  
  # random, whether random in super learner
  # bootID, random ID if random is true
  # nboot, number of nboot
  # whether randomization in SL
  
  # get family
  family  = obj$family
  weights = obj$weights
  y = obj$y
  n = length(y)
  
  if(ifRandomSL) {
    if (is.null(bootID)){
      sampleWeights <- ifelse(weights != 0, 1, 0)
      bootID = list()
      for(j in seq_len(nboot)) {
        
        if(family == "binomial") {
          boot <- c()
          # two groups
          for(k in unique(y)) {
            ind <- which(y == k)
            tmp <- sample(ind, sum(sampleWeights[ind]), replace = TRUE, prob = sampleWeights[ind] / sum(sampleWeights[ind]))
            boot <- c(boot, tmp)
          }
        } else {
          boot <- sample(seq_len(n), n, replace = TRUE, prob = sampleWeights / sum(sampleWeights))
        }
        bootID[[j]] <- boot
      }
    }
  } else{
    bootID <- list(seq_along(weights))
  }
  
  alphas <- list()
  
  for(i in seq_along(bootID)) {
    
    id <- bootID[[i]]
    objNew <- obj
    objNew$CVpred <- as.matrix(objNew$CVpred[id, ])
    objNew$y <- objNew$y[id]
    objNew$weights <- objNew$weights[id]
    
    alphas[[i]] <- StkSL(objNew,  
                         optimFun = optimFun,
                         moreOptimFun = moreOptimFun, 
                         foldid = foldid,
                         nfold = nfold,
                         limits = limits,
                         ifWeights = ifWeights,
                         ifshrinkage = ifshrinkage, 
                         truncValue = truncValue, 
                         transPFun = transPFun)$alphas
    
  }
  
  if (!ifRandomSL) {
    alphas <- alphas[[1]]
  }
  
  # add attr used in prediction
  attr(alphas, "ifRandomSL") <- ifRandomSL
  obj$alphas <- alphas
  return(obj)
}
  
#-------------------------------------------------------------------------------

random.stack <- function(X, y,
                         bootID = NULL, 
                         randomP = c(ncol(X)), 
                         nboot = 100,
                         fastTune = nboot,
                         foldid = NULL,
                         nfold = 10,
                         blocks = rep(1, ncol(X)),
                         weights = rep(1, nrow(X)),
                         method = c("bycomp", "bypf"),
                         penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                         optimFun = c("lasso", "ridge", "solnp"),
                         moreOptimFun = NULL,
                         foldidSL = NULL,
                         nfoldSL = 10,
                         limits = 0,
                         ifWeights = TRUE,
                         ifshrinkage = TRUE, 
                         truncValue = 0, 
                         transPFun = function(x) x,
                         randomSL = FALSE, 
                         ifRandomSL = FALSE,
                         ...) {
  
  # bootID, and nboot, used for both BL and SL, if `ifRandomSL = TRUE`
  
  obj <- random.StkBL(X, y,
                      bootID = bootID, 
                      randomP = randomP, 
                      nboot = nboot,
                      fastTune = fastTune,
                      foldid = foldid,
                      nfold = nfold,
                      method = method,
                      blocks = blocks,
                      weights = weights,
                      penalty.factors = penalty.factors,
                      ...) 
  
  out <- random.StkSL(obj, 
                      optimFun = optimFun,
                      moreOptimFun = moreOptimFun, 
                      foldid = foldidSL,
                      nfold = nfoldSL,
                      limits = limits,
                      ifWeights = ifWeights,
                      ifshrinkage = ifshrinkage,
                      truncValue = truncValue, 
                      transPFun = transPFun, 
                      randomSL = randomSL, 
                      bootID = bootID,
                      nboot = nboot,
                      ifRandomSL = ifRandomSL,
                      ...)
  return(out)
}

#-------------------------------------------------------------------------------

random.predBL <- function(obj, newdata, type.BL = NULL) {
  predBL(obj, newdata, type.BL = type.BL, ifCVMods = TRUE)
}

random.predSL <- function(BL, alphas, type = c("response", "link")) {
  
  if (!attr(alphas, "ifRandomSL")) {
    alphas <- list(alphas)
  }
  
  for(i in seq_along(alphas)) {
    
    pred0 <- predSL(BL = BL, alphas = alphas[[i]], type = type)
    
    if (i == 1)
      pred = pred0
    else
      pred <- mapply(`+`, pred0, pred, SIMPLIFY = FALSE)
    
  }
  # average
  pred <- lapply(pred, function(x) x / length(alphas))
  return(pred)
}

predict.random.stk <- function(obj, newdata, 
                               type = c("response", "link"), 
                               type.BL = NULL){
  
  BL <- random.predBL(obj = obj, newdata = newdata, type.BL = type.BL) 
  SL <- random.predSL(BL, alphas = obj$alphas, type = type)
  return(SL)
}

#------------------------------------------------------------------------------
# 
# random.StkBL <- function(X, y,
#                          bootID = NULL, 
#                          randomP = table(blocks), 
#                          nboot = 100,
#                          fastTune = 100,
#                          foldid = NULL,
#                          nfold = 10,
#                          method = c("bycomp", "bypf"),
#                          blocks = rep(1, ncol(X)),
#                          weights = rep(1, nrow(X)),
#                          penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
#                          ...) {
#   
#   # length of Y
#   
#   require(dplyr)
#   
#   # argument checking
#   if (fastTune > nboot) {
#     fastTune = nboot
#     warning(paste0("fastTune is reduced to ", nboot))
#   }
#   
#   n = nrow(X)
#   p = ncol(X)
#   b = sort(unique(blocks))
#   
#   if (is.null(bootID)){
#     sampleWeights <- ifelse(weights != 0, 1, 0)
#     for(j in seq_len(nboot)) {
#       boot <- c()
#       tmp <- sample(seq_len(n), n, replace = TRUE, prob = sampleWeights / sum(sampleWeights))
#       bootID[[j]] <- tmp
#     }
#   }
#   
#   # CV prediction
#   CVMods <- list()
#   pred   <- c() 
#   
#   lambdasMat <- c()
#   maxLam  <- c()
#   lambdas <- NULL
#   lambdaSeqs <- NULL
#   
#   for(i in seq_along(bootID)) {  
#     
#     # print(i)
#     idSamples <- bootID[[i]]
#     # select features by blocks
#     idFeatures <- c()
#     for(j in seq_along(b)) {
#       idFeatures <- c(idFeatures, sample(which(blocks == b[j]), randomP[j]))
#     }
#     
#     # get foldid
#     if (is.null(foldid)){
#       fd = sample(rep(1:nfold, length = length(idSamples)), length(idSamples))
#     } 
#     
#     # customized function as sub-models
#     if (is.function(method)) {
#       
#       # default arguments
#       default.args <- list(X = X[idSamples, idFeatures],
#                            y = y[idSamples],
#                            family = "binomial",
#                            weights = weights[idSamples],
#                            foldid = fd)
#       # arguments in method function
#       argsList <- formals(method)
#       # arguments 
#       argsList[names(default.args)] <- default.args
#       mods <- do.call(method, argsList)
#       
#     } else {
#       
#       if (method == "bypf")
#         mods <- pf.mods(X[idSamples, idFeatures], 
#                         y[idSamples],
#                         family = "binomial",
#                         foldid = fd,
#                         blocks = blocks[idFeatures],
#                         weights = weights[idSamples],
#                         penalty.factors = penalty.factors,
#                         lambdas = lambdas,
#                         lambdaSeq = lambdaSeqs,
#                         nlam = 50)
#       
#       if (method == "bycomp")
#         mods <- comp.mods(X[idSamples, idFeatures], 
#                           y[idSamples], 
#                           family = "binomial", 
#                           foldid = fd, 
#                           weights = weights[idSamples],
#                           blocks = blocks[idFeatures],
#                           lambdas = lambdas,
#                           lambdaSeqs = lambdaSeqs,
#                           nlam = 50)
#     }
#     
#     # get lambda
#     if (i <= fastTune) {
#       lambdasMat <- rbind(lambdasMat, unlist(lapply(mods, `[[`, "lambda")))
#       
#       #remove possible Inf and NA
#       maxLam <- rbind(maxLam, unlist(lapply(mods, function(x) {
#         tmp <- x[["lambdaSeq"]]
#         max(tmp[is.finite(tmp)], na.rm = TRUE)
#       })))
#     } else if(i == (fastTune + 1)) {
#       lambdas <- colMeans(lambdasMat)
#       maxLam <- apply(maxLam, 2, max)
#       
#       lambdaSeqs <- mapply(function(maxL, minL) {
#         exp(seq(log(maxL), log(minL), length.out = 50))
#       }, maxL = maxLam, minL = lambdas, SIMPLIFY = FALSE)
#     }
#     
#     # simplify the mods
#     simMod <- lapply(mods, function(m){
#       beta <- coef(m)[,1]
#       beta <- beta[beta != 0]
#       beta
#     })
#     
#     # generate y based on trainID and the saved y in the outcomes
#     CVMods[[i]] <- list(mods = simMod, trainID = idSamples)
#     idSamplesRest <- setdiff(seq_len(n), idSamples)
#     tmp <- Reduce(cbind, predict.mods(mods, X[idSamplesRest,]))
#     colnames(tmp) <- names(mods)
#     pred <- rbind(pred, cbind(tmp, id = idSamplesRest))
#   }
#   
#   # bagging the prediction
#   pred <- pred %>% 
#     as.data.frame() %>%
#     group_by(id) %>%
#     summarise_at(vars(everything()), mean) %>%
#     arrange(-desc(id)) %>%
#     ungroup() 
#   
#   id <- as.numeric(pull(pred, id))
#   if (length(setdiff(seq_len(n), id)) > 0)
#     cat(paste0("Cases were not included for super learner: ", paste0(setdiff(seq_len(n), id), collapse = ",")))
#   pred <- pred %>% select(-id)
#   attr(pred, "id") <- id
#   
#   # colnames(pred) <- names(mods)
#   out <- list(CVpred = pred,
#               y = y,
#               weights = weights[id],
#               method = method,
#               CVMods = CVMods)
#   
#   class(out) <- c("random.stk", "random.stk_glm", "stk", "stk_glm","list")
#   return(out)
# }
# 
# random.StkBinSL <- function(obj, 
#                             optimFun = c("lasso", "ridge", "BS"),
#                             moreOptimFun = NULL, 
#                             foldid = NULL,
#                             nfold = 10,
#                             limits = 0,
#                             ifWeights = TRUE,
#                             ifshrinkage = TRUE,
#                             truncValue = 0, 
#                             transPFun = function(x) x,
#                             randomSL = FALSE, 
#                             bootID = NULL,
#                             nboot = 100,
#                             ifRandomSL = FALSE,
#                             ...) {
#   
#   # random, whether random in super learner
#   # bootID, random ID if random is true
#   # nboot, number of nboot
#   # whether randomization in SL
#   
#   weights = obj$weights
#   y = obj$y
#   
#   n = nrow(X)
#   p = ncol(X)
#   b = sort(unique(blocks))
#   
#   if(ifRandomSL) {
#     if (is.null(bootID)){
#       sampleWeights <- ifelse(weights != 0, 1, 0)
#       bootID = list()
#       for(j in seq_len(nboot)) {
#         boot <- c()
#         # two groups
#         for(k in unique(y)) {
#           ind <- which(y == k)
#           tmp <- sample(ind, sum(sampleWeights[ind]), replace = TRUE, prob = sampleWeights[ind] / sum(sampleWeights[ind]))
#           boot <- c(boot, tmp)
#         }
#         bootID[[j]] <- boot
#       }
#     }
#   } else{
#     bootID <- list(seq_along(weights))
#   }
#   
#   alphas <- list()
#   
#   for(i in seq_along(bootID)) {
#     id <- bootID[[i]]
#     objNew <- obj
#     objNew$CVpred <- objNew$CVpred[id, ]
#     objNew$y <- objNew$y[id]
#     objNew$weights <- objNew$weights[id]
#     
#     alphas[[i]] <- StkBinSL(objNew, 
#                             optimFun = optimFun,
#                             moreOptimFun = moreOptimFun, 
#                             foldid = foldid,
#                             nfold = nfold,
#                             limits = limits,
#                             ifWeights = ifWeights,
#                             ifshrinkage = ifshrinkage, 
#                             truncValue = truncValue, 
#                             transPFun = transPFun)$alphas
#     
#   }
#   
#   if (!ifRandomSL) {
#     alphas <- alphas[[1]]
#   }
#   
#   # add attr used in prediction
#   attr(alphas, "ifRandomSL") <- ifRandomSL
#   obj$alphas <- alphas
#   return(obj)
# }
# 
# 
# 
# 
# 




