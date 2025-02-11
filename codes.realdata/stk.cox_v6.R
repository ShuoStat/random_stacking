################################################################################
# versoin 6
# customize sub-models
# minor update: brierScore. Add case numbers
################################################################################

################################################################################
# historical updates
# version 2
# last update
# data, 30/4/2024
# 1, add track to stk.cox
# 2, output Brier score and Cox data 

# version 3
# date, 4/5/2024
# 1, use aic and bic for submodels

# version 4
# date, 7/5/2024
# 1, wrap super learner and base learner, seperately.
# 2, vectorlize ifWeights to a vector for each optimLoss
# 3, predict.stack.cox, split into two steps 1, get survival P; 2, predict

# version 5
# date, 9/8/2024
# update prediction using resampled sub-models: cross-validation submodels
################################################################################

#  based on the Idea of IPs ----------------------------------------------------

#  import pf.mods, comp.mods, predict.mods in stk.glm
# 
# survP <- function(beta, trainX, trainY, validX, times){
#   
#   #- lp, linear predictors obtained from training
#   #- trainX, trainY, training data
#   #- validX, validation data
#   #- times,  survival probabalities for the times 
#   
#   s <- names(beta)
#   lpFitted <- trainX[,s] %*% beta
#   lpPred   <- validX[,s] %*% beta
#   
#   datFitted <- data.frame(x = lpFitted)
#   lpCox     <- coxph(trainY ~ x, dat = datFitted, 
#                      init = 1, iter.max = 0, x = T)
#   
#   survP <- riskRegression::predictCox(object = lpCox,
#                                       newdata = data.frame(x = lpPred),
#                                       times = times,
#                                       se = FALSE,
#                                       iid = FALSE,
#                                       keep.times = FALSE,
#                                       type = "survival")$survival
#   survP
# }

#- calculate baseline hazard 

cumhz <- function(lp, y, times){
  
  require(survival)
  if (!is.Surv(y))
    stop("y should be a Surv object")
  
  n <- length(y)
  ord <- order(y[,1], y[,2])
  
  #- time should be ordered
  #- warning: Strata will not be considered
  
  riskRegression::baseHaz_cpp(starttimes = rep(0, n),
                              stoptimes  = y[ord ,1],
                              status = y[ord ,2],
                              eXb = exp(lp)[ord],
                              strata = rep(0, n),
                              nPatients = n,
                              nStrata = 1,
                              emaxtimes = max(y[,1]),
                              predtimes = times,
                              cause = 1,
                              Efron = T)$cumhazard
  
}

#- survP, calculate survival probability of validation data

survP <- function(beta, trainX, trainY, validX, times) {
  
  #- lp, linear predictors obtained from training
  #- trainX, trainY, training data
  #- validX, validation data
  #- times,  survival probabalities for the times 
  
  s <- names(beta)
  
  if (length(s) == 0) {
    
    lpFitted <- rep(0, nrow(trainX))
    lpPred   <- rep(0, nrow(validX))
  } else {
    lpFitted <- trainX[,s, drop = F] %*% beta
    lpPred   <- validX[,s, drop = F] %*% beta
  }
  chz <- cumhz(lpFitted, trainY, times)
  mapply(function(x) {exp(-x * exp(lpPred))}, x = chz, SIMPLIFY = T)
  
}


brierScore <- function(Xclin.t, Xclin.v, y.t, y.v, 
                       survProb, time, cens.method = "KM") {
  
  #- survP
  #- Xclin, clinical variables, for calculation inverse censoring probability
  #- y.t, for the calculation of censoring P
  #- y.v, for test data
  #- survProb, survival probability for validation data
  #- time, the time of BS, a single time point
  #- cens.method, a method for censoring: "Breslow", "KM", "equal", "medianSurv", 
  #- "clinicalAdjusted"
  #- order survival time
  #- if y.t == y.v and survProb denotes probability of the training, the function 
  #- returns the BS for the training data, which can be used in the optimzation of BS
  
  require(survival)
  cens.method <- match.arg(cens.method, 
                           c("Breslow", "KM", "equal", "clinicalAdjusted", "medianSurv"))
  
  #- check
  if (length(y.v) != length(survProb))
    stop("y.v and survProb cannot match")
  
  # order the training 
  ord <- order(y.t[, 1])
  Xclin.t <- Xclin.t[ord, ]
  ys  <- y.t[ord, 1]
  ds  <- y.t[ord, 2]
  
  # order the validation
  ord <- order(y.v[,1])
  Xclin.v <- Xclin.v[ord, ]
  ys.v  <- y.v[ord ,1]
  ds.v  <- y.v[ord ,2]
  survProb <- survProb[ord]
  
  # calculate censoring probablites
  if (cens.method == "Breslow") {
    y.cens <- survival::Surv(ys, I(1 - ds))
    modG   <- cumhz(rep(0, length(y.cens)), y.cens, ys)
    censP  <- exp(-modG * exp(0)) # censP, censoring probabality
  }
  
  if (cens.method == "KM") {
    modG  <- survfit(Surv(ys, I(1 - ds)) ~ 1, type = "kaplan-meier")							
    # Adjusted Delta Function at time point t
    censP <- rep(modG$surv, modG$n.event + modG$n.censor)								
    # Compute censoring survival at each observed time point (allowing for possible ties)
  }
  
  if (cens.method == "equal"){
    censP <- rep(1, length(ys))
  }
  
  if (cens.method == "clinicalAdjusted") {
    
    surv.c <- Surv(ys, I(1 - ds))
    modG   <- coxph(surv.c ~ ., data = as.data.frame(Xclin.t))
    lp     <- modG$linear.predictors
    lp.v   <- predict(modG, newdata = as.data.frame(Xclin.v))
    
    # get cumulative censoring probability at time
    chz    <- cumhz(lp, surv.c, times = time)
    censP  <- exp(-chz * exp(lp.v)) # censP, censoring probabality
    
    # get cumulative censoring probabilty at event or censoring times
    chzEvent  <- cumhz(lp, surv.c, times = ys.v)
    censPEvent <- mapply(function(lp.v, chz) {
      exp(-chz * exp(lp.v))
    }, lp.v = lp.v, chz = chzEvent) # censP, censoring probabality
    
    # if events before time, use the event time
    Gc <- ifelse(ys.v < time, censPEvent, censP)
  }
  
  if (cens.method == "medianSurv") {
    
    # get median survival time
    medtime <- median(survfit(y.t ~ 1))
    if (time <= medtime) {
      Gc <- rep(1, length(ys.v))
    } else {
      surv.c <- Surv(ys, I(1 - ds))
      modG   <- coxph(surv.c ~ ., data = as.data.frame(Xclin.t))
      lp     <- modG$linear.predictors
      lp.v   <- predict(modG, newdata = as.data.frame(Xclin.v))
      
      # get cumulative censoring probability at time
      chz    <- cumhz(lp, surv.c, times = medtime)
      censP  <- exp(-chz * exp(lp.v)) # censP, censoring probabality
      Gc     <- censP
    }
  }
  
  Zi <- 1 * {ys.v > time}		
  # Response in the weighted LS problem
  Di <- 1 * {{ds.v == 1} | {ys.v >= time}} 
  
  if (cens.method %in% c("Breslow", "KM", "equal")) {
    
    #- get 
    if (identical(y.t, y.v)){
      censP.v = censP
    } else {
      censP.v <- c()
      for(i in seq_along(y.v)){
        # get censP at time ys.v
        censP.v <- c(censP.v, censP[ys > ys.v[i]][1])
      }
    }
    
    # If time > ys.v, censoring at time t, censP.v at ys.v 
    # If time <= ys.v, not censored, censp.v at time 
    Gc <- ifelse(ys.v <= time, censP.v, censP.v[ys.v > time][1])
  }
  
  Wi <- Di / Gc
  
  # Compute weights for the weighted LS problem
  booRemove <- ! Wi %in% c(0, NaN, Inf, NA)
  
  # Remove any observations with weight 0 or infinity...
  # YS <- ys.v[booRemove == TRUE]
  # YS[YS > time] <- time
  
  out <- data.frame(weight = Wi[booRemove],
                    Zind   = Zi[booRemove],
                    survProb = survProb[booRemove])
  attr(out, "sampleID") <- seq_len(nrow(Xclin.v))[ord][booRemove]
  
  return(out)
}


#     # fit cox model
#     mod  <- coxph(surv ~ ., data = as.data.frame(X))
#     lp   <- mod$linear.predictors
#     
#     # high, high risk groups. 1, censoring
#     s <- surv[,1] >= bestCut | surv[,2] != 1
#     # if (any(!s)) {
#     #   warning("Several early censored cases will have zero weights")
#     # }
#     # 
#     # high risk group
#     high <- surv[, 1] < bestCut
#     
#     # get times: high time + bestCut
#     # time should be in increasing order
#     times <- sort(c(surv[, 1], bestCut), decreasing = F)
#     
#     # get cumulative hazard
#     chz <- cumhz(lp, surv, times = times)
#     probs <- mapply(function(x) {exp(-x * exp(lp))}, x = chz, SIMPLIFY = T)
#     
#     # get probabilities
#     timeUse <- ifelse(high, surv[, 1], bestCut)
#     ind     <- match(timeUse, times)
#     
#     p <- c()
#     for(i in seq_along(timeUse)) {
#       p <- c(p, probs[i, ind[i]])
#     }
#     
#     weights <- ifelse(s, 1 / p, 0)
#     setNames(weights, rownames(X))
#   }
#   
#   
# }
#     


#- test for brierScore
# # 
# # set.seed(130971)
# library(prodlim)
# library(survival)
# dat <- SimSurv(100)
# 
# dat.t <- dat[1:70, ]
# dat.v <- dat[1:30, ]
# 
# fit <- coxph(Surv(time, event) ~ X1 + X2, data = dat.t)
# #
# beta <- coef(fit)
# trainX = as.matrix(dat.t[,names(beta), drop = F])
# trainY = Surv(dat.t$time, dat.t$event)
# validX = as.matrix(dat.v[,names(beta), drop = F])
# validY = Surv(dat.v[,"time"], dat.v[,"event"])
# times  = trainY[,1]
# re <- survP(beta = beta, trainX = trainX, trainY = trainY, validX = trainX, times)
# #
# # for training data
# p1 <- brierScore(trainY, trainY, re[,1], times[1], cens.method = "KM")
# p2 <- brierScore1(trainY,  re[,1], times[1], cens.method = "KM")
# cbind(p1, p2)
# 
# # for validation data
# re <- survP(beta = beta, trainX = trainX, trainY = trainY, validX = validX, times)
# p3 <- brierScore(trainY, validY, re[,1], times[1], cens.method = "KM")

#- oringal version
# brierScore1 <- function(Surv, survProb, time, cens.method = "Breslow"){
#   
#   #- survP
#   
#   #- cens.method, a method for censoring: "Breslow", "KM", "equal"
#   #- order survival time
#   ord <- order(Surv[,1])
#   survProb <- survProb[ord]
#   ys  <- Surv[ord, 1]
#   ds  <- Surv[ord, 2]
#   
#   #- 
#   Zi <- 1 * {ys > time}											
#   # Response in the weighted LS problem
#   Di <- 1 * {{ds == 1} | {ys >= time}} 
#   
#   if (cens.method == "Breslow") {
#     y.cens <- survival::Surv(ys, 1 - ds)
#     modG   <- cumhz(rep(0, length(y.cens)), y.cens, ys)
#     censP  <- exp(-modG * exp(0)) #- censP, censoring probabality
#   }
#   
#   if (cens.method == "KM") {
#     modG <- survfit(Surv(ys, I(1 - ds)) ~ 1, type = "kaplan-meier")							
#     # Adjusted Delta Function at time point t
#     censP <- rep(modG$surv, modG$n.event + modG$n.censor)								
#     # Compute censoring survival at each observed time point (allowing for possible ties)
#   }
#   
#   if (cens.method == "equal")
#     censP <- rep(1, length(ys))
#   
#   Gc <- ifelse(ys <= time, censP, censP[ys > time][1])
#   Wi <- Di / Gc		
#   
#   # Compute weights for the weighted LS problem
#   booRemove <- ! Wi %in% c(0, NaN, Inf)  
#   
#   # Remove any observations with weight 0 or infinity...
#   YS <- ys[booRemove == TRUE]
#   YS[YS > time] <- time
#   
#   survProb = survProb[booRemove == TRUE]	
#   
#   return(data.frame(weight = Wi[booRemove],
#                     Zind   = Zi[booRemove],
#                     survProb = survProb))
# }
# 
# 

# it should be optimIBS(), anyway, I do not want to rename the function
optimIPS <- function(X, y, weights = rep(1, nrow(X))){
  
  #- p, number of parameters
  optfun = function(x){ 
    sum(weights * (y - X %*% x)^2)
  }			
  
  sm1 <- function(x){																		
    fnl <- sum(x)
    fnl
  }
  
  #- m, number of models included
  m <- ncol(X)
  stk <- Rsolnp::solnp(rep(1/m, m), fun = optfun, eqfun = sm1, eqB = 1, 
                       LB = rep(0, m),
                       UB = rep(1, m),
                       control = list(trace = 0))$pars
  return(stk)
}


# it should be optimIBS(), anyway, I do not want to rename the function
# optimLogistic <- function(obj){
#   
#   obj <- as.matrix(obj)
#   #- p, number of parameters
#   # optfun = function(x){ 
#   #   sum(obj[,"weight"] * (obj[,"Zind"] - obj[, -c(1, 2)] %*% x) ^ 2)
#   # }			
#   
#   optfun <- function(x){
#     - sum(obj[,"weight"] * (obj[,"Zind"] * (obj[, -c(1, 2)] %*% x) -
#             log(1 + exp(obj[, -c(1, 2)] %*% x))))
#   }
#   
#   sm1 <- function(x){																		
#     fnl <- sum(x)
#     fnl
#   }
#   
#   #- m, number of models included
#   m <- ncol(obj) - 2
#   stk <- Rsolnp::solnp(rep(1/m, m), fun = optfun, eqfun = sm1, eqB = 1, 
#                        LB = rep(0, m),
#                        UB = rep(1, m),
#                        control = list(trace = 0))$pars
#   return(stk)
#   
# }

#' Base Learner for Stacked Cox 
#'
#' Stacking is a two-stage procedure. StkCoxBL fits the base learner of the survival stacking.
#' 
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
#' @param track Logical to track intermediate results.
#'
#' @return A list containing model fits, predictions, and possibly Cox predictions if `ifCox` is TRUE.
#' @examples
#'
#' @export StkCoxBL
#' 
# tmp

StkCoxBL <- function(X, y, 
                     times = NULL,
                     foldid = NULL,
                     foldid.internal = NULL,
                     blocks, 
                     block.clin = 1,
                     nfold = 10,
                     ifCox = TRUE, 
                     cox.time = NULL,
                     method = c("bycomp", "bypf"),
                     cens.method = "equal",
                     penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                     track = FALSE,
                     ...) {
  
  require(dplyr)
  require(survival)
  # cens.method, a method for IBS censoring: "Breslow", "KM", "equal", "clinicalAdjusted", 
  # "medianSurv"
  # intercept, if intercept is required for lasso and ridge
  # method, either c("bycomp", "bypf") or function. When method is a function, X and y should be the mandatory arguments, and the output should be a list of glmnet object. 
  # if survival, add linear predictor `lp` in the output object
  
  n = length(y)
  
  if (is.null(foldid)){
    foldid = sample(rep(1:nfold, length = n), n)
  }
  nfold = length(unique(foldid))
  
  # get foldid internal
  if (is.null(foldid.internal)){
    foldid.internal = list()
    for (i in 1:nfold){
      ns <- sum(foldid != i)
      foldid.internal[[i]] <-  sample(rep(1:nfold, length = ns), ns)
    }
  }
  
  #- Get times
  if (is.null(times)){
    time   <- y[ ,1]
    status <- y[ ,2]
    times <- quantile(c(0, max(time[status == 1])), 
                      probs = seq(0.1, 0.9, length = 10))
  }
  
  #- Get data for brier score
  #- brierDat for all folds and methods
  
  brierDat <- list()
  coxDat   <- c()
  CVMods   <- list()
  
  for(i in seq_len(nfold)) {  
    
    s <- foldid != i
    
    if (is.function(method)) {
      
      # default arguments
      default.args <- list(X = data.matrix(X[s,]),
                           y = y[s],
                           foldid = foldid.internal[[i]])
      # arguments in method function
      argsList <- formals(method)
      # arguments 
      argsList[names(default.args)] <- default.args
      mods <- do.call(method, argsList)
      
    } else {
      
      if (method == "bypf")
        mods <- pf.mods(data.matrix(X[s,]), y[s], 
                        family = "cox", 
                        foldid = foldid.internal[[i]], 
                        blocks = blocks, 
                        penalty.factors = penalty.factors)
      
      if (method == "bycomp")
        mods <- comp.mods(data.matrix(X[s,]), y[s], 
                          family = "cox", 
                          foldid = foldid.internal[[i]], 
                          blocks = blocks)
    }
    
    if (track) {
      # simplify the mods
      simMod <- lapply(mods, function(m){
        beta <- coef(m)[,1]
        beta <- beta[beta != 0]
        list(beta = beta, lp = as.vector(m$lp))
      })
      
      # generate y based on trainID and the saved y in the outcomes
      CVMods[[i]] <- list(mods = simMod, trainID = s)
    }
    
    #- should be more than one models
    if (length(mods) <= 1)
      stop("More than one models should be included")
    
    #- Get Surv Prob
    betas <- lapply(mods, function(x) {
      b <- coef(x)[,1]
      b[b != 0]
    })
    
    mats <- c()
    
    for(j in seq_along(betas)){
      
      b     <- betas[[j]]
      probs <- survP(b, X[s,], y[s], X[!s,], times)
      mat <- c()
      #- Get data for Brier score
      for(k in seq_along(times)) {
        mat <- rbind(mat, 
                     brierScore(X[s, blocks == block.clin],
                                X[!s, blocks == block.clin],
                                y[s,], y[!s,], 
                                probs[,k], 
                                times[k],
                                cens.method = cens.method))
      }
      
      if (j == 1)
        mats = mat
      else
        mats = cbind(mats, mat[ ,"survProb"])
    }
    
    # brierDat will be created, it does not take much time anyway
    brierDat <- rbind(brierDat, mats)
    
    # get survP at cox.time
    if (ifCox){
      
      if (is.null(cox.time))
        cox.time <- quantile(c(0, max(time[y[,2] == 1])), probs = 0.5)
      
      coxP <- c()
      for(j in seq_along(betas)){
        b     <- betas[[j]]
        probs <- survP(b, X[s,], y[s], X[!s,], cox.time)
        coxP = cbind(coxP, probs)
      }
      coxDat <- rbind(coxDat, coxP)
    }
    #- get track information 
  }
  
  # re-order coxDat according to y
  ord <- match(seq(n), order(foldid, decreasing = F))
  coxDat <- coxDat[ord, ]
  
  # refit the models
  if (is.function(method)) {
    
    argsList$X <- X
    argsList$y <- y
    argsList$foldid <- foldid
    mods <- do.call(method, argsList)
    
  } else {
    
    if (method == "bypf")
      mods <- pf.mods(X, y, 
                      family = "cox", 
                      foldid = foldid, 
                      blocks = blocks, 
                      penalty.factors = penalty.factors)
    
    if (method == "bycomp")
      mods <- comp.mods(X, y, 
                        family = "cox", 
                        foldid = foldid, 
                        blocks = blocks)
  }
  
  beta <- mapply(function(obj){
    b <- coef(obj)[ ,1]
    b[b != 0]
  }, obj = mods, SIMPLIFY = F)
  
  lp <- lapply(mods, `[[`, "lp")
  
  subMod <- list(beta = beta, 
                 lp = lp, 
                 y = y, 
                 times = times,
                 cens.method = cens.method)
  
  # CV prediction
  CVpred <- brierDat
  
  # output brierDat, anyway
  out <- list(subMod = subMod,
              CVpred = brierDat)
  
  # output CVpredCox
  if(ifCox) {
    out$CVpredCox <- coxDat
    out$cox.time  <- cox.time
  }
  
  if(track)
    out$CVMods <- CVMods
  
  class(out) <- c("stk", "stk_cox", "list") 
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


StkCoxSL <- function(obj, 
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
  
  # change to base environment, important to save memory size
  environment(transPFun) <- baseenv()
  
  # make ifWeights align with optimLoss
  if (length(ifWeights) == 1)
    ifWeights <- rep(ifWeights, length(optimLoss))
  
  ifWeights <- setNames(ifWeights, optimLoss)
  
  if(!inherits(obj, "stk"))
    stop("obj should be the output of StkCoxBL")
  
  ## argument checking
  # optimLoss <- match.arg(optimLoss, c("ibsLoss", "logLik", "PaLogLik"))
  # optimFun <- match.arg(optimFun, c("lasso", "ridge"))
  
  # argument checking
  if ("PaLogLik" %in% optimLoss) {
    if (!"CVpredCox" %in% names(obj))
      stop("CVpredCox not found")
  }
  
  #
  if (is.null(foldid))
    foldid <- sample(rep(1:nfold, length = nrow(obj$CVpred)), 
                     nrow(obj$CVpred))
  
  if (is.null(foldid.cox))
    foldid.cox <- sample(rep(1:nfold, length = nrow(obj$CVpredCox)), 
                         nrow(obj$CVpredCox))
  
  wts <- list()
  if ("ibsLoss" %in% optimLoss){
    
    if (ifWeights["ibsLoss"]) {
      weights = obj$CVpred[,1]
    } else {
      weights = rep(1, nrow(obj$CVpred))
    }
    
    X <- as.matrix(obj$CVpred[, -c(1, 2), drop = FALSE])
    y <- obj$CVpred[, 2]
    w <- optimIPS(X, y, weights)
    attr(w, "optimLoss") <- "ibsLoss"
    attr(w, "ifWeights") <- ifWeights["ibsLoss"]
    wts[["ibsLoss"]] <- w
  }
  
  ## super learner
  SL <- function(X, y, family, alpha, intercept, 
                 lower.limits, upper.limits,
                 foldid = NULL, nfold = 10, 
                 weights = NULL){
    
    if(is.null(weights)) 
      weights = rep(1, nrow(X))
    
    if (is.null(foldid))
      foldid = sample(rep(1:nfold, length = nrow(X)), nrow(X))
    
    cvMod <- cv.glmnet(X, y,
                       family = family,
                       weights = weights,
                       foldid = foldid,
                       lower.limits = lower.limits,
                       upper.limits = upper.limits,
                       intercept = intercept, 
                       alpha = alpha)
    
    mod <- glmnet(X, y,
                  family = family,
                  weights = weights,
                  intercept = intercept,
                  lambda = cvMod$lambda.min,
                  lower.limits = lower.limits,
                  upper.limits = upper.limits,
                  alpha = alpha)
    
    coef(mod)[ ,1]
  }
  
  ## weights
  optimLoss <- setdiff(optimLoss, "ibsLoss")
  
  if (length(optimLoss) > 0) {
    for(i in optimLoss) {
      for(j in optimFun) {
        
        alpha <- switch(j, "lasso" = 1, "ridge" = 0)
        
        if (i == "logLik") {
          
          family  <- "binomial"
          upper.limits = Inf
          lower.limits = limits
          f = foldid
          # truncate and transform
          X <- obj$CVpred[ , -c(1, 2), drop = FALSE] %>%
            as.data.frame() %>%
            mutate(across(everything(),
                          ~ case_when((.x) < truncValue ~ truncValue,
                                      (.x) > (1 - truncValue) ~ 1 - truncValue,
                                      .default = .x))) %>%
            mutate(across(everything(), transPFun)) %>%
            as.matrix()
          
          y = obj$CVpred[ ,2]
          
          # if weights
          if (ifWeights["logLik"])
            weights = obj$CVpred[ ,1]
          else
            weights = rep(1, nrow(obj$CVpred))
          
        } else if(i == "PaLogLik") {
          
          family  <- "cox"
          upper.limits = limits
          lower.limits = -Inf
          f = foldid.cox
          X <- obj$CVpredCox %>%
            as.data.frame() %>%
            mutate(across(everything(),
                          ~ case_when((.x) < truncValue ~ truncValue,
                                      (.x) > (1 - truncValue) ~ 1 - truncValue,
                                      .default = .x))) %>%
            mutate(across(everything(), transPFun)) %>%
            as.matrix()
          
          y = obj$subMod$y
          weights = rep(1, nrow(X))
          
          if (ifWeights["PaLogLik"])
            warning("Weights not possible for PaLogLik")
          
        }
        
        w <- SL(X, y, family = family, alpha = alpha, 
                intercept = intercept, 
                lower.limits = lower.limits,
                upper.limits = upper.limits,
                foldid = f,
                nfold = nfold,
                weights = weights)
        
        attr(w, "optimLoss") <- i
        attr(w, "optimFun")  <- j
        attr(w, "ifWeights") <- ifelse(i == "PaLogLik", FALSE, ifWeights[i])
        attr(w, "truncValue") <- truncValue
        attr(w, "transPFun") <- transPFun
        
        if (family == "cox") {
          lp <- X %*% w
          attr(w, "lp") <- lp
        }
        wts[[paste0(i,":", j)]] <- w
      }
    }
  }
  
  if (!is.null(moreOptimFun)){
    if (is.function(moreOptimFun)) {
      moreOptimFun <- list(moreOptimFun)
    }
    
    if (is.null(names(moreOptimFun)))
      names(moreOptimFun) <- paste0("addOptim", seq_along(moreOptimFun))
    
    # run moreOptimFun
    for(i in seq_along(moreOptimFun)){
      fun <- moreOptimFun[[i]]
      wts[[names(moreOptimFun)[i]]] <- do.call(fun, list(obj = obj))
    }
  }
  
  # reduce wts, if length is 1
  if (length(wts) == 1) 
    wts <- wts[[1]]
  
  # output
  obj$alpha <- wts
  return(obj)
}

#-------------------------------------------------------------------------------

stack.cox <- function(X, y,
                      times = NULL,
                      foldid = NULL,
                      foldid.internal = NULL,
                      nfold = 10,
                      blocks = rep(1, ncol(X)),
                      block.clin = 1,
                      ifCox = TRUE, 
                      cox.time = NULL,
                      method = c("bycomp", "bypf"),
                      cens.method = "equal",
                      penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                      track = FALSE, 
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
  obj <- StkCoxBL(X, y, 
                  times = times,
                  foldid = foldid,
                  foldid.internal = foldid.internal,
                  blocks = blocks, 
                  block.clin = block.clin,
                  nfold = nfold,
                  ifCox = ifCox, 
                  cox.time = cox.time,
                  method = method,
                  cens.method = cens.method,
                  penalty.factors = penalty.factors,
                  track = track)
  
  # stk super learner
  out <- StkCoxSL(obj, 
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
                  transPFun = transPFun)
  
  # if (!track)
  #   out[["CVMods"]] <- NULL
  
  return(out)
}

# default output the prediction probability
# predBLSurvP, prediction base learner of SurvP
# Internal function for predict.stack.cox

# internal function for predBLSurvP and predBLSurvPCox
# prediction for single Cox model
predBeta.cox <- function(beta, newdata, y, lp, times) {
  
  if (length(beta) == 0) {
    pred = rep(0, nrow(newdata)) 
  } else {
    pred <- newdata[, names(beta), drop = F] %*% beta
  }
  # cumulative hazard
  chz <- cumhz(lp, y, times)
  
  survP <- c()
  for(j in seq_along(times)){
    survP <- cbind(survP, exp(-chz[j] * exp(pred)))
  }
  
  colnames(survP) <- times
  return(survP)
}

# predict multiple models
predBetas.cox <- function(betas, newdata, y, lps, times) {
  
  mapply(predBeta.cox, beta = betas, lp = lps, 
         MoreArgs = list(newdata = newdata, times = times, y = y), SIMPLIFY = FALSE)
}

predBLSurvP <- function(obj, newdata, times = NULL, ifCVMods = FALSE) {
  
  # Stacking object from stack.cox or StkCoxSL
  # newdata, new data
  
  y = obj$subMod$y
  
  if (is.null(times))
    times = obj$subMod$times
  
  # if prediction based on CV sub-models
  if (ifCVMods) {
    
    cvMods <- obj$CVMods
    
    if (is.null(cvMods)) 
      stop("ifCVMods = TRUE: need CV models. Try set `track = TRUE` in stack.cox or StkCoxBL")
    
    # initial preds
    BLSurvP = 0
    for(i in seq_along(cvMods)){
      
      betas <- lapply(cvMods[[i]]$mods, `[[`, "beta")
      lps   <- lapply(cvMods[[i]]$mods, `[[`, "lp")
      y_sub <- y[cvMods[[i]]$trainID]
      
      pred <- predBetas.cox(betas = betas, newdata = newdata, y = y_sub, lps = lps, times = times)
      
      # add predictions, get average when divided by n
      BLSurvP <- mapply(`+`, x = pred, y = BLSurvP, SIMPLIFY = FALSE)
    }
    
    # get the averaged results
    BLSurvP <- lapply(BLSurvP, `/`, length(cvMods))
  }
  
  # if prediction based on universal model 
  if (!ifCVMods) {
    
    betas = obj$subMod$beta
    lps = obj$subMod$lp
    
    if (length(betas) != length(lps))
      stop("Number of beta and lp are different")
    
    # SurvP
    BLSurvP <- predBetas.cox(betas, newdata, y, lps, times)
  }

  class(BLSurvP) <- c(class(BLSurvP), "BLSurvP")
  attr(BLSurvP, "times") <- times
  
  return(BLSurvP)
}


predBLSurvPCox <- function(obj, newdata, times = NULL, ifCVMods = F) {
  
  # corresponding to lasso, adalasso
  cox.time <- obj$cox.time
  # lp <- obj$subMod$lp
  y  <- obj$subMod$y
  if (is.null(times))
    times = obj$subMod$times
  
  if (ifCVMods) {
    
    # get cvMods
    cvMods <- obj$CVMods
    
    if (is.null(cvMods)) 
      stop("ifCVMods = TRUE: need CV models. Try set `track = TRUE` in stack.cox or StkCoxBL")
    
    # initial preds
    BLSurvPCox = 0
    for(i in seq_along(cvMods)){
      
      betas <- lapply(cvMods[[i]]$mods, `[[`, "beta")
      lps   <- lapply(cvMods[[i]]$mods, `[[`, "lp")
      y_sub <- y[cvMods[[i]]$trainID]
      pred  <- predBetas.cox(betas = betas, newdata = newdata, y = y_sub, lps = lps, times = cox.time)
      
      # add predictions, get average when divided by n
      BLSurvPCox <- mapply(`+`, x = pred, y = BLSurvPCox, SIMPLIFY = FALSE)
    }
    
    # get the averaged results
    BLSurvPCox <- lapply(BLSurvPCox, `/`, length(cvMods))
  }
  
  # if prediction based on universal model 
  if (!ifCVMods) {
    
    betas <- obj$subMod$beta
    lps   <- obj$subMod$lp
    
    if (length(betas) != length(lps))
      stop("Number of beta and lp are different")
    
    # SurvP
    BLSurvPCox <- predBetas.cox(betas, newdata, y, lps, times = cox.time)
  }
  
  # list to data.frame
  BLSurvPCox <- as.data.frame(BLSurvPCox)
  colnames(BLSurvPCox) <- names(obj$subMod$beta)
  class(BLSurvPCox) <- c(class(BLSurvPCox), "BLSurvPCox")
  attr(BLSurvPCox, "cox.time") <- cox.time
  attr(BLSurvPCox, "times") <- times
  attr(BLSurvPCox, "y") <- y
  
  return(BLSurvPCox)
}


predSLSurvP <- function(BLSurvP, alphas) {
  
  if (!inherits(BLSurvP, "BLSurvP"))
    stop("BLSurvP should be object from predBLSurvP()")
  
  # alphas should be a list
  if (!is.list(alphas))
    alphas <- list(alphas)
  
  # define the out, should be a list
  out <- list()
  
  # get times
  times <- attr(BLSurvP, "times")
  
  # loop for the method 
  for(j in seq_along(alphas)){
    
    #- get predprob super model 
    w <- alphas[[j]]
    metaSurv <- matrix(NA, nrow = nrow(BLSurvP[[1]]), ncol = length(times),
                       dimnames = list(rownames(BLSurvP),
                                       as.character(times)))
    
    for(k in seq_along(times)){
      
      Z <- mapply(subset, x = BLSurvP,
                  MoreArgs = list(subset = T, 
                                  select = k), SIMPLIFY = T)
      
      if (attr(w, "optimLoss") == "logLik") {
        # truncate and transform
        truncValue <- attr(w, "truncValue")
        transPFun  <- attr(w, "transPFun")
        
        # set truncValue to be zero
        Z <- truncTrans(Z, truncValue, transPFun) 
        
        # if intercept
        if (length(w) == (ncol(Z) + 1)) {
          metaSurv[ ,k] <- 1 / (exp(-(cbind(1, Z) %*% w)) + 1)
        } else {
          metaSurv[ ,k] <- 1 / (exp(-(Z %*% w)) + 1)
        }
        
      } else {
        # if intercept
        if (length(w) == (ncol(Z) + 1)) {
          metaSurv[ ,k] <- cbind(1, Z) %*% w
        } else {
          metaSurv[ ,k] <- Z %*% w
        }
      }
    }
    out[[j]] <- metaSurv
  }
  
  if (length(out) == 1)
    out <- out[[1]]
  else
    names(out) <- names(alphas)
  
  return(out)
}

predSLSurvPCox <- function(BLSurvPCox, alphas) {
  
  if (!inherits(BLSurvPCox, "BLSurvPCox"))
    stop("BLSurvPCox should be object from predBLSurvPCox()")
  
  # alphas should be a list
  if (!is.list(alphas))
    alphas <- list(alphas)
  
  times <- attr(BLSurvPCox, "times")
  y <- attr(BLSurvPCox, "y")
  
  out <- list()
  # remove colnames, incase duplicated name
  colnames(BLSurvPCox) <- NULL
  
  
  for(j in seq_along(alphas)){
    
    alpha <- alphas[[j]]
    truncValue <- attr(alpha, "truncValue")
    transPFun  <- attr(alpha, "transPFun")
    coxLp      <- attr(alpha, "lp")
    
    ##!! need trunction in pred?
    BLSurvPCoxUse <- truncTrans(as.matrix(BLSurvPCox), truncValue, transPFun) 
    
    # chz in super learner
    chz <- cumhz(coxLp, y, times)
    # predict in super learner
    pred <- BLSurvPCoxUse %*% alpha
    
    # predict the survival P 
    predP <- c()
    for(h in chz) {
      predP <- cbind(predP, exp(-h * exp(pred)))
    }
    
    colnames(predP) <- times
    out[[j]] <- predP
  }
  
  names(out) <- names(alphas)
  return(out)
}


# BLSurvPCox <- predBLSurvPCox(obj, X.v)
# predict.stack.cox(mod$stk.pf, X.v)

predict.stack.cox <- function(obj, newdata, times = NULL, ifCVMods = FALSE) {
  
  if(is.null(times))
    times <- obj$subMod$times
  
  # get submodel prediction  
  alphas = obj$alpha
  
  # when length(alpha) == 1, alpha is not a lit
  # when one optim and one cen.method, then alpha is not a list
  if (!is.list(alphas))
    alphas <- list(alphas)
  
  # NOTE: alpha is a list in the new version
  lp = obj$subMod$lp
  y  = obj$subMod$y
  
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
    BLSurvP <- predBLSurvP(obj, newdata, times, ifCVMods = ifCVMods )
    out[alphaNames[ind]] <- predSLSurvP(BLSurvP, alphas[ind])
  }
  
  if (length(indCox) > 0) {
    BLSurvPCox <- predBLSurvPCox(obj, newdata, times = times, ifCVMods = ifCVMods)
    out[alphaNames[indCox]] <- predSLSurvPCox(BLSurvPCox, alphas[indCox])
  }
  
  notPredAlphas <- alphaNames[-c(ind, indCox)]
  if (length(notPredAlphas) > 0)
    warning(paste0("Prediction not support: ", 
                   paste0(notPredAlphas, collapse = ",")))
  
  if (length(out) == 1)
    out <- out[[1]]
  return(out)
}

#-------------------------------------------------------------------------------






































































