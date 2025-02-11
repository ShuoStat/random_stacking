################################################################################
## version 5
## date, 10/09/2024
## update pf.mods and comp.mods, using a sequence of lambdas instead of a single one
################################################################################


################################################################################
# History updates
## version 4
## date, 12/08/2024
## allow track and save CV models
## enable prediction based on CV models

# version 3
# date, 04/05/2024
# submodel based on aic and bic

# version 2
# date, 04/05/2024

# version 1
# no record left

#- based on the Idea of IPs ----------------------------------------------------

get.pcrs <- function(X,
                     ncomp = NULL,
                     blocks = rep(1, ncol(X)),
                     blocks.pc = c(1)){

  cts <- rep(0, length(blocks)) #- centers
  sds <- rep(1, length(blocks)) #- sds
  new.blocks <- c()

  M <- c() #- new data
  eigin.vectors <- list()

  for(i in unique(blocks)){

    s <- blocks == i
    Z <- X[, s]
    if (i %in% blocks.pc) {

      #- pcr, normalization and scaling
      Z <- scale(Z, scale = T, center = T)
      cts[s] <- attr(Z, "scaled:center")
      sds[s] <- attr(Z, "scaled:scale")
      SVD <- svd(Z)
      tmp <- Z %*% SVD$v[,1:ncomp]
      colnames(tmp) <- paste0("pc", i, "_", seq_len(ncol(tmp)))

      eigin.vectors[[i]] <- SVD$v[,1:ncomp]
      rownames(eigin.vectors[[i]]) <- colnames(Z)

    } else {
      tmp <- Z
    }
    M <- cbind(M, tmp)
    new.blocks <- c(new.blocks, rep(i, ncol(tmp)))
  }

  list(centers = cts,
       scales  = sds,
       blocks  = new.blocks,
       M = M,
       eigin.vectors = eigin.vectors)
}

#- transform survival Y to binary Y.
#- 5 year survival

biY <- function(y, cutTimes = NULL) {

  # times: choice for time cutoffs (years), if NULL, all times will be tried.
  if (is.null(cutTimes))
    cutTimes <- unique(sort(y[ ,1], decreasing = F))
  
  dif = 99999
  bestCut <- 99999
  
  for(i in cutTimes){
    
    time <- y[ ,1]
    status <- y[ ,2]

    status <- ifelse(time > i, 1, status)
    time   <- time[status != 0]

    #- rm zero
    tmp   <- abs(sum(time > i) - sum(time <= i))
    dif   <- ifelse(dif < tmp, dif, tmp)
    bestCut <- ifelse(dif < tmp, bestCut, i)
  }

  s <- y[,1] >= bestCut | y[,2] != 0
  # s <- y[,1] < bestCut & y[,2] != 0
  
  list(biY = y[s, 1] < bestCut,
       sel = s,
       bestCut = bestCut)
}


ic_glmnet <- function(fit, ic = c("AIC", "BIC", "AICc")) {
  
  ic = match.arg(ic, c("AIC", "BIC", "AICc"))
  
  # fit is the glmnet model
  lambda <- fit$lambda
  
  tLL <- fit$nulldev *  fit$dev.ratio
  k <- fit$df
  n <- fit$nobs
  
  out <- switch(ic, 
                "AIC" = - tLL + 2 * k,
                "AICc" = - tLL + 2 * k + 2 * k * (k + 1) / (n - k - 1),
                "BIC" = log(n) * k - tLL)
  
  return(out)
}


# function to extract glmnet object at a single lambda
extractGlmnet <- function(fit, lambda, tol = 1e-5) {
  
  # allow extract multiple lambdas
  
  lamSeq <- fit$lambda
  ind <- c()
  notFound <- c()
  
  for(i in lambda) {
    
    whichLam <- which(abs(i - lamSeq) < tol)[1]
    if (is.na(whichLam)){ 
      notFound <- c(notFound, i)
    } else {
      ind <- c(ind, whichLam) # use the first one
    }
  }
  
  
  if (length(ind) == 0){
    stop("Lambda not find in glmnet object")
  } else if(length(notFound) > 0) {
    warning("The folowing lambda is not found:", paste0(notFound, collapse = ","))
  }
  
  fit$beta <- fit$beta[,ind, drop = FALSE]
  fit$df   <- fit$df[ind]
  fit$lambda <- fit$lambda[ind]
  fit$dim[2] <- length(lambda) 
  fit$dev.ratio <-  fit$dev.ratio[ind]
  
  if (!"coxnet" %in% class(fit))
    fit$a0 <- fit$a0[ind]
  
  return(fit)
}

# update here in v5
# Do not supply a single value for lambda (for predictions after CV 
# use predict() instead). Supply instead a decreasing sequence of lambda 
# values. glmnet relies on its warms starts for speed, and its often faster to 
# fit a whole path than compute a single fit.

pf.mods <- function(X, y, 
                    family, 
                    foldid, 
                    blocks = rep(1, ncol(X)),
                    weights = rep(1, nrow(X)),
                    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                    lambdas = NULL,
                    lambdaSeq = NULL,
                    ...){
  
  if (is.null(lambdas)) {
    
    if (!is.null(lambdaSeq))
      lambdaSeq <- sort(lambdaSeq, decreasing = TRUE)
    
      cvs <- cv.ipflas(X, y, 
                       family = family,
                       foldid = foldid,
                       blocks = blocks,
                       weights = weights,
                       penalty.factors = penalty.factors,
                       lambda = lambdaSeq,
                       ...)
    
    lambdas = c()
    
    for(i in penalty.factors) {
      pf  <- paste0(i, collapse = "-")
      lambdas <- c(lambdas, cvs$statlist[[pf]]$lambda[which.min(cvs$statlist[[pf]]$cvm)])
      lambdaSeq[[pf]] <- cvs$statlist[[pf]]$lambda
      # lambdaSeq should be null or list for each penalty.factors
    }
  }
  
  #- get lambda
  mods <- list()
  
  if (is.null(lambdaSeq))
    lambdaSeq <- sort(unlist(lambdas), decreasing = TRUE)
  # if lambdaSeq is null, use lambdas as a short sequence
  
  for(i in seq_along(penalty.factors)){
    
    tmp <- ipflas(X, y,
                  family = family,
                  blocks = blocks,
                  weights = weights,
                  pf = penalty.factors[[i]],
                  lambda = lambdaSeq[[i]],
                  ...)
    
    tmp <- extractGlmnet(tmp, lambda = lambdas[i])
    tmp$lambdaSeq <- lambdaSeq[[i]]
    
    #- output lp for cox, for calculation of cumulative hazard
    if (family == "cox")
      tmp$lp <- predict(tmp, X)
    mods[[paste0(penalty.factors[[i]], collapse = "-")]] <- tmp
  }
  
  class(mods) <- "pf.mods"
  mods
}


comp.mods <- function(X, y,
                      family,
                      foldid,
                      weights = rep(1, nrow(X)),
                      blocks = rep(1, ncol(X)), 
                      lambdas = NULL,
                      lambdaSeqs = NULL,
                      ...){
  
  b <- sort(unique(blocks))
  
  # check lambdaSeqs
  if (is.null(lambdaSeqs))
    lambdaSeqs <- rep(list(NULL), length(b))
  
  if (length(lambdaSeqs) != length(b))
    stop("lambdaSeqs should be the same length with the blocks")
  
  mods <- list()
  
  # get lambdas
  if (is.null(lambdas)) {
    
    lambdas <- c()
    for(i in seq_along(b)) {
      s <- blocks == b[i]
      cvs <- cv.glmnet(X[ ,s], y, family = family, foldid = foldid, 
                       weights = weights, lambda = lambdaSeqs[[i]], ...)
      lambdaSeqs[[i]] <- cvs$lambda
      lambdas <- c(lambdas, cvs$lambda.min)
    }
  }
  
  # fit models
  for(i in seq_along(b)){
    
    s <- blocks == b[i]
    tmp <- glmnet(X[ ,s], y, family = family, weights = weights, lambda = lambdaSeqs[[i]], ...)
    
    # extract lambda
    tmp <- extractGlmnet(tmp, lambda = lambdas[i])
    tmp$block <- i
    tmp$lambdaSeq <- lambdaSeqs[[i]]
    
    if (family == "cox")
      tmp$lp <- predict(tmp, X[,s])
    
    mods[[paste0(i)]] <- tmp
  }

  class(mods) <- "comp.mods"
  mods
}

# tmp
# fit relaxed lasso sub-models
relax.mods <- function(X, y, 
                       family, 
                       pf = c(1, 1), 
                       blocks = rep(1, ncol(X)), 
                       gamma = c(0, 0.25, 0.5, 0.75, 1),
                       foldid = NULL,
                       ...) {
  
  cvs <- cv.ipfrelax(X, y, family = family, 
                     foldid = foldid, 
                     blocks = blocks, 
                     pfs = list(pf), 
                     gamma = gamma,
                     ...)
  
  # find lambda.min for each gamma
  
  mins <- lapply(cvs$statlist[[1]], function(obj) {
    obj$lambda[which.min(obj$cvm)]
  })
  
  mods <- list()
  
  for(i in seq_along(gamma)) {
    
    # set name
    nam <- paste0(paste0(pf, collapse = "-"), ":", gamma[i])
    tmpMod <- ipfrelax(X, y, family = family,
                       blocks = blocks,
                       pfs = list(pf),
                       lambda = mins[[i]],
                       gamma = i,
                       ...)
    
    # if cox, add lp
    if (family == "cox")
      tmpMod$lp <- predict(tmpMod, X)
    
    mods[[nam]] <- tmpMod
  }
  
  names(mods) <- gamma
  return(mods)
}


#- prediction --

predict.mods <- function(obj, newdata, type = "response"){

  if (length(obj) > 1){
    
    pred <- lapply(obj, function(x) {
      if ("coxnet" %in% class(obj))
        s <- rownames(coef(x))
      else
        s <- rownames(coef(x))[-1]
      predict(x, newdata[,s], type = type)})
  } else{
    pred <- predict(obj, newdata, type = type)
  }
  pred
}


predict.mods <- function(obj, newdata, type = "response"){
  
  if (length(obj) > 1){
    pred <- lapply(obj, function(x) {
      if ("coxnet" %in% class(x))
        s <- rownames(coef(x))
      else
        s <- rownames(coef(x))[-1]
      predict(x, newdata[ ,s, drop = FALSE], type = type)})
    
  } else{
    
    if (!"glmnet" %in% class(obj)) {
      obj <- obj[[1]]
    }
    
    if ("coxnet" %in% class(obj))
      s <- rownames(coef(obj))
    else
      s <- rownames(coef(obj))[-1]
    pred <- predict(obj, newdata[,s, drop = FALSE], type = type)
    pred <- list(pred)
  }
  pred
}

optimBi <- function(Z, y, weights){
  
  #- Z, predictions of sum-models
  Z <- as.matrix(Z)
  
  optfun <- function(x){
    sum(weights * (y - Z %*% x)^2)
  }
  
  m <- ncol(Z)
  sm1 <- function(x) sum(x)
  
  stk <- Rsolnp::solnp(rep(1/m, m), fun = optfun, eqfun = sm1, eqB = 1, 
                       LB = rep(0, m),
                       UB = rep(1, m),
                       control = list(trace = 0))$pars
  return(stk)  
}

#- stacking for binomial --

truncTrans <- function(X, truncValue, transPFun) {
  
  if (!is.function(transPFun))
    stop ("transPFun should be a function")
  
  X %>% 
    as.data.frame() %>%
    mutate(across(everything(),
                  ~ case_when((.x) < truncValue ~ truncValue,
                              (.x) > (1 - truncValue) ~ 1 - truncValue,
                              .default = .x))) %>%
    mutate(across(everything(), transPFun)) %>%
    as.matrix()
}

#-

StkBinBL <- function(X, y,
                     foldid = NULL,
                     foldid.internal = NULL,
                     seed = NULL,
                     nfold = 10,
                     method = c("bycomp", "bypf"),
                     blocks = rep(1, ncol(X)),
                     weights = rep(1, nrow(X)),
                     penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                     track = FALSE, 
                     ...) {
  
  # length of Y
  n = length(y)
  
  # seeds
  if (!is.null(seed))
    set.seed(seed)
  
  # foldid
  if (is.null(foldid)){
    foldid = sample(rep(1:nfold, length = n), n)
  }
  
  # nfold
  nfold = length(unique(foldid))
  
  # internal foldid
  if (is.null(foldid.internal)) {
    foldid.internal = list()
    for (i in 1:nfold){
      ns <- sum(foldid != i)
      foldid.internal[[i]] <-  sample(rep(1:10, length = ns), ns)
    }
  }
  
  # method
  # if (method == "bypf")
  #   pred <- matrix(NA, nrow = n, ncol = length(penalty.factors))
  # 
  # if (method == "bycomp")
  #   pred <- matrix(NA, nrow = n, ncol = length(unique(blocks)))
  
  # CV prediction
  CVMods <- list()
  pred   <- c() 
  
  for(i in seq_len(nfold)) {
    
    s <- foldid != i
    
    # customized function as sub-models
    if (is.function(method)) {
      
      # default arguments
      default.args <- list(X = X[s,],
                           y = y[s],
                           family = "binomial",
                           foldid = foldid.internal[[i]])
      # arguments in method function
      argsList <- formals(method)
      # arguments 
      argsList[names(default.args)] <- default.args
      mods <- do.call(method, argsList)
      
    } else {
      
      if (method == "bypf")
        mods <- pf.mods(X[s, ], y[s],
                        family = "binomial",
                        foldid = foldid.internal[[i]],
                        blocks = blocks,
                        weights = weights[s],
                        penalty.factors = penalty.factors, ...)
      
      if (method == "bycomp")
        mods <- comp.mods(X[s, ], y[s],
                          family = "binomial",
                          foldid = foldid.internal[[i]],
                          blocks = blocks,
                          weights = weights[s], ...)
    }
    
    
    if (track) {
      # simplify the mods
      simMod <- lapply(mods, function(m){
        beta <- coef(m)[,1]
        beta <- beta[beta != 0]
        beta
      })
      
      # generate y based on trainID and the saved y in the outcomes
      CVMods[[i]] <- list(mods = simMod, trainID = s)
    }
    
    tmp  <- Reduce(cbind, predict.mods(mods, X[!s,]))
    pred <- rbind(pred, tmp)
    # pred[!s,] <- tmp
  }
  
  # re-order the pred matrix 
  pred <- pred[rank(foldid, ties.method = "first"), ]
  colnames(pred) <- names(mods)
  
  # refit the models
  if (is.function(method)) {
    
    mods <- do.call(method, argsList)
    
  } else {
    
    if (method == "bypf")
      mods <- pf.mods(X, y,
                      family = "binomial",
                      foldid = foldid,
                      blocks = blocks,
                      weights = weights,
                      penalty.factors = penalty.factors, ...)
    
    if (method == "bycomp")
      mods <- comp.mods(X, y,
                        family = "binomial",
                        foldid = foldid,
                        weights = weights,
                        blocks = blocks, ...)
  }
  
  beta = mapply(function(obj){
    b <- coef(obj)[,1]
    b[b != 0]
  }, obj = mods, SIMPLIFY = F)
  
  # return object
  out <- list(CVpred = pred,
              y = y,
              beta = beta,
              weights = weights,
              method = method)
  
  if(track)
    out$CVMods <- CVMods
  
  class(out) <- c("stk", "stk_glm", "list")
  return(out)
}


StkBinSL <- function(obj, 
                     optimFun = c("lasso", "ridge", "BS"),
                     moreOptimFun = NULL, 
                     foldid = NULL,
                     nfold = 10,
                     limits = 0,
                     ifWeights = TRUE,
                     ifshrinkage = TRUE, 
                     truncValue = 0, 
                     transPFun = function(x) x,
                     ...){
  
  y = obj$y
  # change to base environment, important to save memory size
  environment(transPFun) <- baseenv()
  
  # make ifWeights align with optimLoss
  if (length(ifWeights) == 1)
    ifWeights <- rep(ifWeights, length(optimFun))
  
  ifWeights <- setNames(ifWeights, optimFun)
  
  if(!inherits(obj, "stk"))
    stop("obj should be the output of StkBinBL")
  
  if (is.null(foldid))
    foldid <- sample(rep(1:nfold, length = nrow(obj$CVpred)), 
                     nrow(obj$CVpred))
  
  wts <- list()
  
  for(i in optimFun) {
    
    if (i %in% c("lasso", "ridge")) {
      
      # truncate and transform
      pred <- truncTrans(obj$CVpred, truncValue, transPFun)

      # weights
      if (!ifWeights[i])
        weights = rep(1, nrow(pred))
      else
        weights = obj$weights
      
      alpha <- ifelse(i == "lasso", 1, 0)
      if (ifshrinkage) {
        
        cv.stack <- cv.glmnet(pred, 
                              y, 
                              family = "binomial", 
                              lower.limits = limits,
                              foldid = foldid,
                              weights = weights, 
                              alpha = alpha,
                              ...)
        lambda = cv.stack$lambda.min
        
      } else {
        lambda = 0
      }
      
      stack <- glmnet(pred, 
                      y, 
                      family = "binomial", 
                      lower.limits = limits, 
                      lambda = lambda,
                      weights = weights,
                      alpha = alpha,
                      ...)
      w <- coef(stack)[ ,1]
    }
    
    if (i == "BS") {
      pred <- obj$CVpred 
      w <- optimBi(pred, y, weights = weights)
    }
    
    attr(w, "optimFun")  <- i
    attr(w, "ifWeights") <- ifWeights[i]
    
    # no truncation and transformation using BS 
    if (i == "BS") {
      
      transPfunBS <- function(x) x
      environment(transPfunBS) <- baseenv()
      
      attr(w, "truncValue") <- 0
      attr(w, "transPFun")  <- transPfunBS
    } else {
      attr(w, "truncValue") <- truncValue
      attr(w, "transPFun")  <- transPFun
    }
    
    wts[[i]] <- w
  } 
  
  # moreOptimFun
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
  
  # return obj
  obj$alphas <- wts
  
  class(obj) <- c("StkBinSL", class(obj))
  return(obj)
  
}

stk.bin <- function(X, y,
                    foldid = NULL,
                    foldid.internal = NULL,
                    seed = NULL,
                    nfold = 10,
                    method = c("bycomp", "bypf"),
                    blocks = rep(1, ncol(X)),
                    weights = rep(1, nrow(X)),
                    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                    track = FALSE, 
                    optimFun = c("lasso", "ridge", "BS"),
                    moreOptimFun = NULL, 
                    foldidSL = NULL,
                    nfoldSL = 10,
                    limits = 0,
                    ifWeights = TRUE,
                    ifshrinkage = TRUE, 
                    truncValue = 0, 
                    transPFun = function(x) x,
                    ...) {
  
  # Base learner
  obj <- StkBinBL(X, y,
                  foldid = foldid,
                  foldid.internal = foldid.internal,
                  seed = seed,
                  nfold = nfold,
                  method = method,
                  blocks = blocks,
                  weights = weights,
                  penalty.factors = penalty.factors,
                  track = track,
                  ...)
  
  # Super learner
  out <- StkBinSL(obj, 
                  optimFun = optimFun,
                  moreOptimFun = moreOptimFun, 
                  foldid = foldidSL,
                  nfold = nfoldSL,
                  limits = limits,
                  ifWeights = ifWeights,
                  ifshrinkage = ifshrinkage, 
                  truncValue = truncValue, 
                  transPFun = transPFun,
                  ...)
  
  class(out) <- c("stk.bin", class(out))
  return(out)
  
}

# For merging sub-models

# 
# moreOptimFun <- function(obj) {
#   
#   stack <- glmnet(obj$CVpred, y,
#                   family = "binomial",
#                   weights = obj$weights,
#                   lambda = 0,
#                   lower.limits = 0)
#   alpha <- coef(stack)
#   alpha
#   
# }
#   
#   stack.binomial <- function(X, y,
#                              foldid = NULL,
#                              foldid.internal = NULL,
#                              seed = NULL,
#                              nfold = 10,
#                            method = c("bycomp", "bypf"),
#                            blocks = rep(1, ncol(X)),
#                            weights = rep(1, nrow(X)),
#                            optim = "glmnet",
#                            penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8))){
# 
#   n = length(y)
# 
#   if (!is.null(seed))
#     set.seed(seed)
# 
#   if (is.null(foldid)){
#     foldid = sample(rep(1:nfold, length = n), n)
#   }
# 
#   nfold = length(unique(foldid))
# 
#   if (is.null(foldid.internal)) {
#     foldid.internal = list()
#     for (i in 1:nfold){
#       ns <- sum(foldid != i)
#       foldid.internal[[i]] <-  sample(rep(1:10, length = ns), ns)
#     }
#   }
# 
#   if (method == "bypf")
#     pred <- matrix(NA, nrow = n, ncol = length(penalty.factors))
# 
#   if (method == "bycomp")
#     pred <- matrix(NA, nrow = n, ncol = length(unique(blocks)))
# 
#   for(i in seq_len(nfold)) {
# 
#     s <- foldid != i
#     if (method == "bypf")
#       mods <- pf.mods(X[s,], y[s],
#                       family = "binomial",
#                       foldid = foldid.internal[[i]],
#                       blocks = blocks,
#                       weights = weights[s],
#                       penalty.factors = penalty.factors)
# 
#     if (method == "bycomp")
#       mods <- comp.mods(X[s,], y[s],
#                         family = "binomial",
#                         foldid = foldid.internal[[i]],
#                         blocks = blocks,
#                         weights = weights[s])
# 
#     tmp <- Reduce(cbind, predict.mods(mods, X[!s,]))
#     pred[!s,] <- tmp
#   }
#   
#   if (optim == "glmnet"){
#     stack <- glmnet(pred, y, 
#                     family = "binomial", 
#                     weights = weights,
#                     lambda = 0, 
#                     lower.limits = 0)
#     alpha <- coef(stack)
#   }
#   
#   if (optim == "L1"){
#     cv.stack <- cv.glmnet(pred, y, 
#                           family = "binomial", 
#                           lower.limits = 0,
#                           foldid = foldid,
#                           weights = weights)
#     stack    <- glmnet(pred, y, 
#                        family = "binomial", 
#                        lower.limits = 0, 
#                        lambda = cv.stack$lambda.min,
#                        weights = weights)
#     alpha <- coef(stack)
#   }
#   
#   if (optim == "solnp"){
#     alpha <- optimBi(pred, y)
#   }
#   
#   #- refit the models
#   if (method == "bypf")
#     mods <- pf.mods(X, y,
#                     family = "binomial",
#                     foldid = foldid,
#                     blocks = blocks,
#                     weights = weights,
#                     penalty.factors = penalty.factors)
# 
#   if (method == "bycomp")
#     mods <- comp.mods(X, y,
#                       family = "binomial",
#                       foldid = foldid,
#                       weights = weights,
#                       blocks = blocks)
# 
#   beta = mapply(function(obj){
# 
#     b <- coef(obj)[,1]
#     b[b != 0]
#   }, obj = mods, SIMPLIFY = F)
# 
#   return(list(beta = beta, alpha = alpha, optim = optim))
# }
# 

# Post-CV method----------------------------------------------------------------

# stack.binomial.postCV <- function(X, y,
#                                   foldid = NULL,
#                                   seed = NULL,
#                                   nfold = 10,
#                                   method = c("bycomp", "bypf"),
#                                   blocks = rep(1, ncol(X)),
#                                   weights = rep(1, nrow(X)),
#                                   penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8))){
#   
#   n = length(y)
#   
#   if (!is.null(seed))
#     set.seed(seed)
#   
#   if (is.null(foldid)){
#     foldid = sample(rep(1:nfold, length = n), n)
#   }
#   
#   nfold = length(unique(foldid))
#   pred <- c()
#   mods <- list()
#   
#   if (method == "bycomp") {
#     
#     nblocks <- length(unique(blocks))
#     for(i in seq_len(nblocks)) {
#       
#       s   <- blocks == i
#       cvs <- cv.glmnet(X[ ,s], y, 
#                        family = "binomial", 
#                        foldid = foldid,
#                        weights = weights)
#       
#       mod <- glmnet(X[ ,s], y, 
#                     family = "binomial", 
#                     lambda = cvs$lambda.min,
#                     weights = weights)
#       
#       pred <- cbind(pred, predict(mod, X[,s], type = "response"))
#       mods[[paste0("comp", i)]] <- mod
#     }
#   }
#   
#   if (method == "bypf"){
#     
#     cvs <- cv.ipflas(X, y, 
#                      family = "binomial",
#                      foldid = foldid,
#                      blocks = blocks,
#                      weights = weights,
#                      penalty.factors = penalty.factors)
#     
#     for(i in seq_along(penalty.factors)){
#       
#       pf  <- paste0(penalty.factors[[i]], collapse = "-")
#       lam <- cvs$statlist[[pf]]$lambda[which.min(cvs$statlist[[pf]]$cvm)]
#       mod <- ipflas(X, y,
#                     family = "binomial",
#                     blocks = blocks,
#                     weights = weights, 
#                     pf = penalty.factors[[i]],
#                     lambda = lam)
#       
#       pred <- cbind(pred, predict(mod, X, type = "response"))
#       mods[[pf]] <- mod
#     }
#   }
#   
#   #- stacking
#   cvs <- cv.glmnet(pred, y, 
#                    family = "binomial", 
#                    foldid = foldid,
#                    weights = weights)
#   stk <- glmnet(pred, y, 
#                 family = "binomial", 
#                 lambda = cvs$lambda.min,
#                 weights = weights)
#   
#   alpha <- coef(stk)
#  
#   beta = mapply(function(obj){
#     b <- coef(obj)[,1]
#     b[b != 0]
#   }, obj = mods, SIMPLIFY = F)
#   
#   return(list(beta = beta, alpha = alpha))
# }

#-------------------------------------------------------------------------------
#- default output the prediction probability
# prediction for single Cox model

predBeta.glm <- function(beta, newdata, type = "response") {
  
  if (length(beta) == 0)
    link <- rep(1, nrow(newdata))
  else
    link <- cbind(1, newdata[, names(beta[-1]), drop = F]) %*% beta
  
  if (type == "response")
    out <-  1 / (1 + exp(-link))
  if (type == "link")
    out <- link
  
  return(out)
}

# predict multiple models
predBetas.glm <- function(betas, newdata, type = "response") {
  
  mapply(predBeta.glm, beta = betas, 
         MoreArgs = list(newdata = newdata, type = type), SIMPLIFY = FALSE)
  
}


predBLBin <- function(obj, newdata, ifCVMods = FALSE, type.BL = "response") {
  
  #
  if (ifCVMods) {
    
    cvMods <- obj$CVMods
    
    if (is.null(cvMods)) 
      stop("ifCVMods = TRUE: need CV models. Try set `track = TRUE` in stack.cox or StkCoxBL")
    
    # initial preds
    BLBin = 0
    for(i in seq_along(cvMods)){
      
      betas <- cvMods[[i]]$mods
      pred <- predBetas.glm(betas = betas, newdata = newdata, type = type.BL)
      
      # add predictions, get average when divided by n
      # is it reasonable to average the probabilites 
      BLBin <- mapply(`+`, x = pred, y = BLBin, SIMPLIFY = FALSE)
    }
    
    # get the averaged results
    BLBin <- lapply(BLBin, `/`, length(cvMods))
  }
  
  
  if (!ifCVMods){
    BLBin <-  predBetas.glm(betas = obj$beta, newdata = newdata, type = type.BL)
  }

  BLBin <- as.data.frame(BLBin)
  
  class(BLBin) <- c("predBLBin", "data.frame")
  return(BLBin)
}

  
predSLBin <- function(BLBin, alphas, type = c("link", "response")) {
  
  
  if (!inherits(BLBin, "predBLBin"))
    stop("obj should be object from predBLBin()")
  
  # alphas should be a list
  if (!is.list(alphas))
    alphas <- list(alphas)
  
  # define the out, should be a list
  predList <- list()
  BLBin <- as.matrix(BLBin)
  
  for(j in seq_along(alphas)){
    
    #- get predprob super model 
    w <- alphas[[j]]
    
    truncValue <- attr(w, "truncValue")
    transPFun  <- attr(w, "transPFun")
    
    # set truncValue to be zero
    Z <- truncTrans(BLBin, truncValue, transPFun)
      
    if (length(w) == (ncol(Z) + 1)) {
      pred <- cbind(1, Z) %*% w
    } 
    
    if (length(w) == (ncol(Z))) {
      pred <- Z %*% w
    }
    
    if (type == "response")
      pred <- 1 / (1 + exp(-pred))
    
    predList[[names(alphas)[j]]] <- pred
  } 
  
  return(predList)
}


predict.stk.bin <- function(obj, newdata, 
                            type = c("response", "link"), 
                            type.BL = "response",
                            ifCVMods = FALSE){

  BL <- predBLBin(obj = obj, newdata = newdata, ifCVMods = ifCVMods, type.BL = type.BL) 
  
  # super learner
  SL <- predSLBin(BL, alphas = obj$alphas, type = type)
  return(SL)
}

# merge submodels for both stk.cox and stk.glm

merge.stk_glm <- function(..., subset = c()) {
  
  mods <- list(...)
  # if vector, trans to list
  if (!is.list(subset)) {
    subset = replicate(length(mods), subset, simplify = FALSE) 
  }
  
  CVpred <- c()
  beta   <- list()
  method <- c()
  
  ifCVMods <- "CVMods" %in% names(mods[[1]])
  
  for(i in seq_along(mods)){
    
    ind <- subset[[i]]
    if(is.null(ind))
      ind <- seq_len(length(mods[[i]][["beta"]]))
    
    # if duplicated name
    CVpred <- cbind(CVpred, mods[[i]][["CVpred"]][ ,ind, drop = FALSE])
    beta   <- c(beta, mods[[i]][["beta"]][ind])
    method <- c(method, mods[[i]][["method"]])
    
    if (ifCVMods) {
      if(i == 1) {
        CVMods <- mods[[1]][["CVMods"]]
        
        # selection
        for(j in seq_along(CVMods)){
          CVMods[[j]][["mods"]][-ind] <- NULL
        }
        
      } else {
        # tmp
        CVMods_tmp <- mods[[i]][["CVMods"]]
        
        # selection
        for(j in seq_along(CVMods_tmp)){
          CVMods_tmp[[j]][["mods"]][-ind] <- NULL
        }
        
        # merge
        CVMods <- mapply(function(x, y) {
          x$mods <- c(x$mods, y$mods)
          return(x)
        }, x = CVMods, y = CVMods_tmp, SIMPLIFY = FALSE)
      }
    }
  }
  
  m <- mods[[1]]
  
  # rename to aviod duplicated names
  nams <- colnames(CVpred)
  if (any(duplicated(nams)))
    nams <- make.names(nams, unique = TRUE)
  
  colnames(CVpred) <- nams
  m$CVpred <- CVpred 
  
  names(beta) <- nams
  m$beta   <- beta

  m$method <- method
  if (ifCVMods){
    names(CVMods) <- nams
    m$CVMods <- CVMods
  }
  
  return(m)
}


merge.stk_cox <- function(..., subset = c()) {
  
  mods <- list(...)
  # if vector, trans to list
  if (!is.list(subset)) {
    subset = replicate(length(mods), subset, simplify = FALSE) 
  }
  
  # if using partial cox model
  ifCox <- "CVpredCox" %in% names(mods[[1]])
  ifCVMods <- "CVMods" %in% names(mods[[1]])
  
  CVpred <- mods[[1]][["CVpred"]][ ,c(1, 2), drop = FALSE]
  beta   <- list()
  lp     <- list()
  method <- c()
  
  CVpredCox <- c()
  
  for(i in seq_along(mods)){
    
    ind <- subset[[i]]
    if(is.null(ind))
      ind <- seq_len(length(mods[[i]][["subMod"]][["beta"]]))
    
    CVpred <- cbind(CVpred, mods[[i]][["CVpred"]][ ,ind + 2, drop = FALSE])
    beta   <- c(beta, mods[[i]][["subMod"]][["beta"]][ind])
    lp     <- c(lp, mods[[i]][["subMod"]][["lp"]][ind])
    
    if (ifCox)
      CVpredCox <- cbind(CVpredCox, mods[[i]][["CVpredCox"]][ ,ind, drop = FALSE])
    
    if (ifCVMods) {
      if(i == 1) {
        CVMods <- mods[[1]][["CVMods"]]
        
        # selection
        for(j in seq_along(CVMods)){
          CVMods[[j]][["mods"]][-ind] <- NULL
        }
        
      } else {
        # tmp
        CVMods_tmp <- mods[[i]][["CVMods"]]
        
        # selection
        for(j in seq_along(CVMods_tmp)){
          CVMods_tmp[[j]][["mods"]][-ind] <- NULL
        }
        
        # merge
        CVMods <- mapply(function(x, y) {
          x$mods <- c(x$mods, y$mods)
          return(x)
        }, 
        x = CVMods, y = CVMods_tmp, SIMPLIFY = FALSE)
      }
    }
  }
  
  
  mod <- mods[[1]]
  mod$subMod$beta <- beta
  mod$subMod$lp   <- lp
  mod$CVpred      <- CVpred
  if (ifCox)
    mod$CVpredCox  <- CVpredCox
  
  if (ifCVMods)
    mod$CVMods  <- CVMods
  
  return(mod)
}

# merge submodels, will be replaced with S3 function c.stk()

mergeSubMods <- function(..., subset = c(), ifCox = FALSE) {
  
  message("Consider using the new S3 funcition: merge()")
  
  # for Cox, the first two columns are time and events
  # subset, vector or list, the submodels to be included
  
  mods <- list(...)
  # if vector, trans to list
  if (!is.list(subset)) {
    subset = replicate(5, subset, simplify = FALSE) 
  }
  
  CVpred <- c()
  beta   <- list()
  method <- c()
  
  for(i in seq_along(mods)){
    ind <- subset[[i]]
    if(is.null(ind))
      ind <- seq_len(length(mods[[i]][["beta"]]))
    if (ifCox)
      ind <- setdiff(ind, c(1, 2))
    
    CVpred <- cbind(CVpred, mods[[i]][["CVpred"]][ ,ind, drop = FALSE])
    beta   <- c(beta, mods[[i]][["beta"]][ind])
    method <- c(method, mods[[i]][["method"]])
  }
  
  m <- mods[[1]]
  m$CVpred <- CVpred 
  m$beta   <- beta
  m$method <- method
  
  return(m)      
}

#   beta = obj$beta
#   alpha = obj$alpha
# 
#   pred <- mapply(function(b){
#     s <- names(b)[-1]
#     p <- cbind(1, newdata[,s, drop = FALSE]) %*% b
#     1 / (1 + exp(-p))
#   }, b = beta, SIMPLIFY = T)
# 
#   # pred
#   
#   pred <- cbind(1, pred) %*% alpha
# 
#   if (type == "link")
#     return(pred)
#   if (type == "response")
#     return(1 / (1 + exp(-pred)))
# }

# #- use stacking combine pcr and pf
# 
# stack.pcr.binomial <- function(X, y,
#                                foldid = NULL,
#                                foldid.internal = NULL,
#                                seed = NULL,
#                                nfold = 10,
#                                blocks = rep(1, ncol(X)),
#                                ncomp,
#                                blocks.pc = 1,
#                                penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8))){
# 
#   n = length(y)
# 
#   if (!is.null(seed))
#     set.seed(seed)
# 
#   if (is.null(foldid)){
#     foldid = sample(rep(1:nfold, length = n), n)
#   }
# 
#   nfold = length(unique(foldid))
# 
#   if (is.null(foldid.internal)){
#     foldid.internal = list()
#     for (i in 1:nfold){
#       ns <- sum(foldid != i)
#       foldid.internal[[i]] <-  sample(rep(1:10, length = ns), ns)
#     }
#   }
# 
#   pred1 <- pred2 <- matrix(NA, nrow = n, ncol = length(penalty.factors))
# 
#   obj <- get.pcrs(X,
#                   ncomp = min(table(y)),
#                   blocks = blocks,
#                   blocks.pc = blocks.pc)
#   M  <- obj$M
# 
#   for(i in seq_len(nfold)) {
# 
#     s <- foldid != i
#     mods1 <- pf.mods(X[s,], y[s],
#                      family = "binomial",
#                      foldid = foldid.internal[[i]],
#                      blocks = blocks,
#                      penalty.factors = penalty.factors)
# 
#     tmp <- Reduce(cbind, predict.mods(mods1, X[!s,]))
#     pred1[!s,] <- tmp
# 
#     mods2 <- pf.mods(M[s,], y[s],
#                      family = "binomial",
#                      foldid = foldid.internal[[i]],
#                      blocks = obj$blocks,
#                      penalty.factors = penalty.factors)
# 
#     tmp <- Reduce(cbind, predict.mods(mods2, M[!s,]))
#     pred2[!s,] <- tmp
#     pred <- cbind(pred1, pred2)
#   }
# 
#   stack <- glmnet(pred, y, family = "binomial", lambda = 0, lower.limits = 0)
#   alpha <- setNames(coef(stack)[,1], c("(Intercept)",
#                                        paste0("pf(", names(mods1), ")"),
#                                        paste0("pcr(", names(mods2), ")")))
# 
#   #- refit the models
# 
#   #- ipf
#   mods1 <- pf.mods(X, y,
#                    family = "binomial",
#                    foldid = foldid,
#                    blocks = blocks,
#                    penalty.factors = penalty.factors)
# 
#   #- pcr
#   mods2 <- pf.mods(M, y,
#                    family = "binomial",
#                    foldid = foldid,
#                    blocks = obj$blocks,
#                    penalty.factors = penalty.factors)
# 
#   mods1 = mapply(function(obj){
#     b <- coef(obj)[,1]
#     b[b != 0]
#   }, obj = mods1, SIMPLIFY = F)
# 
#   mods2 <- mapply(function(mod, e, newblocks){
# 
#     a0 = mod$a0
#     beta = mod$beta[,1, drop = T]
#     beta = split(beta, newblocks)
# 
#     for(i in seq_along(beta)){
#       if (!is.null(e[[i]]))
#         beta[[i]] <- drop(e[[i]] %*% beta[[i]])
#     }
# 
#     beta <- Reduce(c, beta)
#     c(a0, beta[colnames(X)])
# 
#   }, mod = mods2,
#   MoreArgs = list(e = obj$eigin.vectors,
#                   newblocks = obj$blocks), SIMPLIFY = F)
# 
#   mods <- c(mods1, mods2)
# 
#   return(list(beta = mods, alpha = alpha))
# }

#-------------------------------------------------------------------------------
#- Stacking: clin, mod, las, ipf, pcr.las
#-------------------------------------------------------------------------------
# remove "pcrlas" "ipflas"
# default clin is i1 in blocks
# pcr, denotes the pcr only on molecular part
# fs, denotes the fs only on molecular part
# 
# sub.mods <- function(X, y,
#                      methods = c("clin", "mol", "naive", "stk.pf", 
#                                  "pcrmol", "fsmol"),
#                      ncomp,
#                      family,
#                      foldid,
#                      blocks = rep(1, ncol(X)),
#                      weights = rep(1, nrow(X)),
#                      block.clin,
#                      block.mol,
#                      penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8))) {
#   
# 
#   b <- unique(blocks)
#   
#   #- independent models
#   mods <- list()
#   
#   if ("clin" %in% methods) {
#     
#     s   <- blocks == block.clin
#     cvs <- cv.glmnet(X[,s], y, 
#                      family = family, 
#                      foldid = foldid,
#                      weights = weights)
#     tmp <- glmnet(X[,s], y, 
#                   family = family, 
#                   lambda = cvs$lambda.min,
#                   weights = weights)
#     beta <- coef(tmp)[,1]
#     beta <- beta[beta != 0]
#     
#     mods[["clin"]] <- list(beta = beta,
#                            lp   = predict(tmp, X[ ,s]),
#                            block = block.clin)
#   }
#   
#   # molecular variables
#   if ("mol" %in% methods) {
#     
#     s   <- blocks == block.mol
#     cvs <- cv.glmnet(X[ ,s], y, 
#                      family = family, 
#                      foldid = foldid,
#                      weights = weights)
#     tmp <- glmnet(X[,s], y, 
#                   family = family, 
#                   lambda = cvs$lambda.min,
#                   weights = weights)
#     beta <- coef(tmp)[,1]
#     beta <- beta[beta != 0]
#     
#     mods[["mol"]] <- list(beta = beta,
#                           lp   = predict(tmp, X[ ,s]),
#                           block = block.mol)
#   }
#   
#   
#   #- naive lasso
#   
#   if ("naive" %in% methods){
#     
#     cvs <- cv.glmnet(X, y, 
#                      family = family, 
#                      foldid = foldid,
#                      weights = weights)
#     tmp <- glmnet(X, y, 
#                   family = family, 
#                   lambda = cvs$lambda.min,
#                   weights = weights)
#     beta <- coef(tmp)[,1]
#     beta <- beta[beta != 0]
#     mods[["naive"]] <- list(beta = beta, 
#                             lp   = predict(tmp, X),
#                             block = b)
#   }
#   
#   # pcrmol
#   
#   if ("pcrmol" %in% methods){
#     
#     s = blocks == block.mol
#     cv.pcr <- cv.drpcr(X[ ,s], y,
#                        alpha = 1,
#                        family = family,
#                        ncomp = ncomp, # number of comp for principle components
#                        blocks = blocks[s],
#                        blocks.pc = block.mol,
#                        weights = weights, 
#                        penalty.factors = list(c(1)),
#                        foldid = foldid)
#     
#     pcr.mol <- drpcr(X[ ,s], y,
#                      alpha = 1,
#                      ncomp = ncomp,
#                      family = family,
#                      blocks = blocks[s],
#                      blocks.pc = block.mol,
#                      weights = weights, 
#                      pf = cv.pcr$pf.min,
#                      lambda = cv.pcr$lambda.min)
#     
#     if (family == "cox"){
#       beta = setNames(drop(pcr.mol$back.mod$beta), pcr.mol$variables)
#       lp   = X[,s] %*% beta
#     } else {
#       beta = setNames(drop(rbind(pcr.mol$back.mod$a0, 
#                                  pcr.mol$back.mod$beta)),
#                       c("(Intercept)", pcr.mol$variables))
#       lp  = cbind(1, X[,s]) %*% beta 
#     }
#     
#     mods[["pcrmol"]] <- list(beta = beta,
#                              beta.pcr = c(pcr.mol$a0, pcr.mol$beta[,1]),
#                              centers  = pcr.mol$centers,
#                              scales   = pcr.mol$scales,
#                              variables = pcr.mol$variables,
#                              lp = lp)
#   }
#   
#   
#   # fsmol
#   
#   a <- Sys.time()
#   
#   if ("fsmol" %in% methods){
#     
#     s = blocks == block.mol
#     cv.fs <- cv.drfs(X[ , s], y, 
#                      family = family, 
#                      foldid = foldid,
#                      blocks = rep(1, ncol(X[,s])), 
#                      blocks.screen = 1,
#                      cuts = seq(0.1, 1, 0.1), 
#                      weights = weights,
#                      penalty.factors = list(c(1)))
#     
#     fs.mol <- drfs(X[ , s], y,
#                    # alpha = 1,
#                    family = family,
#                    feature.order = cv.fs$feature.order,
#                    cut = cv.fs$cut.min,
#                    blocks = rep(1, ncol(X[,s])),
#                    blocks.screen = 1,
#                    weights = weights, 
#                    pf = cv.fs$pf.min,
#                    lambda = cv.fs$lambda.min)
#     
#     beta <- coef(fs.mol)[ ,1]
#     beta <- beta[beta != 0]
#     mods[["fsmol"]] <- list(beta = beta, 
#                             lp   = predict(fs.mol, 
#                                            newx = as.matrix(X[ ,rownames(fs.mol$beta)])),
#                             block = block.mol,
#                             cut = cv.fs$cut.min)
#   }
#   
#   
#   # stk.pf
#   
#   if ("stk.pf" %in% methods) {
#     
#     pf.mod <- pf.mods(X, y,
#                       family = family,
#                       foldid = foldid,
#                       blocks = blocks,
#                       weights = weights,
#                       penalty.factors = penalty.factors)
#     
#     tmp <- lapply(pf.mod, function(obj) {
#       
#       beta <- coef(obj)[,1]
#       beta <- beta[beta != 0]
#       
#       lp <- predict(obj, newx = X)
#       list(beta = beta,
#            lp = lp)
#     })
#     
#     mods <- append(mods, tmp) 
#   }
#   
#   return(list(mods = mods, family = family, y = y))
# }
# 
# 
# predict.sub.mods <- function(obj, newdata, type = c("link", "response")){
#   
#   pred.beta <- function(mod, newdata, type){
#     
#     sel  <- intersect(names(mod$beta), colnames(newdata))
#     newdata <- newdata[,sel, drop = F]
#     
#     #- center for pcr
#     if ("centers" %in% names(mod)){
#       cts  <- mod$centers[match(sel, mod$variables)]
#     } else {
#       cts <- rep(0, length(sel))
#     }
#     
#     if ("scales" %in% names(mod)){
#       sds <-  mod$scales[match(sel, mod$variables)]
#     } else {
#       sds <- rep(1, length(sel))
#     }
#     
#     #- in case beta is NULL
#     if(length(mod$beta) == 0)
#       beta = 0
#     else
#       beta = mod$beta
#     pred <- cbind(1, newdata) %*% beta
#     
#     if (type == "link") {
#       pred 
#     }
#     
#     if (type == "response"){
#       1 / (1 + exp(-pred))
#     }
#   }
#   
#   re <- mapply(pred.beta, mod = obj$mods, 
#                MoreArgs = list(newdata = newdata, 
#                                type = type))
#   return(re)
# }
# 
# stack.super.binomial <- function(X, y,
#                                  methods = c("clin", "mol", "naive", "stk.pf", 
#                                              "pcrmol", "fsmol"),
#                                  foldid = NULL,
#                                  foldid.internal = NULL,
#                                  nfold = 10,
#                                  blocks = rep(1, ncol(X)),
#                                  weights = rep(1, nrow(X)),
#                                  ncomp,
#                                  block.clin,
#                                  block.mol,
#                                  optim = "glmnet",
#                                  penalty.factors = list(c(1, 1), c(1, 2),
#                                                         c(1, 4), c(1, 8))){
#   
#   n = length(y)
#   # if(is.null(ncomp))
#   #   ncomp <- pmin(table(y))
# 
#   if (is.null(foldid)){
#     foldid = sample(rep(1:nfold, length = n), n)
#   }
# 
#   nfold = length(unique(foldid))
# 
#   if (is.null(foldid.internal)) {
#     foldid.internal = list()
#     for (i in 1:nfold){
#       ns <- sum(foldid != i)
#       foldid.internal[[i]] <- sample(rep(1:10, length = ns), ns)
#     }
#   }
#   
#   npred <- ifelse("stk.pf" %in% methods, 
#                   length(penalty.factors) + length(methods) - 1,
#                   length(methods))
#   pred <- matrix(NA, nrow = n, ncol = npred)
#   
#   for(i in seq_len(nfold)) {
#     
#     s <- foldid != i
#     mods <- sub.mods(X[s, ], y[s],
#                      family = "binomial",
#                      methods = methods,
#                      foldid = foldid.internal[[i]],
#                      blocks = blocks,
#                      weights = weights[s], 
#                      ncomp = ncomp,
#                      block.clin = block.clin,
#                      block.mol = block.mol,
#                      penalty.factors = penalty.factors)
#     
#     tmp <- predict.sub.mods(mods, X[!s,], type = "response")
#     pred[!s,] <- tmp
#   }
#   colnames(pred) <- names(mods$mods)
#   
#   if (optim == "glmnet"){
#     stack <- glmnet(pred, y, 
#                     family = "binomial", 
#                     weights = weights, 
#                     lambda = 0, 
#                     lower.limits = 0)
#     alpha <- coef(stack)
#   }
#   
#   if (optim == "L1"){
#     cv.stack <- cv.glmnet(pred, y, 
#                           family = "binomial", 
#                           lower.limits = 0,
#                           foldid = foldid,
#                           weights = weights)
#     stack    <- glmnet(pred, y, 
#                        family = "binomial", 
#                        weights = weights,
#                        lower.limits = 0, 
#                        lambda = cv.stack$lambda.min)
#     alpha <- coef(stack)
#   }
#   
#   # note, weights not support
#   if (optim == "solnp"){
#     alpha <- optimBi(pred, y)
#   }
#   
#   #- refit the models
#   mods <- sub.mods(X, y,
#                    family = "binomial",
#                    methods = methods,
#                    foldid = foldid,
#                    blocks = blocks,
#                    weights = weights, 
#                    ncomp = ncomp,
#                    block.clin = block.clin,
#                    block.mol = block.mol,
#                    penalty.factors = penalty.factors)
#   
#   mods$alpha = alpha
#   mods$optim = optim
#   return(mods)
# }
# 
# 
# predict.stack.super.binomial <- function(obj, newdata,
#                                          type =  c("link", "response")){
# 
#   dat <- predict.sub.mods(obj, newdata, type = "response")
#   
#   # if (obj$optim == "solno")
#   #   pred <- dat %*% obj$alpha
#   # else
#     pred <- cbind(1, dat) %*% obj$alpha
# 
#   if (type == "link"){
#     return(pred)
#   }
# 
#   if (type == "response"){
#     return( 1 / (1 + exp(-pred)) )
#   }
# }

# #- Stk super with postCV
# stack.super.binomial.postCV <- function(X, y,
#                                         methods = c("comp-wise", "naive",
#                                                     "pcr.las", "ipflasso"),
#                                         foldid = NULL,
#                                         nfold = 10,
#                                         blocks = rep(1, ncol(X)),
#                                         ncomp,
#                                         blocks.pc = 1,
#                                         penalty.factors = list(c(1, 1), c(1, 2),
#                                                                c(1, 4), c(1, 8))){
#   
#   
#   n = length(y)
#   
#   if (is.null(foldid)){
#     foldid = sample(rep(1:nfold, length = n), n)
#   }
#   
#   mods <- sub.mods(X, y,
#                    family = "binomial",
#                    methods = methods,
#                    foldid = foldid,
#                    blocks = blocks,
#                    ncom = ncomp,
#                    penalty.factors = penalty.factors)
#   
#   pred <- predict.sub.mods(mods, X, type = "response")
#   
#   colnames(pred) <- names(mods)
#   cvs <- cv.glmnet(pred, y, family = "binomial", 
#                    foldid = foldid)
#   
#   stack <- glmnet(pred, y, family = "binomial", 
#                   lambda = cvs$lambda.min)
#   alpha <- coef(stack)
#   
#   return(list(beta = mods, alpha = alpha))
# }
# 
#-------------------------------------------------------------------------------


















