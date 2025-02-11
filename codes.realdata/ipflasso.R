
# debug in get.lambda, 28/09/2024

get.lambda <- function(X, y, family, alpha,
                       intercept = T,
                       is.offset = F,
                       offset = rep(0, nrow(X)),
                       weights = rep(1, nrow(X)),
                       exclude = c(),
                       pfs = rep(1, ncol(X)),
                       nlam = 100,
                       lambda.min.ratio = 0.01) {

  # check and remove variables with equal values
  # check by rows to save time
  indCols <- seq_len(ncol(X))
  for(i in seq_len(nrow(X) - 1)){
    # ind <- order(X[1,])
    indCols <- indCols[X[i, indCols] == X[i + 1, indCols]]
    if (length(indCols) == 0)
      break()
  }
  
  # system.time(apply(X, 2, function(x) length(unique(x)) == 1))
  # remove columns of equal values
  if (length(indCols) > 0) {
    X <- X[ ,-indCols]
    if (is.list(pfs))
      pfs <- lapply(pfs, function(x) x[-indCols])
    else
      pfs <- pfs[-indCols]
    exclude <- setdiff(exclude, indCols)
  }

  family = do.call(family, list())
  nobs <- nrow(X); nvars <- ncol(X)

  mysd <- function(x) sqrt(sum((x-mean(x))^2) / length(x))
  sx <- scale(X, scale = apply(X, 2, mysd))
  sx <- as.matrix(sx)

  # compute mu and null deviance
  # family = binomial() gives us warnings due to non-integer weights
  # to avoid, suppress warnings
  
  if (intercept) {
    if (is.offset) {
      suppressWarnings(tempfit = glm(y ~ 1, family = family,
                                     weights = weights, offset = offset))
      mu <- tempfit$fitted.values
    } else {
      mu <- rep(weighted.mean(y, weights), times = nobs)
    }
  } else {
    mu <- family$linkinv(offset)
  }
  
  # dev_function <- function(y, mu, weights, family) {
  #   sum(family$dev.resids(y, mu, weights))
  # }
  # # 
  # nulldev <- dev_function(y, mu, weights, family)
  # cat(nulldev)
  #- standardize pfs, pfs can be a list
  
  if (is.list(pfs)){
    pfs <- lapply(pfs, function(x){
      x = pmax(0, x)
      if (max(x) <= 0) stop("All penalty factors are <= 0")
      as.double(x * nvars / sum(x))
    })} else {
      pfs = pmax(0, pfs)
      if (max(pfs) <= 0) stop("All penalty factors are <= 0")
      pfs <- list(as.double(pfs * nvars / sum(pfs)))
    }

  lambda_max <- list()
  for(i in seq_along(pfs)){
    pf <- pfs[[i]]
    # if some penalty factors are zero, we have to recompute mu
    pf_zero <- setdiff(which(pf == 0), exclude)
    if (length(pf_zero) > 0) {
      tempx <- sx[ , pf_zero, drop = FALSE]

      if (intercept) {
        tempx <- cbind(1,tempx)
      }
      tempfit <- glm.fit(tempx, y, family = family, weights = weights,
                         offset = offset)
      mu <- tempfit$fitted.values
    }
    # compute lambda max
    ju <- rep(1, nvars)
    ju[exclude] <- 0 # we have already included constant variables in exclude
    r <- y - mu
    eta <- family$linkfun(mu)
    v <- family$variance(mu)
    m.e <- family$mu.eta(eta)
    weights <- weights / sum(weights)
    rv <- r / v * m.e * weights
    g <- abs(drop(t(rv) %*% sx))

    g <- g * ju / ifelse(pf > 0, pf, 1)
    g <- g[is.finite(g)]
    lambda_max[[i]] <- max(g, na.rm = TRUE) / max(alpha, 1e-3)
    
  }
  

  lambda_max <- unlist(lambda_max)
  lambda_max <- max(lambda_max, na.rm = TRUE)
  exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
          length.out = nlam))
}

#- recoded from get_cox_lambda_max is a function in glmnet package
get.lambda.cox <- function(X, y,
                           nlam = 100,
                           lambda.min.ratio = 0.01,
                           alpha,
                           weights = rep(1, nrow(X)),
                           offset = rep(0, nrow(X)),
                           exclude = c(),
                           pfs = rep(1, ncol(X))){

  # check and remove variables with equal values
  # check by rows to save time
  indCols <- seq_len(ncol(X))
  for(i in seq_len(nrow(X) - 1)){
    indCols <- indCols[X[i, indCols] == X[i + 1, indCols]]
    if (length(indCols) == 0)
      break()
  }
  # system.time(apply(X, 2, function(x) length(unique(x)) == 1))
  # remove columns of equal values
  # remove columns of equal values
  if (length(indCols) > 0) {
    X <- X[ ,-indCols]
    pfs <- pfs[-indCols]
    exclude <- setdiff(exclude, indCols)
  }

  mysd <- function(x) sqrt(sum((x-mean(x))^2) / length(x))
  sx <- scale(X, scale = apply(X, 2, mysd))
  sx <- as.matrix(sx)

  nobs <- nrow(X); nvars <- ncol(X)
  # extract strata (if any)
  if ("strata" %in% names(attributes(y))) {
    strata <- attr(y, "strata")
  } else {
    strata <- rep(1, nobs)
  }
  if (length(strata) != nobs) stop("length of strata != nobs")

  # if some penalty factors are zero, we need to compute eta

  if (is.list(pfs)){
    pfs <- lapply(pfs, function(x){
      x = pmax(0, x)
      if (max(x) <= 0) stop("All penalty factors are <= 0")
      as.double(x * nvars / sum(x))
    })}else{
      pfs = pmax(0, pfs)
      if (max(pfs) <= 0) stop("All penalty factors are <= 0")
      pfs <- list(as.double(pfs * nvars / sum(pfs)))
    }

  lambda_max <- list()

  for(i in seq_along(pfs)){

    pf <- pfs[[i]]
    pf_zero <- setdiff(which(pf == 0), exclude)

    if (length(pf_zero) > 0) {

      tempx <- sx[, pf_zero, drop = FALSE]
      eps <- glmnet.control()$epsnr
      if (length(unique(strata)) == 1) {
        fit <- survival::coxph(y ~ offset(offset) + tempx,
                               weights = weights,
                               eps = eps)
      } else {
        fit <- survival::coxph(y ~ offset(offset) + tempx + strata(strata),
                               weights = weights,
                               eps = eps)
      }
      eta <- predict(fit, reference = "sample") ## Coxph can do strata-specific centering
    } else {
      eta <- offset
    }

    # keep numbers small; partial likelihood independent of centering
    eta <- eta - mean(eta)
    ju  <- rep(1, nvars)
    # we have already included constant variables in exclude
    ju[exclude] <- 0

    # Get cox gradient at "null" point
    # Note that coxgrad already includes weights, so no need to include them
    # In subsequent computations
    null_grad <- glmnet::coxgrad(eta = eta, y = y, w = weights)
    g <- abs(drop(t(null_grad) %*% sx))
    g <- g * ju / ifelse(pf > 0, pf, 1)
    lambda_max[[i]] <- max(g) / max(alpha, 1e-3)
  }

  lambda_max <- max(unlist(lambda_max))
  
  exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio),
          length.out = nlam))
}

#- pool statlist

pool.stat <- function(statlist.ncv, nfolds, ncv){

  #- cvm
  cvms <- lapply(statlist.ncv, `[[`, "cvm")
  cvm <- Reduce(`+`, cvms) / length(cvms)

  #- cvsd
  pool.sd <- function(sds, nfolds, ncv, dev){

    ss <- (nfolds - 1) * sds^2
    pool.ss <- (sum(ss) + sum(dev^2 * nfolds)) / (nfolds * ncv - 1)
    sqrt(pool.ss)
  }

  cvsds <- lapply(statlist.ncv, `[[`, "cvsd")
  devs <- lapply(cvms, function(x) x - cvm)

  cvsd <- c()
  for (i in 1:length(cvm)) {
    sds <- unlist(lapply(cvsds, `[`, i))
    dev <- unlist(lapply(devs, `[`, i))
    cvsd <- c(cvsd, pool.sd(sds, nfolds, ncv, dev))
  }

  lambda = statlist.ncv[[1]]$lambda
  nz =  statlist.ncv[[1]]$nzero
  nas = is.na(cvsd)

  if (any(nas)) {
    lambda = lambda[!nas]
    cvm    = cvm[!nas]
    cvsd   = cvsd[!nas]
    nz     = nz[!nas]
  }

  list(lambda = lambda,
       cvm = cvm, cvsd = cvsd,
       cvup = cvm + cvsd, cvlo = cvm - cvsd,
       nzero = nz)
}

#' cv.ipflas
#'
#' @description
#'
#' Favoring strategy for lasso
#'
#' @param X matrix
#' @param y vector of response
#' @param family defaulted as gaussian
#' @param foldid list or vector depends on ncv
#' @param nfolds defaulted as 10
#' @param ncv defaulted as 5
#' @param alpha defaulted as 1
#' @param blocks vector of the same length of ncol(X)
#' @param penalty.factors list of penalty factors
#' @param type.measure "default"
#' @param lambda NULL
#' @param nlam 100
#' @param lambda.min.ratio 0.001
#' @param weights 1
#'
#' @examples
#'
#' @import glmnet
#' @export cv.ipflas


cv.ipflas <- function(X, y, family = "gaussian",
                      foldid = NULL,
                      nfolds = 10,
                      ncv = 1,
                      alpha = 1,
                      blocks = rep(1, ncol(X)),
                      penalty.factors = NULL,
                      type.measure = "default",
                      lambda = NULL,
                      nlam= 100,
                      lambda.min.ratio = 0.001,
                      weights = rep(1, nrow(X)),
                      ...){


  #- foldid
  if (is.null(foldid))
    foldid = sapply(1:ncv, function(x) sample(rep(1:nfolds, length = n), n),
                    simplify = F)

  if (!is.list(foldid))
    foldid <- list(foldid)

  #- updat nfolds, ncv
  nfolds = length(unique(foldid[[1]]))
  ncv = length(foldid)
  nblocks = length(unique(blocks))

  #- using default penalty factors
  if (is.null(penalty.factors)) {
    l <- list()
    for(i in 1:nblocks){
      l[[i]] <- c(1, 2, 4, 8, 16)
    }
    penalty.factors <- do.call("expand.grid", list(... = l))
    penalty.factors <- penalty.factors[apply(penalty.factors, 1, function(x) any(x == 1)),, drop = F]
    colnames(penalty.factors) <- paste0("B", 1:nblocks)
    #- transform to list
    penalty.factors <- split(penalty.factors,
                             rep(1:nrow(penalty.factors)),
                             ncol(penalty.factors))
  }

  if (any(unlist(lapply(penalty.factors, function(x) length(x) != nblocks))))
    stop("penalty.factors and blocks are not equvalent")

  if (is.null(lambda)){

    pfs <- list()
    for(i in seq_along(penalty.factors)){
      pfs[[i]] <- as.numeric(as.character(factor(blocks, unique(blocks),
                                                 penalty.factors[[i]])))
    }

    if (family == "cox"){
      lambda = get.lambda.cox(X, y,
                              nlam = nlam,
                              alpha = alpha,
                              lambda.min.ratio = lambda.min.ratio,
                              weights = weights,
                              pfs = pfs)
    } else {
      lambda = get.lambda(X, y,
                          family = family,
                          alpha = alpha,
                          nlam = nlam,
                          lambda.min.ratio = lambda.min.ratio,
                          weights = weights,
                          pfs = pfs)
    }
  }
  lambda <- sort(lambda, decreasing = T)

  statlist <- list()

  for(i in seq_along(penalty.factors)){
    pf <- as.numeric(as.character(factor(blocks, unique(blocks),
                                         penalty.factors[[i]])))
    statlist.ncv <- list()
    for(j in 1:ncv){
      cv <- cv.glmnet(x = X, y = y,
                      family = family,
                      foldid = foldid[[j]],
                      type.measure = type.measure,
                      alpha = alpha,
                      weights = weights,
                      penalty.factor = pf,
                      lambda = lambda, ...)

      statlist.ncv[[j]] <- cv[1:6]
    }
    #- summarize of ncv
    nam <- paste0(penalty.factors[[i]], collapse = "-")
    if (ncv > 1)
      statlist[[nam]] <- pool.stat(statlist.ncv, nfolds, ncv)
    else
      statlist[[nam]] <- statlist.ncv[[1]]
  }

  cvm <- as.matrix(Reduce(cbind, lapply(statlist, `[[`, "cvm")))
  cvsd <- as.matrix(Reduce(cbind, lapply(statlist, `[[`, "cvsd")))
  nzero <- as.matrix(Reduce(cbind, lapply(statlist, `[[`, "nzero")))

  if (type.measure %in% c("auc", "C"))
    cvm <- -cvm

  #- get ordered penalty.factors
  if (length(penalty.factors) == 1)
    message("Penalty factors should be ordered decresingly or increasingly")

  id.pf <- unlist(lapply(penalty.factors, function(x) min(x) / max(x)))
  id.pf <- order(id.pf, decreasing = T)

  cvm.pf <- cvm[, id.pf, drop = F]
  cvsd.pf <- cvsd[, id.pf, drop = F]
  nzero.pf <- nzero[,id.pf, drop = F]

  #- get position
  pos <- which(cvm.pf <= min(cvm.pf), arr.ind = T)

  #- favoring large lambda and equal penalty factor
  id.min <- pos[order(pos[,1], pos[,2], decreasing = F)[1], ]
  lambda.min <- lambda[id.min[1]]
  pf.min <- penalty.factors[[id.pf[id.min[2]]]]
  nzero.min <- nzero.pf[id.min[1], id.min[2]]

  #- 1sd rule
  cvm.1se <- cvm.pf[id.min[1], id.min[2]] + cvsd.pf[id.min[1], id.min[2]]
  pos <- which(cvm.pf <= cvm.1se, arr.ind = T)
  id.1se <- pos[order(pos[,1], pos[,2], decreasing = F)[1], ]
  lambda.1se <- lambda[id.1se[1]]
  pf.1se <- penalty.factors[[id.pf[id.1se[2]]]]
  nzero.1se <- nzero.pf[id.1se[1], id.1se[2]]

  #- summary index
  index = matrix(NA, ncol = 5, nrow = 2,
                 dimnames = list(c("min", "1se"),
                                 c("Lambda", "PF", "Measure", "SE", "Nonzero")))

  index[,1] <- round(c(lambda.min, lambda.1se), 2)
  index[,2] <- c(paste0(pf.min, collapse = "-"), paste0(pf.1se, collapse = "-"))
  index[,3] <- round(c(cvm.pf[id.min[1], id.min[2]], cvm.pf[id.1se[1], id.1se[2]]), 2)
  index[,4] <- round(c(cvsd.pf[id.min[1], id.min[2]], cvsd.pf[id.1se[1], id.1se[2]]), 2)
  index[,5] <- c(nzero.pf[id.min[1], id.min[2]], nzero.pf[id.1se[1], id.1se[2]])

  fit <- list(statlist = statlist,
              pf = penalty.factors,
              lambda.min = lambda.min,
              lambda.1se = lambda.1se,
              pf.min = pf.min,
              pf.1se = pf.1se,
              nzero.min = nzero.min,
              nzero.1se = nzero.1se,
              ncv = ncv,
              index = index)

  return(fit)
}


#' ipflas
#'
#' @description
#' Favoring strategy for lasso
#'
#' @param X matrix
#' @param y vector of response
#' @param family defaulted as gaussian
#' @param pf a vector
#' @param blocks vector of the same length of ncol(X)
#' @param lambda NULL
#' @param nlam 100
#' @param lambda.min.ratio 0.01
#' @param weights 1
#'
#' @examples
#'
#' @import glmnet
#' @export ipflas


ipflas <- function(X, y,
                   family = "gaussian",
                   alpha = 1,
                   blocks = rep(1, ncol(X)),
                   pf = NULL,
                   lambda = NULL,
                   nlam= 100,
                   lambda.min.ratio = 0.01,
                   weights = rep(1, nrow(X)),
                   ...){

  if (length(unique(blocks)) != length(pf))
    stop("Types of blocks should be identical to length of pf")

  pf <- as.numeric(as.character(factor(blocks, unique(blocks), pf)))

  if (is.null(lambda)){

    if (family == "cox"){
      lambda = get.lambda.cox(X, y,
                              alpha =alpha,
                              nlam = nlam,
                              lambda.min.ratio = lambda.min.ratio,
                              weights = weights,
                              pfs = pf)
    }else{
      lambda = get.lambda(X, y,
                          family = family,
                          alpha = alpha,
                          nlam = nlam,
                          lambda.min.ratio = lambda.min.ratio,
                          weights = weights,
                          pfs = pf)
    }
  }
  lambda <- sort(lambda, decreasing = T)

  fit <- glmnet(x = X, y = y,
                family = family,
                alpha = alpha,
                lambda = lambda,
                penalty.factor = pf,
                weights = weights,
                ...)

  return(fit)
}

#-end --------------------------------------------------------------------------




