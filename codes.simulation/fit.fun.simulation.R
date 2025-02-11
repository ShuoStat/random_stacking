
simData <- function(n.obs, t.gene, r, sd_sett, beta){
    
    # 1, generate data
    n.clin <- c(3, 3, 2, 2, 5)
    n.gene <- rep(10, 49)
    
    rho.c <- rep(r, 5)
    rho.b <- rep(r, 5)
    rho.g <- rep(r, 49)
    
    d <- generaClinMol(n.obs = n.obs, tot.genes = t.gene, n.groups = 49, 
                       n.clin = n.clin, n.gene = n.gene, mean.n.gene = 15, 
                       mu.g = 6, sigma.g = 0.65, mu.c = 1, sigma.c = 0.5, 
                       rho.c = rho.c, rho.b = rho.b, rho.g = rho.g, phi = 0.1, 
                       nu = 10, tau = 20)
    
    # 2, 
    X <- cbind(d$clin, d$gene)
    y <- X %*% beta + rnorm(n.obs, sd = sd_sett)
    y <- as.vector(y)
    blocks = as.numeric(grepl("gene", colnames(X))) + 1
    
    return(list(X = X, y = y, blocks = blocks))
}

# 3, stacking
fit.basic.sim <- function(X, y, blocks) {
    
    n <- nrow(X)
    nfold <- 10
    foldid = sample(rep(1:nfold, length = n), n)
    
    # naive
    t1 <- Sys.time()
    cvs <- cv.glmnet(X, y, family = "gaussian",
                     foldid = foldid,
                     weights = rep(1, nrow(X)))
    naive <- glmnet(X, y, family = "gaussian",
                    lambda = cvs$lambda.min,
                    weights = rep(1, nrow(X)))
    
    t2 <- Sys.time()
    beta <- coef(naive)[,1]
    beta <- beta[beta != 0]
    naive <- list(beta = beta, time = t2 - t1)
    
    
    # bagging
    t1 <- Sys.time()
    baglas <- baglasso(X = X, y = y, 
                       family = "gaussian",
                       weights = rep(1, nrow(X)),
                       bootID = NULL, 
                       randomP = ncol(X), 
                       nboot = 100,
                       fastTune = 20, 
                       aggFun = mean) 
    t2 <- Sys.time()
    baglas$time <- t2 - t1
    
    
    # ipflasso
    t1 <- Sys.time()
    cv.ipf <- cv.ipflas(X, y, family = "gaussian",
                        foldid = foldid,
                        blocks = blocks,
                        weights = rep(1, nrow(X)),
                        alpha = 1,
                        penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)))
    
    ipf.las <- ipflas(X, y, family = "gaussian",
                      weights = rep(1, nrow(X)),
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
    
    list(naive = naive,
         baglas = baglas,
         ipf.las = ipf.las)
    
}

# stacking

fit.stk.sim <- function(X, y, blocks) {
    
    n <- nrow(X)
    nfold <- 10
    
    foldid = sample(rep(1:nfold, length = n), n)
    foldid.internal = list()
    for (i in 1:nfold){
        ns <- sum(foldid != i)
        foldid.internal[[i]] <-  sample(rep(1:10, length = ns), ns)
    }
    
    t1 <- Sys.time()
    stk.comp <- stk(X, y,
                    family = "gaussian", 
                    foldid = foldid,
                    foldid.internal = foldid.internal,
                    seed = NULL,
                    nfold = 10,
                    method = "bycomp",
                    blocks = blocks,
                    weights = rep(1, nrow(X)),
                    penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                    track = FALSE, 
                    optimFun = c("lasso", "ridge", "solnp"),
                    moreOptimFun = NULL, 
                    foldidSL = NULL,
                    nfoldSL = 10,
                    limits = 0,
                    ifWeights = TRUE,
                    ifshrinkage = TRUE)
    t2 <- Sys.time()
    stk.comp$time <- t2 - t1
    
    
    t1 <- Sys.time()
    stk.pf <- stk(X, y,
                  family = "gaussian", 
                  foldid = foldid,
                  foldid.internal = foldid.internal,
                  seed = NULL,
                  nfold = 10,
                  method = "bypf",
                  blocks = blocks,
                  weights = rep(1, nrow(X)),
                  penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                  track = TRUE, 
                  optimFun = c("lasso", "ridge", "solnp"),
                  moreOptimFun = NULL, 
                  foldidSL = NULL,
                  nfoldSL = 10,
                  limits = 0,
                  ifWeights = TRUE,
                  ifshrinkage = TRUE)
    t2 <- Sys.time()
    stk.pf$time <- t2 - t1
    
    return(list(stk.comp = stk.comp,
                stk.pf = stk.pf))
    
}


fit.random.stk.sim <- function(X, y, blocks) {
    
    # random stacking by comp
    t1 <- Sys.time()
    Rstk.comp <- random.stack(X, y,
                              family = "gaussian",
                              bootID = NULL, 
                              randomP = table(blocks) / c(1, 3), 
                              nboot = 100,
                              fastTune = 20,
                              foldid = NULL,
                              nfold = 10,
                              blocks = blocks,
                              weights = rep(1, nrow(X)),
                              method = "bycomp",
                              penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                              optimFun = c("lasso", "ridge", "solnp"),
                              moreOptimFun = NULL,
                              foldidSL = NULL,
                              nfoldSL = 10,
                              limits = 0,
                              ifWeights = TRUE,
                              ifshrinkage = TRUE, 
                              # truncValue = 0, 
                              # transPFun = function(x) x,
                              randomSL = FALSE, 
                              ifRandomSL = FALSE)
    t2 <- Sys.time()
    Rstk.comp$time <- t2 - t1
    
    # random stacking by pf
    t1 <- Sys.time()
    Rstk.pf <- random.stack(X, y,
                            family = "gaussian",
                            bootID = NULL, 
                            randomP = table(blocks) / c(1, 3),  
                            nboot = 100,
                            fastTune = 20,
                            foldid = NULL,
                            nfold = 10,
                            blocks = blocks,
                            weights = rep(1, nrow(X)),
                            method = "bypf",
                            penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                            optimFun = c("lasso", "ridge", "solnp"),
                            moreOptimFun = NULL,
                            foldidSL = NULL,
                            nfoldSL = 10,
                            limits = 0,
                            ifWeights = TRUE,
                            ifshrinkage = TRUE, 
                            # truncValue = 0, 
                            # transPFun = function(x) x,
                            randomSL = FALSE, 
                            ifRandomSL = FALSE)
    t2 <- Sys.time()
    Rstk.pf$time <- t2 - t1
    
    list(Rstk.comp = Rstk.comp,
         Rstk.pf = Rstk.pf)
    
}


























