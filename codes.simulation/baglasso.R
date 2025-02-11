

# bagging lasso

baglasso <- function(X, y, 
                     family,
                     weights = rep(1, nrow(X)),
                     bootID = NULL, 
                     randomP = ncol(X), 
                     nboot = 100,
                     fastTune = 100, 
                     aggFun = mean,
                     ...) {
    require(dplyr)
    require(survival)
    
    # argument checking
    if (fastTune > nboot) {
        fastTune = nboot
        warning(paste0("fastTune is reduced to ", nboot))
    }
    
    n = nrow(X)
    p = ncol(X)
    
    if (is.null(bootID)) {
        bootID <- lapply(seq_len(nboot), function(x) sample(seq_len(n), n, replace = TRUE))
    }
    
    CVMods <- list()
    maxLam <- c()
    lambdas <- NULL
    lambdaSeq <- NULL
    pred <- c()
    
    for (i in seq_along(bootID)) {
        
        idSamples <- bootID[[i]]
        idFeatures <- sample(seq_len(p), randomP)
        
        if (is.null(lambdaSeq)) {
            cvs <- cv.glmnet(x = X, y = y, family = family, weights = weights, ...)
            lambdas <- c(lambdas, cvs$lambda.min)
        }
        
        las <- glmnet(X, y, family = family, weights = weights, lambda = lambdaSeq, ...)
        tmplambdaSeq <- las$lambda
        las <- extractGlmnet(las, lambda = rev(lambdas)[1])
        
        if (i <= fastTune) {
            lambdas <- c(lambdas, las$lambda)
            maxLam <- c(maxLam, max(tmplambdaSeq[is.finite(tmplambdaSeq)], na.rm = TRUE))
        } else if (i == (fastTune + 1)) {
            lambdas <- mean(lambdas)
            lambdaSeq <- exp(seq(log(max(maxLam)), log(lambdas), length.out = 50))
        }
        
        tmpBeta <- coef(las)[, 1]
        simMod <- tmpBeta[tmpBeta != 0]
        CVMods[[i]] <- list(mods = simMod, trainID = idSamples)
        
        idSamplesRest <- setdiff(seq_len(n), idSamples)
        pred <- rbind(pred, cbind(predict(las, X[idSamplesRest, ]), id = idSamplesRest))
    }
    
    pred <- as.data.frame(pred) %>%
        group_by(id = as.numeric(pred[, "id"])) %>%
        summarise_all(aggFun) %>%
        arrange(desc(id)) %>%
        ungroup()
    
    id <- as.numeric(pred$id)
    if (length(setdiff(seq_len(n), id)) > 0) {
        cat(paste0("Cases were not included for super learner: ", paste0(setdiff(seq_len(n), id), collapse = ",")))
    }
    pred <- pred %>% dplyr::select(-id)
    attr(pred, "id") <- id
    
    out <- list(pred = pred,
                y = y,
                weights = weights[id],
                CVMods = CVMods,
                family = family,
                aggFun = aggFun)
    
    class(out) <- c("bagLasso", "list")
    return(out)
}

pred.baglasso <- function(obj, newdata, type = "link") {
    
    cvMods <- obj$CVMods
    
    BL <- sapply(cvMods, function(mod) {
        pred <- predBeta.glm(beta = mod$mods, newdata = newdata, type = type)
        return(pred)
    })
    
    BL <- apply(BL, 1, obj$aggFun)
    attr(BL, "type") <- type
    
    BL <- as.data.frame(BL)
    class(BL) <- c("predbaglasso", "data.frame")
    
    return(BL)
}


#------------------------------------------------------------------------------


































