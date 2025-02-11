
# import functions and packages
library(glmnet)
library(tidyverse)
library(doParallel)
library(doRNG)
library(doSNOW)

source("./simData.R")
source("./fit.fun.simulation.R")
source("./stk.glm_v6.R")
source("./random.stk_v6.R")
source("../codes.realdata/ipflasso.R")
source("./baglasso.R")

# data generating mechniasm parameters
# distribution of clinical variables by blocks
n.clin <- c(3, 3, 2, 2, 5)
# distribution of molecular variables by blocks, the 50-th block containing rest variables
n.gene <- rep(10, 49)

# coefficients in original simulation design: De Bin, Riccardo, et al. "Combining clinical and molecular data in regression prediction models: insights from a simulation study." Briefings in Bioinformatics 21.6 (2020): 1904-1919.
beta.clin <- c(3, -3, 0, 3, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0)
beta.gene <- c(3, -3, 1, -1, 1, rep(0, 5), 2, -2, 2, rep(0, 7), 
               rep(0, 10), -3, -2, -1, rep(0, 7),
               rep(-1, 5), rep(0, 5), -2, rep(0, 9), 
               rep(1, 11), rep(0, (10000 - 71)))

# use the absolute coefficients in the present study
beta.clin <- abs(beta.clin)
beta.gene <- abs(beta.gene)
beta <- c(beta.clin, beta.gene)

#-------------------------------------------------------------------------------
# calculate the residual variance for each scenario
para <- c()
set.seed(123)

for(r2 in c(0.3, 0.5, 0.7)){
    for(r in c(0, 0.5, 0.8)){
      
      rho.c <- rep(r, 5)
      rho.b <- rep(r, 5)
      rho.g <- rep(r, 49)
      
      d <- generaClinMol(n.obs = 10000, tot.genes = 1000, n.groups = 49, 
                         n.clin = n.clin, n.gene = n.gene, mean.n.gene = 15, 
                         mu.g = 6, sigma.g = 0.65, mu.c = 1, sigma.c = 0.5, 
                         rho.c = rho.c, rho.b = rho.b, rho.g = rho.g, 
                         phi = 0.1, nu = 10, tau = 20)
      
      d <- cbind(d$clin, d$gene)
      sigma <- cov(d)
      
      # bb <- beta[[b]][1:1015]
      bb <- c(beta.clin, beta.gene)[1:1015]
      v  <- t(bb) %*% sigma %*% bb * (1 - r2) / r2
      gc()        
      v <- c(r2, r, v)
      names(v) <- c("r2", "r", "var")
      para <- rbind(para, v)
    }
}

save(list = "para", file = "../output/para.RData")
para[ ,"var"] <- round(para[ ,"var"])
d <- as.data.frame(para)
write.csv(d, "../results/para.csv")
load("../output/para.RData")
para <- as.data.frame(para)

#-------------------------------------------------------------------------------

#- get the vars
method = c("basic", "stk", "random.stk")
# method = "basic"

nCores = 25
# sample size, 

load("../output/para.RData")
para <- as.data.frame(para)

#-------------------------------------------------------------------------------
#- parallel runs

# set parallel running
nCores <- pmin(detectCores() - 2, nCores)
cl <- makeCluster(nCores)
registerDoSNOW(cl)

# sample size = {100, 500}
# t.gene = {1000, 10000}
# r2 = {0.3, 0.5, 0.7}
# r = {0, 0.5, 0.8}

n <- 2 * 2 * 3 * 3  * 400
# seeds for each duplicates
rng <- RNGseq(n, 1384760)

for(n.n in 1:2){
  for(n.t.gene in 1:2){
    for(n.r2 in 1:3){
      for(n.r in 1:3){
        
        n.obs = c(100, 500)[n.n]
        t.gene = c(1000, 10000)[n.t.gene]
        r2 = c(0.3, 0.5, 0.7)[n.r2]
        r = c(0, 0.5, 0.8)[n.r]
        
        print(Sys.time())
        print(paste0("nobs = ", n.obs))
        print(paste0("N genes = ", t.gene))
        print(paste0("R2 = ", r2))
        print(paste0("r = ", r))
        
        pb <- txtProgressBar(min = 1, max = 400, style = 3)
        progress <- function(k) setTxtProgressBar(pb, k)
        opts <- list(progress = progress)
        
        out <- foreach(times = 1:400, 
                       .packages = c("MASS", "glmnet"),
                       .options.snow = opts,
                       seed = rng[(n.r2 * n.t.gene * n.n * n.r - 1) * 400 + 1:400]) %dopar% {
                         
                         # set seed
                         rngtools::setRNG(seed)
                         which <- which(para$r2 == r2 & para$r == r)
                         sd_sett = sqrt(para[which, "var"])
                         bb <- beta[1:(t.gene + 15)]
                         
                         # generate data
                         dat <- simData(n.obs, t.gene, r, sd_sett, bb)
                         nam <- paste0(c(n.obs, t.gene, r2, r, times), collapse = "_")
                         
                         # export basic models
                         if ("basic" %in% method) {
                           mod <- fit.basic.sim(X = dat$X, y = dat$y, blocks = dat$blocks)              
                           save(list = "mod", file = paste0("../output/basic_", nam, ".RData"))
                         }
                         
                         # export stacking
                         if ("stk" %in% method) {
                           mod <- fit.stk.sim(X = dat$X, y = dat$y, blocks = dat$blocks)
                           save(list = "mod", file = paste0("../output/stk_", nam, ".RData"))
                         }
                         
                         # export random stacking
                         if ("random.stk" %in% method) {
                           mod <- fit.random.stk.sim(X = dat$X, y = dat$y, blocks = dat$blocks)
                           save(list = "mod", file = paste0("../output/random.stk_", nam, ".RData"))
                         }
                       }
      }
    }
  }
}

stopCluster(cl)

#-------------------------------------------------------------------------------
# Prediction
#-------------------------------------------------------------------------------

# calculate model predictions
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
  indRanStk <- grepl("^Rstk", names(obj))
  indBag <- grepl("^bag", names(obj))
  indBasic <- !(indStk | indRanStk | indBag)
  # indBasic <- !(indStk | indRanStk)
  
  # generate results
  re <- list()
  if (any(indBasic)) {
    reBasic<- mapply(predFun, m = obj[indBasic], 
                     MoreArgs = list(newdata = X.v, type = type), 
                     SIMPLIFY = FALSE)
    re <- c(re, reBasic)
  }
  
  if (any(indBag)) {
    reBag<- mapply(pred.baglasso, obj = obj[indBag],
                   MoreArgs = list(newdata = X.v, type = type),
                   SIMPLIFY = FALSE)
    re <- c(re, unlist(reBag, recursive = FALSE))
  }
  # 
  # prediction of stacking
  if (any(indStk)) {
    reStk <- mapply(predict.stk, 
                    obj = obj[indStk],
                    MoreArgs = list(newdata = X.v, type = type, ifCVMods = ifCVMods),
                    SIMPLIFY = FALSE)
    re <- c(re, unlist(reStk, recursive = FALSE))
  }
  
  if (any(indRanStk)) {
    reRanStk <- mapply(predict.random.stk, 
                       obj = obj[indRanStk],
                       MoreArgs = list(newdata = X.v, type = type, type.BL = NULL),
                       SIMPLIFY = FALSE)
    re <- c(re, unlist(reRanStk, recursive = FALSE))
  }
  return(re)
}


set.seed(123)
mes <- list()
nCores <- 25
#- parallel runs
nCores <- pmin(detectCores() - 2, nCores)
cl <- makeCluster(nCores)
registerDoSNOW(cl)

for(n.n in 1:2){
  for(n.t.gene in 1:2){
    for(n.r2 in 1:3){
      for(n.r in 1:3){
        
        n.obs = c(100, 500)[n.n]
        t.gene = c(1000, 10000)[n.t.gene]
        r2 = c(0.3, 0.5, 0.7)[n.r2]
        r = c(0, 0.5, 0.8)[n.r]
        
        print(Sys.time())
        print(paste0("nobs = ", n.obs))
        print(paste0("N genes = ", t.gene))
        print(paste0("R2 = ", r2))
        print(paste0("r = ", r))
        
        which <- which(para$r2 == r2 & para$r == r)
        sd_sett = sqrt(para[which, "var"])
        # generate large data
        bb <- beta[1:(t.gene + 15)]
        dat.v <- simData(10000, t.gene, r, sd_sett, bb)
        # yy, predicted yy
        yy <- dat.v$X %*% bb
        
        pred <- foreach(times = 1:400, 
                        .packages = c("MASS", "glmnet"),
                        .combine = "c") %dopar% {
                          
                          # generate data
                          nam <- paste0(c(n.obs, t.gene, r2, r, times), collapse = "_")
                          
                          mods <- c()
                          
                          if ("basic" %in% method) {
                            load(paste0("../output/basic_", nam, ".RData"))
                            mods <- c(mods, mod)
                          }
                          
                          if ("stk" %in% method) {
                            load(paste0("../output/stk_", nam, ".RData"))
                            mods <- c(mods, mod)
                          }
                          
                          if ("random.stk" %in% method) {
                            load(paste0("../output/random.stk_", nam, ".RData"))
                            mods <- c(mods, mod)
                          }
                          
                          pred <- predMods(dat.v$X, mods, type = "link", ifCVMods = FALSE)
                        }
        
        meDat <- c()
        for(m in unique(names(pred))) {
          
          ind <- grep(m, names(pred))
          mpred <- as.data.frame(pred[ind])
          bias <- mean((rowMeans(mpred) - yy)^2)
          haty <- rowMeans(mpred)
          variance <- mean(rowMeans(sweep(mpred, 1, haty)^2))
          
          me <- mean(colMeans((mpred - yy)^2))
          mese <- sd(colMeans((mpred - yy)^2)) / sqrt(ncol(mpred) - 1)
          
          re <- data.frame(method = m, me = me, mese = mese, bias = bias, variance = variance)
          meDat <- rbind(meDat, re)
        }
        nam <- paste0(c("p", n.obs, t.gene, r2, r), collapse = "_")
        mes[[nam]] <- meDat
      }
    }
  }
}

stopCluster(cl)
save("mes", file = paste0("../results/sim_mes.RData"))

# visualization prediction performance
load(paste0("../results/sim_mes.RData"))

# generate stacked barplot 
barPlot <- function(obj, type = c("stacked", "percent")) {
  
  # rename
  newNames <- c("naive" = "Lasso", "ipf.las" = "IPFlasso", "baglas.BL" = "Bagging", 
                "stk.comp.lasso" = "Stk(CP, Lasso)", "stk.comp.ridge" = "Stk(CP,Rdg)", "stk.comp.solnp" = "Stk(CP,RLM)",
                "stk.pf.lasso" = "Stk(PF, Lasso)", "stk.pf.ridge" = "Stk(PF,Rdg)", "stk.pf.solnp" = "Stk(PF,RLM)",
                "Rstk.comp.lasso" = "RStk(CP, Lasso)", "Rstk.comp.ridge" = "RStk(CP,Rdg)", "Rstk.comp.solnp" = "RStk(CP,RLM)",
                "Rstk.pf.lasso" = "RStk(PF, Lasso)", "Rstk.pf.ridge" = "RStk(PF,Rdg)", "Rstk.pf.solnp" = "RStk(PF,RLM)")
  
  plotDat <- obj %>%
    mutate(method = newNames[obj$method],
           method = factor(method, levels = method),
           grp = case_when(grepl("^Stk", method) ~ "Stk",
                           grepl("^RStk", method) ~ "RStk",
                           .default = method)) %>%
    pivot_longer(cols = c("bias", "variance"), names_to = "measures") %>%
    mutate(percent = value / me * 100,
           super = case_when(grepl("\\(", method) ~ sub(".*,([^)]*)\\)", "\\1", method),
                             .default = ""))
  
  if (type == "stacked") {
    plotDat$percent <- plotDat$value
  }
  
  accuracy <- ifelse(max(plotDat$me) > 90, 1, 0.1)
  
  # if all Rdg, RLM, or Lasso, remove Rdg and RLM is the parentheses
  if (sum(c("Lasso", "Rdg", "RLM") %in% plotDat$super) <= 1) {
    plotDat$method <- gsub(",[^)]*\\)", ")", plotDat$method)
  }
  
  plotDat$method = factor(plotDat$method, levels = unique(plotDat$method))
  plotDat$grp = factor(plotDat$grp, levels = unique(plotDat$grp))
  
  ylab <- ifelse(type == "stacked", "ME", "Percentage, %")
  
  # Create the plot
  p <- ggplot(plotDat, aes(x = method, y = percent, fill = measures)) +
    geom_bar(stat = "identity", position = "stack") +   # Stacked bar plot
    scale_fill_manual(values = c("bias" = "coral1", "variance" = "darkcyan"),
                      labels = c(expression(Bias^2), expression(Variance))) + 
    geom_text(aes(label = format(percent, digits = 1, nsmall = 1)), 
              position = position_stack(vjust = 0.5),
              angle = 0, size = 3, color = "white") +
    theme_bw() +
    ylab(ylab) + xlab(NULL) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_blank(),
          legend.position = "bottom") +
    scale_y_continuous(labels = scales::number_format(accuracy = accuracy))
  
  return(p)
}

# generate mean and error bar plot
meanSEplot <- function(obj, legend = FALSE) {
  
  newNames <- c("naive" = "Lasso", "ipf.las" = "IPFlasso", "baglas.BL" = "Bagging", 
                "stk.comp.lasso" = "Stk(CP,Lasso)", "stk.comp.ridge" = "Stk(CP,Rdg)", "stk.comp.solnp" = "Stk(CP,RLM)",
                "stk.pf.lasso" = "Stk(PF,Lasso)", "stk.pf.ridge" = "Stk(PF,Rdg)", "stk.pf.solnp" = "Stk(PF,RLM)",
                "Rstk.comp.lasso" = "RStk(CP,Lasso)", "Rstk.comp.ridge" = "RStk(CP,Rdg)", "Rstk.comp.solnp" = "RStk(CP,RLM)",
                "Rstk.pf.lasso" = "RStk(PF,Lasso)", "Rstk.pf.ridge" = "RStk(PF,Rdg)", "Rstk.pf.solnp" = "RStk(PF,RLM)")
  
  plotDat <- obj %>%
    mutate(method = newNames[obj$method],
           grp = case_when(grepl("^Stk", method) ~ "Stk",
                           grepl("^RStk", method) ~ "RStk",
                           .default = method),
           super = case_when(grepl("\\(", method) ~ sub(".*,([^)]*)\\)", "\\1", method),
                             .default = ""))
  
  # if all Rdg, RLM, or Lasso, remove Rdg and RLM is the parentheses
  if (sum(c("Lasso", "Rdg", "RLM") %in% plotDat$super) <= 1) {
    plotDat$method <- gsub(",[^)]*\\)", ")", plotDat$method)
  }
  # rank
  plotDat$method = factor(plotDat$method, levels = unique(plotDat$method))
  plotDat$grp = factor(plotDat$grp, levels = unique(plotDat$grp))
  
  accuracy <- ifelse(max(plotDat$me) > 90, 1, 0.1)
  
  p <- ggplot(plotDat, aes(x = method, y = me))+ 
    geom_pointrange(aes(ymin = me - mese, ymax = me + mese, colour = grp), size = 0.2) +
    theme_bw() + 
    ylab("ME") + xlab(NULL) + 
    theme(
      axis.text.x = element_blank(),
      # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      legend.title = element_blank(),
      legend.position = "bottom") + 
    scale_y_continuous(labels = scales::number_format(accuracy = accuracy)) +
    scale_color_brewer(palette = "Set1")
  # coord_flip()
  
  if (!legend){
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

N <- c(100, 500)
P <- c(1000, 10000)
R <- c(0, 0.5, 0.8)
R2 <- c(0.3, 0.5, 0.7)

selMethod <- c("naive", "ipf.las", "stk.comp.ridge", "Rstk.comp.ridge", "stk.pf.ridge", "Rstk.pf.ridge")

# Get all single plots for plots selection
barPlots <- list()
meanPlots <- list()

for(n in N) {
  for(p in P) {
    for(r in R) {
      for(r2 in R2) {
        
        nam <- paste0(c("p", n, p, r2, r), collapse = "_")
        obj <- mes[[nam]]
        obj <- obj[match(selMethod, obj$method),]
        
        # Define your lines and colors
        text <- paste0(c("N=", "P=", "R=", "R2="), c(n, p, r, r2), ";")
        cols <- rep("gray", 4)
        
        # x = min(obj$me)
        x = nrow(obj) - 1.65
        y = Inf
        # trans r2 to snr
        snr <- format(round(r2 / (1- r2), 2), nsmall = 2)
        
        plotObj <- meanSEplot(obj, legend = TRUE) + 
          annotate("text", x = x, y = y, label = paste0("N = ", n), color = "gray", size = 3, hjust = 0, vjust = 2, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("P = ", p), color = "gray", size = 3, hjust = 0, vjust = 3.5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("SNR = ", snr), color = "gray", size = 3, hjust = 0, vjust = 5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("r = ", r), color = "gray", size = 3, hjust = 0, vjust = 6.5, fontface = "bold")
        
        meanPlots[[nam]] <- plotObj
        barPlots[[nam]]  <- barPlot(obj, type = "stacked")
          
        # # get legend
        # plotPercent <- barPlot(obj, type = "stacked")
        # legend.grob <- ggpubr::get_legend(plotPercent)
        # # remove legend
        # tmp_plotPercent <- plotPercent + theme(legend.position = "none")
        # tmp <- ggpubr::ggarrange(plotlist = list(plotObj, tmp_plotPercent), nrow = 2, ncol = 1, heights = c(1, 1.5))
        # attr(tmp, "legend.grob") <- legend.grob
        # class(tmp) <- c(class(tmp), "simPlot")
        # simPlots[[nam]] <- tmp
      }
    } 
  }
}

save(meanPlots, barPlots, file = "../output/simPlots.RData")

load("../output/simPlots.RData")

# function to export selected plots
selPlot <- function(meanPlots, barPlots, selList, nrow = 1, ncol = length(selList), titleList = list(rep(NULL, length(selList)))) {
  
  sels <- unlist(lapply(selList, function(x) paste0(c("p", x), collapse = "_"))) 
  mPlots <- meanPlots[sels]
  # add title
  for(i in seq_along(mPlots)) {
    if (!is.null(titleList[[i]])) {
      mPlots[[i]] <- mPlots[[i]] + ggtitle(titleList[[i]])
    }
  }
  
  bPlots <- barPlots[sels]
  orderPlot <- list()
  for(i in 1:nrow){
    orderPlot <- c(orderPlot, mPlots[(ncol * (i - 1) + 1) : (i * (ncol))]) 
    orderPlot <- c(orderPlot, bPlots[(ncol * (i - 1) + 1) : (i * (ncol))]) 
  }
  
  # ggpubr::ggarrange(plotlist = plots, legend = "right", nrow = nrow, ncol = ncol, legend.grob = attr(plots[[1]], "legend.grob"))
  ggpubr::ggarrange(plotlist = orderPlot, legend = "bottom", nrow = nrow * 2, ncol = ncol, 
                    heights = rep(c(1, 1.2), times = nrow), common.legend = TRUE, legend.grob = ggpubr::get_legend(bPlots[[1]]))
  
}

# selected scenarios
selList <- list(c(100, 10000, 0.3, 0.5),
                c(100, 10000, 0.5, 0.5),
                c(100, 10000, 0.7, 0.5))
# g1 <- selPlot(meanPlots, selList, nrow = 1, ncol = length(selList))
# g2 <- selPlot(barPlots, selList, nrow = 1, ncol = length(selList))
g <-  selPlot(meanPlots, barPlots, selList, nrow = 1, ncol = 3, titleList = list("A,SNR = 0.43", "B,SNR = 1.00", "C,SNR = 2.33"))
g

# by noise to single ratio -----------------------------------------------------

# output plots by SNR
pdf("../results/snr.pdf", width = 9, height = 4.5)

for(n in N) {
  for(p in P) {
    for(r in R) {
      
      fileName <- paste0(c("N", "P", "R"),  c(n, p, r), collapse = "_")
      meanList <- list()
      percentList <- list()
      
      for(r2 in R2) {
        
        nam <- paste0(c("p", n, p, r2, r), collapse = "_")
        obj <- mes[[nam]]
        obj <- obj[match(selMethod, obj$method),]
        
        # Define your lines and colors
        text <- paste0(c("N=", "P=", "R=", "R2="), c(n, p, r, r2), ";")
        cols <- rep("gray", 4)
        
        # x = min(obj$me)
        x = nrow(obj) - 1.6
        y = Inf
        # trans r2 to snr
        snr <- format(round(r2 / (1- r2), 2), nsmall = 2)
        
        plotObj <- meanSEplot(obj, legend = TRUE) + 
          annotate("text", x = x, y = y, label = paste0("N = ", n), color = "gray", size = 3, hjust = 0, vjust = 2, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("P = ", p), color = "gray", size = 3, hjust = 0, vjust = 3.5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("SNR = ", snr), color = "black", size = 3, hjust = 0, vjust = 5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("r = ", r), color = "gray", size = 3, hjust = 0, vjust = 6.5, fontface = "bold")
        # coord_cartesian(clip = "off")
        
        # output all gray annotations 
        meanList <- c(meanList, list(plotObj))
        plotPercent <- barPlot(obj, type = "percent")
        percentList <- c(percentList, list(plotPercent))
        
      }
      
      out1 <- ggpubr::ggarrange(plotlist = meanList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out2 <- ggpubr::ggarrange(plotlist = percentList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out  <- ggpubr::ggarrange(plotlist = list(out1, out2), nrow = 2, heights = c(1, 1.4), common.legend = TRUE)
      
      # ggsave(file = paste0("../../../results/", fileName, ".tiff"),  plot = out, width = 8, height = 5, 
      #        units = "in", dpi = 600, compression = "lzw")
      print(out)
    } 
  }
}

dev.off()

# output figures used in manuscript

load("../output/simPlots.RData")
selList <- list(c(100, 10000, 0.3, 0.5),
                c(100, 10000, 0.5, 0.5),
                c(100, 10000, 0.7, 0.5))
g <-  selPlot(meanPlots, barPlots, selList, nrow = 1, ncol = 3, titleList = list("A,SNR = 0.43", "B,SNR = 1.00", "C,SNR = 2.33"))

ggsave(file = "../results/snrInMain.tiff",  plot = g, width = (3 * 2.5), height = 4.5, units = "in", dpi = 300, compression = "lzw")

# to appendix
selList <- list(c(100, 10000, 0.3, 0),
                c(100, 10000, 0.5, 0),
                c(100, 10000, 0.7, 0),
                c(100, 10000, 0.3, 0.8),
                c(100, 10000, 0.5, 0.8),
                c(100, 10000, 0.7, 0.8))
g <-  selPlot(meanPlots, barPlots, selList, nrow = 2, ncol = 3, 
              titleList = list("A,SNR = 0.43, r = 0", "B,SNR = 1.00, r = 0", "C,SNR = 2.33, r = 0", 
                               "D,SNR = 0.43, r = 0.8", "E,SNR = 1.00, r = 0.8", "F,SNR = 2.33, r = 0.8"))

ggsave(file = "../results/snrInAppendix.tiff",  plot = g, width = (3 * 2.5), height = 9, units = "in", dpi = 300, compression = "lzw")


# by correlations r

selMethod <- c("naive", "ipf.las", "stk.comp.ridge", "Rstk.comp.ridge", "stk.pf.ridge", "Rstk.pf.ridge")
pdf("../results/r.pdf", width = 9, height = 4.5) #

for(n in N) {
  for(p in P) {
    for(r2 in R2) {
      
      fileName <- paste0(c("N", "P", "R2"),  c(n, p, r2), collapse = "_") #
      meanList <- list()
      percentList <- list()
      for(r in R) {
        
        nam <- paste0(c("p", n, p, r2, r), collapse = "_")
        obj <- mes[[nam]]
        obj <- obj[match(selMethod, obj$method),]
        # Define your lines and colors
        text <- paste0(c("N=", "P=", "R=", "R2="), c(n, p, r, r2), ";")
        cols <- rep("gray", 4)
        
        x = length(selMethod) - 1.6
        # trans r2 to snr
        snr <- format(round(r2 / (1- r2), 2), nsmall = 2)
        
        plotObj <- meanSEplot(obj) + 
          annotate("text", x = x, y = y, label = paste0("N = ", n), color = "gray", size = 3, hjust = 0, vjust = 2, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("P = ", p), color = "gray", size = 3, hjust = 0, vjust = 3.5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("SNR = ", snr), color = "gray", size = 3, hjust = 0, vjust = 5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("r = ", r), color = "black", size = 3, hjust = 0, vjust = 6.5, fontface = "bold")
        # coord_cartesian(clip = "off")
        
        plotPercent <- barPlot(obj, type = "stacked") 
        meanList <- c(meanList, list(plotObj))
        percentList <- c(percentList, list(plotPercent))
      }
      
      out1 <- ggpubr::ggarrange(plotlist = meanList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out2 <- ggpubr::ggarrange(plotlist = percentList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out  <- ggpubr::ggarrange(plotlist = list(out1, out2), nrow = 2, heights = c(1, 1.4), common.legend = FALSE)
      
      print(out)
    } 
  }
}

dev.off()

# output figures used in main and appendix text

load("../output/simPlots.RData")
selList <- list(c(100, 10000, 0.5, 0),
                c(100, 10000, 0.5, 0.5),
                c(100, 10000, 0.5, 0.8))

g <-  selPlot(meanPlots, barPlots, selList, nrow = 1, ncol = 3, titleList = list("A,r = 0", "B,r = 0.5", "C,r = 0.8"))
ggsave(file = "../results/rInMain.tiff",  plot = g, width = 3 * 2.5, height = 4.5, units = "in", dpi = 600, compression = "lzw")

# to appendix
# selList <- list(c(100, 10000, 0.3, 0),
#                 c(100, 10000, 0.3, 0.5),
#                 c(100, 10000, 0.3, 0.8),
#                 c(100, 10000, 0.7, 0.0),
#                 c(100, 10000, 0.7, 0.5),
#                 c(100, 10000, 0.7, 0.8))
# 
# g <-  selPlot(meanPlots, barPlots, selList, nrow = 2, ncol = 3,
#               list("A,SNR = 0.43, r = 0", "B,SNR = 0.43, r = 0.5", "C,SNR = 0.43, r = 0.8", 
#                    "D,SNR = 2.33, r = 0", "E,SNR = 2.33, r = 0.5", "F,SNR = 2.33, r = 0.8"))
# ggsave(file = "../results/rInAppendix.tiff",  plot = g, width = 3 * 2.5, height = 9, units = "in", dpi = 600, compression = "lzw")


# by feature size P

selMethod <- c("naive", "ipf.las", "stk.comp.ridge", "Rstk.comp.ridge", "stk.pf.ridge", "Rstk.pf.ridge")

pdf("../../../results/P.pdf", width = 7, height = 4.5) #

for(n in N) {
  for(r2 in R2) {
    for(r in R) {
      fileName <- paste0(c("N", "R2", "R"),  c(n, r2, r), collapse = "_") #
      meanList <- list()
      percentList <- list()
      for(p in P) {
        
        nam <- paste0(c("p", n, p, r2, r), collapse = "_")
        obj <- mes[[nam]]
        obj <- obj[match(selMethod, obj$method),]
        
        # Define your lines and colors
        text <- paste0(c("N=", "P=", "R=", "R2="), c(n, p, r, r2), ";")
        cols <- rep("gray", 4)
        
        x = length(selMethod) - 1.4
        # trans r2 to snr
        snr <- format(round(r2 / (1- r2), 2), nsmall = 2)
        
        plotObj <- meanSEplot(obj) + 
          annotate("text", x = x, y = y, label = paste0("N = ", n), color = "gray", size = 3, hjust = 0, vjust = 2, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("P = ", p), color = "black", size = 3, hjust = 0, vjust = 3.5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("SNR = ", snr), color = "gray", size = 3, hjust = 0, vjust = 5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("r = ", r), color = "gray", size = 3, hjust = 0, vjust = 6.5, fontface = "bold")
        # coord_cartesian(clip = "off")
        
        plotPercent <- barPlot(obj, type = "stacked") 
        meanList <- c(meanList, list(plotObj))
        percentList <- c(percentList, list(plotPercent))
      }
      
      out1 <- ggpubr::ggarrange(plotlist = meanList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out2 <- ggpubr::ggarrange(plotlist = percentList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out  <- ggpubr::ggarrange(plotlist = list(out1, out2), nrow = 2, heights = c(1, 1.4), common.legend = FALSE)
      
      # ggsave(file = paste0("../../../results/", fileName, ".tiff"),  plot = out, width = 5.5, height = 5, 
      #        units = "in", dpi = 600, compression = "lzw")
      print(out)
    } 
  }
}

dev.off()

# output figures used in main and appendix text

load("../output/simPlots.RData")
selList <- list(c(100, 1000, 0.5, 0),
                c(100, 10000, 0.5, 0),
                c(100, 1000, 0.5, 0.5),
                c(100, 10000, 0.5, 0.5))
g <-  selPlot(meanPlots, barPlots, selList, nrow = 1, ncol = 4, 
              titleList = list("A,P = 1000, r = 0", "B,P = 10,000, r = 0", 
                               "C,P = 1000, r = 0.5", "D,P = 10,000, r = 0.5"))
ggsave(file = "../results/pInMain.tiff",  plot = g, width = length(selList) * 2.5, height = 5, units = "in", dpi = 300, compression = "lzw")

# load("../../../output/simPlots_v1.RData")
# selList <- list(c(100, 1000, 0.5, 0.5),
#                 c(100, 10000, 0.5, 0.5))
# g <- selPlot(simPlots, selList, nrow = 1, ncol = length(selList))
# ggsave(file = "../../../results/pInMain.tiff",  plot = g, width = length(selList) * 2.5, height = 5, units = "in", dpi = 600, compression = "lzw")

# to appendix
# selList <- list(c(100, 1000, 0.3, 0),
#                 c(100, 10000, 0.3, 0),
#                 c(100, 1000, 0.3, 0),
#                 c(100, 10000, 0.7, 0.8),
#                 c(100, 1000, 0.7, 0.8),
#                 c(100, 10000, 0.7, 0.8))
# g <-  selPlot(meanPlots, barPlots, selList, nrow = 2, ncol = 3)
# ggsave(file = "../../../results/pInAppendix.tiff",  plot = g, width = 3 * 2.5, height = 10, units = "in", dpi = 300, compression = "lzw")


# by sample size N

selMethod <- c("naive", "ipf.las", "stk.comp.ridge", "Rstk.comp.ridge", "stk.pf.ridge",  "Rstk.pf.ridge")
pdf("../results/N.pdf", width = 6, height = 4.5) #

for(r2 in R2) {
  for(r in R) {
    for(p in P) {
      fileName <- paste0(c("P", "R2", "R"),  c(p, r2, r), collapse = "_") #
      meanList <- list()
      percentList <- list()
      
      for(n in N) {
        
        nam <- paste0(c("p", n, p, r2, r), collapse = "_")
        obj <- mes[[nam]]
        obj <- obj[match(selMethod, obj$method),]
        
        # Define your lines and colors
        text <- paste0(c("N=", "P=", "R=", "R2="), c(n, p, r, r2), ";")
        cols <- rep("gray", 4)
        
        x = length(selMethod) - 1.4
        # trans r2 to snr
        snr <- format(round(r2 / (1- r2), 2), nsmall = 2)
        
        plotObj <- meanSEplot(obj) + 
          annotate("text", x = x, y = y, label = paste0("N = ", n), color = "black", size = 3, hjust = 0, vjust = 2, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("P = ", p), color = "gray", size = 3, hjust = 0, vjust = 3.5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("SNR = ", snr), color = "gray", size = 3, hjust = 0, vjust = 5, fontface = "bold") +
          annotate("text", x = x, y = y, label = paste0("r = ", r), color = "gray", size = 3, hjust = 0, vjust = 6.5, fontface = "bold")
        # coord_cartesian(clip = "off")
        
        plotPercent <- barPlot(obj, type = "stacked") 
        meanList <- c(meanList, list(plotObj))
        percentList <- c(percentList, list(plotPercent))
      }
      out1 <- ggpubr::ggarrange(plotlist = meanList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out2 <- ggpubr::ggarrange(plotlist = percentList, nrow = 1, ncol = length(meanList), common.legend = TRUE, legend = "right")
      out  <- ggpubr::ggarrange(plotlist = list(out1, out2), nrow = 2, heights = c(1, 1.4), common.legend = FALSE)
      
      ggsave(file = paste0("../../../results/", fileName, ".tiff"),  plot = out, width = 5.5, height = 5, 
             units = "in", dpi = 600, compression = "lzw")
      print(out)
    } 
  }
}

dev.off()

# output to manuscript

load("../output/simPlots.RData")
selList <- list(c(100, 10000, 0.5, 0.5),
                c(500, 10000, 0.5, 0.5))
g <-  selPlot(meanPlots, barPlots, selList, nrow = 1, ncol = 2, titleList = list("A,N = 100", "B,N = 500"))
ggsave(file = "../results/nInMain.tiff",  plot = g, width = length(selList) * 3, height = 5, units = "in", dpi = 300, compression = "lzw")

# in appendix
# load("../../../output/simPlots_v4.RData")
# selList <- list(c(100, 10000, 0.5, 0),
#                 c(500, 10000, 0.5, 0),
#                 c(100, 10000, 0.5, 0.8),
#                 c(500, 10000, 0.5, 0.8))
# g <-  selPlot(meanPlots, barPlots, selList, nrow = 1, ncol = 4)
# 
# ggsave(file = "../../../results/nInAppendix.tiff",  plot = g, width = 4 * 2.5, height = 5, units = "in", dpi = 300, compression = "lzw")





















