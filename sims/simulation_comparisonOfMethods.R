# ******************************************************************************
#                     SIMULATION: VANILLA COMPARISON OF METHODS
# ******************************************************************************

# ------------------------------------------------------------------------------
# PREAMBLE
# ------------------------------------------------------------------------------
# load packages
pkgs <- c("dplyr", "lme4", "nlme", "InvariantCausalPrediction", "nonlinearICP", 
          "pcalg", "ggplot2", "reshape2", "doParallel")
for(pkg in pkgs){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# setting the correct working directory
if("Linux" %in% Sys.info()) {
  wdir <- "~/RICP/"
} else {
  wdir <- "/Users/florianschwarb/Desktop/Master-Thesis/Code/RICP/"
}
setwd(paste0(wdir, "src"))

source("RICP.R") 
source("getpval.R")
source("getpvalwsubenvs.R")
source("lmeFit.R")
source("simDAG.R")
source("simDAGwsubenvs.R")
source("runSimRICP.R")

# ------------------------------------------------------------------------------
# SIMULATION: VANILLA COMPARISON OF METHODS
# ------------------------------------------------------------------------------
# parameters
nsim <- 100

# initializing the cluster
if("Linux" %in% Sys.info()) {
  cl <- makeCluster(50)
} else {
  nCores <- detectCores()
  cl <- makeCluster(nCores - 1)
}
clusterEvalQ(cl, c(library(dplyr), library(lme4), library(nlme),
                   library(InvariantCausalPrediction), library(nonlinearICP),
                   library(pcalg)))
clusterExport(cl, c("nsim"), envir = environment())
clusterExport(cl, c("RICP", "getpvalwsubenvs", "lmeFit", "simDAGwsubenvs", "runSimRICP"),
              envir = environment())

# running simulation in parallel
res <- parLapply(cl, 1:nsim, function(sim) {
  runSimRICP(p = 5, k = 2, nenv = 10, renv = c(80, 100), rBeta = c(-5, 5), tau = 1, 
             alpha = 0.05, interType = "do", interMean = 2, interStrength = 5, 
             subenvs = T, nsubenvs = 30, 
             methods = c("random", "pooled regression", "GES", "LinGAM", "ICP", 
             "nonlinearICP", "RICP"))
  })

# shutting down cluster
stopCluster(cl)

# compute average over all simulation runs
scores <- list()
metrics <- c("FWER", "avgJaccard", "successProbability")
methods <- names(res[[1]]$acceptedSets)
for(metric in metrics) {
  board <- matrix(NA, nrow = length(methods), ncol = nsim + 1)
  colnames(board) <- c(1:nsim, "avg")
  rownames(board) <- methods
  for(method in methods) {
    for(sim in 1:nsim) {
      board[method, sim] <- res[[sim]][[metric]][[method]]
    }
  }
  board[, "avg"] <- sapply(1:length(methods), function(i) {mean(board[i, 1:nsim])})
  scores[[metric]] <- board
}

# saving as .RData-file
setwd(paste0(wdir, "res"))
save(scores, file = "scores_comparison-of-methods.RData")

# PLOTS
# ------------------------------------------------------------------------------
ylabNames <- list(FWER = "FWER", avgJaccard = "AVERAGE JACCARD SIMILARITY", 
                  successProbability = "SUCCESS PROBABILITY")
metrics <- c("FWER", "avgJaccard", "successProbability")
for(metric in metrics) {
  df <- scores[[metric]][, -ncol(scores[[metric]]), drop = FALSE] %>%
    as.data.frame()
  df$method <- rownames(df)
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df <- df %>% slice(match(rowOrder, method))
  df_melted <- melt(df, id = "method")
  df_melted[, "variable"] <- as.double(df_melted[, "variable"])
  iter <- 1
  for(method in rownames(df)) {
    rowInd <- which(df_melted[, "method"] == method)
    xvals <- sample(seq(0.4, 0.6, length.out = 1000), size = length(rowInd))
    df_melted[rowInd, "sim"] <- iter + xvals
    iter <- iter + 1
  }
  quantiles <- sapply(1:nrow(df), function(i) quantile(df[i, -ncol(df)], probs = c(0, 0.25, 0.5, 0.75, 1)) %>% unlist()) %>%
    suppressWarnings()
  p_compMethods <- ggplot(df_melted, aes(x = sim, y = value, group = variable)) +
    geom_point(color = "darksalmon", shape = 21, size = 3) +
    geom_line(color = "grey", size = 0.1) +
    expand_limits(x = nrow(df) + 1, y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    xlab("") +
    ylab(ylabNames[[metric]])
  iter <- 1
  for(method in rownames(df)) {
    if(max(df[iter, -ncol(df)]) < 0.5){
      if(!(method == "pooled regression")) {
        txtStart <- min(df[iter, -ncol(df)]) + 0.05
      } else{
        txtStart <- min(df[iter, -ncol(df)]) + 0.12
      }
    } else {
      if(!(method == "pooled regression")) {
        txtStart <- max(df[iter, -ncol(df)]) - 0.07
      } else{
        txtStart <- max(df[iter, -ncol(df)]) - 0.15
      }
    }
    txt <- paste0("geom_text(label = '", method, "', x = ", 1 + iter, ", y = ", txtStart,
                  ", size = 5, angle = 90, color = 'darkgrey')")
    p_compMethods <- p_compMethods +   
      geom_segment(size = 0.3, aes_(x = iter + 0.3, y = quantiles[1, iter] , xend = iter + 0.7, yend = quantiles[1, iter])) +
      geom_segment(size = 0.3, aes_(x = iter + 0.3, y = quantiles[2, iter] , xend = iter + 0.7, yend = quantiles[2, iter])) +
      geom_segment(size = 0.3, aes_(x = iter + 0.3, y = quantiles[3, iter] , xend = iter + 0.7, yend = quantiles[3, iter])) +
      geom_segment(size = 0.3, aes_(x = iter + 0.3, y = quantiles[4, iter] , xend = iter + 0.7, yend = quantiles[4, iter])) +
      geom_segment(size = 0.3, aes_(x = iter + 0.3, y = quantiles[5, iter] , xend = iter + 0.7, yend = quantiles[5, iter])) +
      eval(parse(text = txt))
    iter <- iter + 1
  }
  
  setwd(paste0(wdir, "fig"))
  ggsave(paste0("compMethods_", metric, ".pdf"), width = 7, height = 5)
}



