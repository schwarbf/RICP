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
if("florianschwarb" %in% Sys.info()){
  wdir <- "/Users/florianschwarb/Desktop/Master-Thesis/Code/RICP/"
} else if("Linux" %in% Sys.info()) {
  wdir <- "~/RICP/"
} else{
  wdir <- getwd()
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
  cl <- makeCluster(100)
} else {
  nCores <- detectCores()
  cl <- makeCluster(nCores - 1)
}
clusterEvalQ(cl, c(library(dplyr), library(lme4), library(nlme),
                   library(InvariantCausalPrediction), library(nonlinearICP),
                   library(pcalg))) %>% invisible()
clusterExport(cl, c("nsim"), envir = environment())
clusterExport(cl, c("RICP", "getpval", "getpvalwsubenvs", "lmeFit", "simDAG",
                    "simDAGwsubenvs", "runSimRICP"),
              envir = environment())

# running simulation in parallel
res <- parLapply(cl, 1:nsim, function(sim) {
  runSimRICP(p = 5, k = 2, nenv = 30, renv = c(80, 100), rBeta = c(-5, 5), tau = 0.5,
             alpha = 0.05, interType = "do", interMean = 2, interStrength = 5,
             nInter = "multiple", subenvs = T, nsubenvs = 100, test = "LRT-lme4",
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
save(scores, file = "scores_comparisonOfMethods.RData")

# PLOTS
# ------------------------------------------------------------------------------
ylabNames <- list(FWER = "FWER", avgJaccard = "AVERAGE JACCARD SIMILARITY", 
                  successProbability = "SUCCESS PROBABILITY")
metrics <- c("FWER", "avgJaccard", "successProbability")
p_compMethods <- list()
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
    xvals <- sample(seq(0.2, 0.8, length.out = 1000), size = length(rowInd))
    df_melted[rowInd, "simSpread"] <- iter - 0.5 + xvals
    df_melted[rowInd, "sim"] <- iter 
    iter <- iter + 1
  }
  quantiles <- sapply(1:nrow(df), function(i) quantile(df[i, -ncol(df)], probs = c(0, 0.25, 0.5, 0.75, 1)) %>% unlist()) %>%
    suppressWarnings()
  p_compMethods[[metric]] <- ggplot(df_melted, aes(x = sim, y = value, group = sim)) +
    geom_point(aes(x = simSpread, y = value), color = "darksalmon", shape = 4, size = 1.5) +
    geom_boxplot(outlier.size = 0) +
    geom_line(aes(x = simSpread, y = value, group = variable), color = "grey", size = 0.02) +
    expand_limits(x = nrow(df) + 1, y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7), labels = rowOrder) +
    ylab(ylabNames[[metric]])
  if(metric == "FWER") {
    p_compMethods[[metric]] <- p_compMethods[[metric]] + geom_hline(yintercept = 0.05, linetype = 3)
  }
  iter <- 1
  for(method in rownames(df)) {
    # p_compMethods[[metric]] <- p_compMethods[[metric]] +   
      # geom_segment(size = 0.3, aes_(x = iter + 0.1, y = quantiles[1, iter] , xend = iter + 0.9, yend = quantiles[1, iter])) +
      # geom_segment(size = 0.3, aes_(x = iter + 0.1, y = quantiles[2, iter] , xend = iter + 0.9, yend = quantiles[2, iter])) +
      # geom_segment(size = 0.3, aes_(x = iter + 0.1, y = quantiles[3, iter] , xend = iter + 0.9, yend = quantiles[3, iter])) +
      # geom_segment(size = 0.3, aes_(x = iter + 0.1, y = quantiles[4, iter] , xend = iter + 0.9, yend = quantiles[4, iter])) +
      # geom_segment(size = 0.3, aes_(x = iter + 0.1, y = quantiles[5, iter] , xend = iter + 0.9, yend = quantiles[5, iter])) 
    iter <- iter + 1
  }
  
  setwd(paste0(wdir, "fig"))
  ggsave(paste0("compMethods_", metric, ".png"), dpi = "retina", width = 13, height = 10, units = "cm", device = "png")
}


