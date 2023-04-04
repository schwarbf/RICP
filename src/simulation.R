# ******************************************************************************
#                                 SIMULATIONS
# ******************************************************************************

# ------------------------------------------------------------------------------
# PREAMBLE
# ------------------------------------------------------------------------------
# load packages
libsre <- .libPaths()
pkgs <- c("dplyr", "lme4", "nlme", "InvariantCausalPrediction", "nonlinearICP", 
          "pcalg", "ggplot2", "reshape2", "doParallel", "tictoc")
for(pkg in pkgs){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

wdir <- "/Users/florianschwarb/Desktop/Master-Thesis/Code/R-ICP"
setwd(wdir)

source("RICP.R") 
source("getpval.R")
source("getpvalwsubenvs.R")
source("lmeFit.R")
source("simDAG.R")
source("simDAGwsubenvs.R")
source("runSimRICP.R")

# ------------------------------------------------------------------------------
# SIMULATION 1: COMPARISON OF METHODS
# ------------------------------------------------------------------------------
# parameters
nsim <- 3

# initializing the cluster
nCores <- detectCores()
cl <- makeCluster(nCores - 1)

# running simulation in parallel
tic()
res <- foreach(sim = 1:nsim) %do% {
  runSimRICP(p = 4, k = 1, nenv = 5, renv = c(80, 100), rBeta = c(-5, 5), tau = 5, 
             alpha = 0.05, interType = "do", interMean = 2, interStrength = 5, 
             subenvs = T, nsubenvs = 30) 
}
toc()

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
scores$successProbability
scores$FWER
scores$avgJaccard

# PLOTTING
# ------------------------------------------------------------------------------



