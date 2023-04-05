# ******************************************************************************
#                                     TEST
# ******************************************************************************

# ------------------------------------------------------------------------------
# PREAMBLE
# ------------------------------------------------------------------------------
# load packages
pkgs <- c("dplyr", "lme4", "nlme", "InvariantCausalPrediction", "nonlinearICP", 
          "pcalg", "ggplot2", "reshape2", "doParallel", "tictoc")
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
# TEST CLUSTER
# ------------------------------------------------------------------------------
# initializing the cluster
nCores <- detectCores()
print(nCores)
cl <- makeCluster(nCores)

ns <- 200
res <- foreach(n = 1:ns) %do% {
  sqrt(n)
}
print(res)

# shutting down cluster
stopCluster(cl)