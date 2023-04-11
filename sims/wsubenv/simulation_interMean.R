# ******************************************************************************
#                          SIMULATION: INCREASING TAU
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
# SIMULATION: INCREASING TAU
# ------------------------------------------------------------------------------
# parameters
interMeans <- c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 50)
nsim <- 50

# initializing the cluster
if("Linux" %in% Sys.info()) {
  cl <- makeCluster(50)
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
scoresAll <- list()
for(interMean in interMeans) {
  # run simulations
  clusterExport(cl, "interMean")
  res <- parLapply(cl, 1:nsim, function(sim) {
    runSimRICP(p = 5, k = 2, nenv = 100, renv = c(80, 100), rBeta = c(-5, 5), tau = 0.5,
               alpha = 0.05, interType = "do", interMean = interMean, interStrength = 10,
               nInter = "multiple", subenvs = F, nsubenvs = 30, test = "lme4",
               methods = c("random", "pooled regression", "GES", "LinGAM", "ICP",
                           "nonlinearICP", "RICP"))
  })
  
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
  scoresAll[[as.character(interMean)]] <- scores
  
  # progress bar
  cat(paste0("*** ", round(100 * which(interMeans == interMean)/length(interMeans)), 
             "% complete: tested ", which(interMeans == interMean), " out of ", length(interMeans), 
             " intervention means \n"))
}

# shutting down cluster
stopCluster(cl)

# saving as .RData-file
setwd(paste0(wdir, "res/wsubenv"))
save(scoresAll, file = "scores_interMean.RData")

# PLOTS
# ------------------------------------------------------------------------------
methods <- rownames(scoresAll[[1]][["FWER"]])
df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(interMeans)))
rownames(df) <- methods
colnames(df) <- names(scoresAll)
for(interMean in names(scoresAll)) {
  for(method in methods) {
    df[method, interMean] <- scoresAll[[interMean]][["successProbability"]][method, "avg"]
  }
}
df$method <- methods
rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
df$method <- factor(df$method, levels = rowOrder)
df_melted <- melt(df, id = "method")
p_interMean <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                               shape = method)) +
  geom_point(size = 4) +
  geom_line(size = 0.3) +
  expand_limits(y = 1) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 1), 
        legend.key=element_blank(), 
        legend.position = "right") +
  guides(color = guide_legend(title = 'intervention mean')) + 
  scale_colour_manual(name = 'intervention mean', 
                      labels = rowOrder, 
                      values = c('yellow3', 'orange', 'mediumpurple1', 'purple4', 
                                          'blue', 'brown', 'red')) + 
  scale_shape_manual(name = 'intervention mean', 
                     labels = rowOrder, 
                     values = c(1, 2, 3, 4, 5, 6, 7)) +
  xlab("INTERVENTION MEAN") +
  ylab("SUCCESS PROBABILITY")

setwd(paste0(wdir, "fig/wsubenv"))
ggsave(paste0("interMean_", metric, ".pdf"), width = 6, height = 5)




