# ******************************************************************************
#                                 JOINT PLOTS
# ******************************************************************************

# ------------------------------------------------------------------------------
# PREAMBLE
# ------------------------------------------------------------------------------
# load packages
pkgs <- c("dplyr", "ggplot2", "reshape2", "ggpubr", "RColorBrewer", "wesanderson")
for(pkg in pkgs){
  if(!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# setting the correct working directory
if("florianschwarb" %in% Sys.info()){
  wdir <- "/Users/florianschwarb/Desktop/Master-Thesis/Code/RICP/"
} else{
  wdir <- getwd()
}

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------
simSimpleDAG <- function(nenv = 10, renv = c(200, 300), tau = 0.25, 
                         interMean = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10),
                         interStrength = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) {
  # sampling DAG
  Beta <- matrix(0, ncol = 4, nrow = 4)
  Beta[1, ] <- c(0, 1, 4, 0)
  Beta[2, ] <- c(0, 0, 2, 2)
  Beta[3, ] <- c(0, 0, 0, 3)
  colnames(Beta) <- rownames(Beta) <- c("X1", "X2", "Y", "X3")
  
  # sample data within environments 
  n <- ceiling(runif(nenv, min = renv[1], max = renv[2]))
  eps <- matrix(0, ncol = 4, nrow = sum(n))
  X <- matrix(0, ncol = 4, nrow = sum(n))
  colnames(X) <- colnames(Beta)
  B <- list()
  interInd <- list()
  for(env in 1:nenv) {
    if(env != 1) {
      numInter <- sample(1:3, 1)
      interInd[[env]] <- sample(c(1, 2, 4), numInter) 
    }
    interMean <- if(length(interMean) > 1) sample(interMean, size = 1) else interMean
    interStrength <- if(length(interStrength) > 1) sample(interStrength, size = 1) else interStrength
    interSig <- interMean + interStrength*runif(1)
    rowInd <- if(env == 1) {1:n[env]} else {(sum(n[1:(env-1)]) +1):(sum(n[1:(env-1)]) + n[env])}
    # sdEps <- sample(c(0.1, 0.2, 0.5, 1, 2), size = 1)
    sdEps <- 0.2
    eps[rowInd, ] <- matrix(rnorm(4*n[env], sd = sdEps), nrow = n[env], ncol = 4)
    B[[env]] <- Beta
    B[[env]][which(Beta > 0, arr.ind = TRUE)] <- rnorm(length(which(Beta > 0)), sd = tau)
    for(j in 1:4) {
      if(env != 1) {
        X[rowInd, j] <- X[rowInd, 1:j, drop = FALSE] %*% 
          (Beta[1:j, j] + B[[env]][1:j, j]) + eps[rowInd, j] 
      } else{
        X[rowInd, j] <- X[rowInd, 1:j, drop = FALSE] %*% Beta[1:j, j] + eps[rowInd, j] 
      }
      if(env != 1 && j %in% interInd[[env]]) {
        X[rowInd, j] <- interSig
      }
    }
  }
  ExpInd <- sapply(1:nenv, function(i){rep(i, n[i])}) %>% unlist()
  
  return(list(X = X, Beta = Beta, interInd = interInd, ExpInd = ExpInd))
}

# ------------------------------------------------------------------------------
# PLOTTING
# ------------------------------------------------------------------------------
# parameters
interMean <- 0.2
nenv <- 15
tau <- 1

# simulating the DAG
tmp <- simSimpleDAG(nenv = nenv, renv = c(300, 400), tau = tau, interMean = interMean, 
                    interStrength = 0)
data <- tmp$X %>% as.data.frame()
ExpInd <- tmp$ExpInd

# creating the regressions
df <- data
envs <- unique(tmp$ExpInd)
preds <- c("X1", "X2", "X3")
df_pred <- data.frame(matrix(NA, nrow = nenv*1000, ncol = 3))
for(env in envs) {
  indEnv <- which(tmp$ExpInd == env)
  indEnv_pred <- if(env == 1) 1:1000 else ((env-1)*1000 + 1):(env*1000)
  for(pred in preds) {
    predNum <- as.numeric(gsub(".*?([0-9]+).*", "\\1", pred))
    tmpPred <- if(predNum == 3) 4 else predNum
    df_pred[indEnv_pred, "ExpInd"] <- env
    if(!(tmpPred %in% tmp$interInd[[env]])) {
      lm.fit <- lm(as.formula(paste0("Y ~ ", pred)), data = data[indEnv, , drop = FALSE])
      newdata <- data.frame(matrix(NA, ncol = 1, nrow = 1000))
      df_pred[indEnv_pred, pred] <- newdata[, pred] <- seq(min(df[, pred]), max(df[, pred]), length.out = 1000)
      df_pred[indEnv_pred, paste0("Y", pred)] <- predict(lm.fit, newdata = newdata)
    }
  }
}
df$ExpInd <- ExpInd

# Y ~ X1
ExpIndX1 <- unique(df[intersect(which(df[, "X1"] != interMean), which(df[, "ExpInd"] != 1)), "ExpInd"])
rowIndX1 <- which(df[, "ExpInd"] %in% ExpIndX1)
rowIndX1_pred <- which(df_pred[, "ExpInd"] %in% ExpIndX1)
p_X1 <- ggplot(df[rowIndX1, ], aes(x = X1, y = Y, group = ExpInd, color = factor(ExpInd))) +
  geom_point(size = 0.5, shape = 4) +
  geom_line(data = df_pred[rowIndX1_pred, , drop = FALSE], aes(x = X1, y = YX1, group = ExpInd, colour = factor(ExpInd)), linewidth = 0.3) +
  geom_point(data = df[which(df$ExpInd == 1),  , drop = FALSE], aes(x = X1, y =Y), color = "royalblue", size = 0.5, shape = 21, alpha = 0.5) +
  geom_line(data = df_pred[which(df_pred$ExpInd == 1),  , drop = FALSE], aes(x = X1, y =YX1), color = "royalblue", linewidth = 1) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 1), 
        legend.key=element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = wes_palette("FantasticFox1", length(ExpIndX1), type = "continuous")) +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 2, label = "observational data", color = "royalblue")
p_X1
setwd(paste0(wdir, "fig/simpleExample"))
ggsave("simpleExample_X1.pdf", width = 3, height = 2.5)

# Y ~ X2
ExpIndX2 <- unique(df[intersect(which(df[, "X2"] != interMean), which(df[, "ExpInd"] != 1)), "ExpInd"])
rowIndX2<- which(df[, "ExpInd"] %in% ExpIndX2)
rowIndX2_pred <- which(df_pred[, "ExpInd"] %in% ExpIndX2)
p_X2 <- ggplot(df[rowIndX2, ], aes(x = X2, y = Y, group = ExpInd, color = factor(ExpInd))) +
  geom_point(size = 0.5, shape = 4) +
  geom_line(data = df_pred[rowIndX2_pred, , drop = FALSE], aes(x = X2, y = YX2, group = ExpInd, colour = factor(ExpInd)), linewidth = 0.3) +
  geom_point(data = df[which(df$ExpInd == 1),  , drop = FALSE], aes(x = X2, y =Y), color = "royalblue", size = 0.5, shape = 21, alpha = 0.5) +
  geom_line(data = df_pred[which(df_pred$ExpInd == 1),  , drop = FALSE], aes(x = X2, y =YX2), color = "royalblue", linewidth = 1) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 1), 
        legend.key=element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = wes_palette("FantasticFox1", length(ExpIndX2), type = "continuous")) +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 2, label = "observational data", color = "royalblue")
p_X2
ggsave("simpleExample_X2.pdf", width = 3, height = 2.5)

# Y ~ X3
ExpIndX3 <- unique(df[intersect(which(df[, "X3"] != interMean), which(df[, "ExpInd"] != 1)), "ExpInd"])
rowIndX3 <- which(df[, "ExpInd"] %in% ExpIndX3)
rowIndX3_pred <- which(df_pred[, "ExpInd"] %in% ExpIndX3)
p_X3 <- ggplot(df[rowIndX3, ], aes(x = X3, y = Y, group = ExpInd, color = factor(ExpInd))) +
  geom_point(size = 0.5, shape = 4) +
  geom_line(data = df_pred[rowIndX3_pred, , drop = FALSE], aes(x = X3, y = YX3, group = ExpInd, colour = factor(ExpInd)), linewidth = 0.3) +
  geom_point(data = df[which(df$ExpInd == 1),  , drop = FALSE], aes(x = X3, y =Y), color = "royalblue", size = 0.5, shape = 21, alpha = 0.5) +
  geom_line(data = df_pred[which(df_pred$ExpInd == 1),  , drop = FALSE], aes(x = X3, y =YX3), color = "royalblue", linewidth = 1) +
  theme(panel.background = element_blank(), 
        panel.grid = element_blank(), 
        panel.border = element_rect(colour = "black", fill =NA, size = 1), 
        legend.key=element_blank(), 
        legend.position = "none") +
  scale_color_manual(values = wes_palette("FantasticFox1", length(ExpIndX3), type = "continuous")) +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.05, vjust = 2, label = "observational data", color = "royalblue")
p_X3
ggsave("simpleExample_X3.pdf", width = 3, height = 2.5)

# joint plot
p_joint <- ggarrange(p_X1, p_X2, p_X3, ncol = 3, nrow = 1)
p_joint
setwd(paste0(wdir, "fig/simpleExample"))
ggexport(p_joint, filename = paste0("jointPlot_RICP_simpleExample.png"), width = 350, height = 125)

