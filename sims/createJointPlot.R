# ******************************************************************************
#                                 JOINT PLOTS
# ******************************************************************************

# ------------------------------------------------------------------------------
# PREAMBLE
# ------------------------------------------------------------------------------
# load packages
pkgs <- c("dplyr", "ggplot2", "reshape2", "ggpubr")
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
# LOADING DATA
# ------------------------------------------------------------------------------
setwd(paste0(wdir, "res/subenv"))
files <- c("scores_interMean.RData", "scores_interTypes.RData", 
           "scores_k.RData", "scores_nenv.RData", "scores_tau.RData")
for(file in files){
  load(file = file)
  fileName <- strsplit(file, ".RData")[[1]]
  assign(fileName, scoresAll)
}
rm(scoresAll)

# ------------------------------------------------------------------------------
# PLOTTING
# ------------------------------------------------------------------------------
plotsReady <- c("interMean", "interType", "nenv", "tau")

# INTERVENTION MEAN
# ------------------------------------------------------------------------------
if("interMean" %in% plotsReady) {
  methods <- rownames(scores_interMean[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_interMean)))
  rownames(df) <- methods
  colnames(df) <- names(scores_interMean)
  for(interMean in names(scores_interMean)) {
    for(method in methods) {
      df[method, interMean] <- scores_interMean[[interMean]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_interMean <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                       shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("INTERVENTION MEAN") +
    ylab("SUCCESS PROBABILITY")
}

# INTERVENTION STRENGTH
# ------------------------------------------------------------------------------
if("interStrength" %in% plotsReady) {
  methods <- rownames(scores_interStrength[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_interStrengths)))
  rownames(df) <- methods
  colnames(df) <- names(scores_interStrength)
  for(interStrength in names(scores_interStrength)) {
    for(method in methods) {
      df[method, interStrength] <- scores_interStrength[[interStrength]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_interStrength <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                           shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("INTERVENTION STRENGTH") +
    ylab("SUCCESS PROBABILITY")
}

# INTERVENTION TYPE
# ------------------------------------------------------------------------------
if("interType" %in% plotsReady) {
  methods <- rownames(scores_interTypes[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_interTypes)))
  rownames(df) <- methods
  colnames(df) <- names(scores_interTypes)
  for(interType in names(scores_interTypes)) {
    for(method in methods) {
      df[method, interType] <- scores_interTypes[[interType]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_interType <- ggplot(df_melted, aes(x = variable, y = value, group = method, 
                                       colour = method, shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    scale_x_discrete(labels = c("do", "soft", "simult. noise")) +
    xlab("INTERVENTION TYPE") +
    ylab("SUCCESS PROBABILITY")
}

# k
# ------------------------------------------------------------------------------
if("k" %in% plotsReady) {
  methods <- rownames(scores_k[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_k)))
  rownames(df) <- methods
  colnames(df) <- names(scores_k)
  for(k in names(scores_k)) {
    for(method in methods) {
      df[method, k] <- scores_k[[k]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_k <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                               shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("k") +
    ylab("SUCCESS PROBABILITY")
}

# n
# ------------------------------------------------------------------------------
if("n" %in% plotsReady) {
  methods <- rownames(scores_n[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_n)))
  rownames(df) <- methods
  colnames(df) <- names(scores_n)
  for(n in names(scores_n)) {
    for(method in methods) {
      df[method, n] <- scores_n[[n]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_n <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                               shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("n") +
    ylab("SUCCESS PROBABILITY")
}

# nenv
# ------------------------------------------------------------------------------
if("nenv" %in% plotsReady) {
  methods <- rownames(scores_nenv[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_nenv)))
  rownames(df) <- methods
  colnames(df) <- names(scores_nenv)
  for(nenv in names(scores_nenv)) {
    for(method in methods) {
      df[method, nenv] <- scores_nenv[[nenv]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_nenv <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                  shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("# ENVIRONMENTS") +
    ylab("SUCCESS PROBABILITY")
}

# nsubenv
# ------------------------------------------------------------------------------
if("nsubenv" %in% plotsReady) {
  methods <- rownames(scores_subnenv[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_nsubenv)))
  rownames(df) <- methods
  colnames(df) <- names(scores_subnenv)
  for(subenv in names(scores_subnenv)) {
    for(method in methods) {
      df[method, subenv] <- scores_subnenv[[subenv]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_subenv <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                    shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("# SUBENVIRONMENTS") +
    ylab("SUCCESS PROBABILITY")
}

# p
# ------------------------------------------------------------------------------
if("p" %in% plotsReady) {
  methods <- rownames(scores_p[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_p)))
  rownames(df) <- methods
  colnames(df) <- names(scores_p)
  for(p in names(scores_p)) {
    for(method in methods) {
      df[method, p] <- scores_p[[p]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_p <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                               shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("p") +
    ylab("SUCCESS PROBABILITY")
}

# tau
# ------------------------------------------------------------------------------
if("tau" %in% plotsReady) {
  methods <- rownames(scores_tau[[1]][["FWER"]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_tau)))
  rownames(df) <- methods
  colnames(df) <- names(scores_tau)
  for(tau in names(scores_tau)) {
    for(method in methods) {
      df[method, tau] <- scores_tau[[tau]][["successProbability"]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_tau <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                 shape = method)) +
    geom_point(size = 2) +
    geom_line(size = 0.3) +
    expand_limits(y = 1) +
    theme(panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "black", fill =NA, size = 1), 
          legend.key=element_blank(), 
          legend.position = "right") +
    guides(color = guide_legend(title = 'Legend')) + 
    scale_colour_manual(name = 'Legend', 
                        labels = rowOrder, 
                        values = c('yellow2', 'orange', 'mediumpurple1',  'green3', 'blue', 'cyan1', 'red')) + 
                                              scale_shape_manual(name = 'Legend', 
                                                                 labels = rowOrder, 
                                                                 values = c(1, 2, 3, 4, 5, 6, 7)) +
    xlab("TAU") +
    ylab("SUCCESS PROBABILITY")
}

# JOINT PLOT
# ------------------------------------------------------------------------------
p_joint <- ggarrange(p_tau, p_nenv, p_interMean, p_interType, p_tau, p_nenv, p_interMean, p_interType, p_interType,
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom")
p_joint
setwd(paste0(wdir, "fig/subenv"))
ggexport(p_joint, filename = "jointPlot_RICP.pdf")


