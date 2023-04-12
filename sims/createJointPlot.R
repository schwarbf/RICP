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
setwd(paste0(wdir, "res/wsubenv/multipleInter"))
files <- c("scores_p.RData", "scores_k.RData", "scores_n.RData", "scores_nenv.RData", 
           "scores_interMean.RData", "scores_interTypes.RData", "scores_k.RData", 
           "scores_nInter.RData", "scores_interStrength.RData", "scores_tau.RData")
for(file in files){
  load(file = file)
  fileName <- strsplit(file, ".RData")[[1]]
  assign(fileName, scoresAll)
}
rm(scoresAll)

# ------------------------------------------------------------------------------
# PLOTTING
# ------------------------------------------------------------------------------
plotsReady <- c("p", "k", "n", "tau", "interType", "interMean", "interStrength", "nInter", "nenv")
metric <- "avgJaccard" # successProbability, avgJaccard, FWER

metricName <- if(metric == "successProbability") {"SUCCESS PROBABILITY"} else if(metric == "avgJaccard") {"JACCARD SIMILARITY"} else {"FWER"}

# INTERVENTION MEAN
# ------------------------------------------------------------------------------
if("interMean" %in% plotsReady) {
  methods <- rownames(scores_interMean[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_interMean)))
  rownames(df) <- methods
  colnames(df) <- names(scores_interMean)
  for(interMean in names(scores_interMean)) {
    for(method in methods) {
      df[method, interMean] <- scores_interMean[[interMean]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_interMean <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                       shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    xlab("INTER. MEAN") +
    ylab(metricName)
  if(metric == "FWER") {
    p_interMean <- p_interMean + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# INTERVENTION STRENGTH
# ------------------------------------------------------------------------------
if("interStrength" %in% plotsReady) {
  methods <- rownames(scores_interStrength[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_interStrength)))
  rownames(df) <- methods
  colnames(df) <- names(scores_interStrength)
  for(interStrength in names(scores_interStrength)) {
    for(method in methods) {
      df[method, interStrength] <- scores_interStrength[[interStrength]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_interStrength <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                           shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    xlab("INTER. STRENGTH") +
    ylab(metricName)
  if(metric == "FWER") {
    p_interStrength <- p_interStrength + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# INTERVENTION TYPE
# ------------------------------------------------------------------------------
if("interType" %in% plotsReady) {
  methods <- rownames(scores_interTypes[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_interTypes)))
  rownames(df) <- methods
  colnames(df) <- names(scores_interTypes)
  for(interType in names(scores_interTypes)) {
    for(method in methods) {
      df[method, interType] <- scores_interTypes[[interType]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_interType <- ggplot(df_melted, aes(x = variable, y = value, group = method, 
                                       colour = method, shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    xlab("INTER. TYPE") +
    ylab(metricName)
  if(metric == "FWER") {
    p_interType <- p_interType + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# NUMBER OF INTERVENTIONS
# ------------------------------------------------------------------------------
if("nInter" %in% plotsReady) {
  methods <- rownames(scores_nInter[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_nInter)))
  rownames(df) <- methods
  colnames(df) <- names(scores_nInter)
  for(numInter in names(scores_nInter)) {
    for(method in methods) {
      df[method, numInter] <- scores_nInter[[numInter]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_nInter <- ggplot(df_melted, aes(x = variable, y = value, group = method, 
                                    colour = method, shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    xlab("# INTERVENTIONS/ENV.") +
    ylab(metricName)
  if(metric == "FWER") {
    p_interType <- p_interType + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# k
# ------------------------------------------------------------------------------
if("k" %in% plotsReady) {
  methods <- rownames(scores_k[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_k)))
  rownames(df) <- methods
  colnames(df) <- names(scores_k)
  for(k in names(scores_k)) {
    for(method in methods) {
      df[method, k] <- scores_k[[k]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_k <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                               shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    ylab(metricName)
  if(metric == "FWER") {
    p_k <- p_k + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# n
# ------------------------------------------------------------------------------
if("n" %in% plotsReady) {
  methods <- rownames(scores_n[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_n)))
  rownames(df) <- methods
  colnames(df) <- names(scores_n)
  for(n in names(scores_n)) {
    for(method in methods) {
      df[method, n] <- scores_n[[n]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_n <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                               shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    ylab(metricName)
  if(metric == "FWER") {
    p_n <- p_n + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# nenv
# ------------------------------------------------------------------------------
if("nenv" %in% plotsReady) {
  methods <- rownames(scores_nenv[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_nenv)))
  rownames(df) <- methods
  colnames(df) <- names(scores_nenv)
  for(nenv in names(scores_nenv)) {
    for(method in methods) {
      df[method, nenv] <- scores_nenv[[nenv]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_nenv <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                  shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    ylab(metricName)
  if(metric == "FWER") {
    p_nenv <- p_nenv + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# nsubenv
# ------------------------------------------------------------------------------
if("nsubenvs" %in% plotsReady) {
  methods <- rownames(scores_nsubenvs[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_nsubenvs)))
  rownames(df) <- methods
  colnames(df) <- names(scores_nsubenvs)
  for(subenv in names(scores_nsubenvs)) {
    for(method in methods) {
      df[method, subenv] <- scores_nsubenvs[[subenv]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_nsubenvs <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                    shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    ylab(metricName)
  if(metric == "FWER") {
    p_nsubenvs <- p_nsubenvs + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# p
# ------------------------------------------------------------------------------
if("p" %in% plotsReady) {
  methods <- rownames(scores_p[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_p)))
  rownames(df) <- methods
  colnames(df) <- names(scores_p)
  for(p in names(scores_p)) {
    for(method in methods) {
      df[method, p] <- scores_p[[p]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_p <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                               shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    ylab(metricName)
  if(metric == "FWER") {
    p_p <- p_p + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# tau
# ------------------------------------------------------------------------------
if("tau" %in% plotsReady) {
  methods <- rownames(scores_tau[[1]][[metric]])
  df <- data.frame(matrix(NA, nrow = length(methods), ncol = length(scores_tau)))
  rownames(df) <- methods
  colnames(df) <- names(scores_tau)
  for(tau in names(scores_tau)) {
    for(method in methods) {
      df[method, tau] <- scores_tau[[tau]][[metric]][method, "avg"]
    }
  }
  df$method <- methods
  rowOrder <- c("random", "pooled regression", "GES", "LinGAM", "ICP", "nonlinearICP", "RICP")
  df$method <- factor(df$method, levels = rowOrder)
  df_melted <- melt(df, id = "method")
  p_tau <- ggplot(df_melted, aes(x = variable, y = value, group = method, colour = method, 
                                 shape = method)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.3) +
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
    ylab(metricName)
  if(metric == "FWER") {
    p_tau <- p_tau + geom_hline(yintercept = 0.05, linetype = 3)
  }
}

# JOINT PLOT
# ------------------------------------------------------------------------------
p_joint <- ggarrange(p_p, p_k, p_n, p_tau, p_nenv, p_nInter, p_interType, p_interMean, p_interStrength, 
          ncol = 3, nrow = 3, common.legend = TRUE, legend = "bottom")
p_joint
setwd(paste0(wdir, "fig/wsubenv/multipleInter"))
ggexport(p_joint, filename = paste0("jointPlot_RICP_", metric, ".pdf"))


