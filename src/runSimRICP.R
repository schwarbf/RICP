#' @title runSimRICP.R
#'
#' @description A helper function to compare the novel method RICP (Invariant 
#'   Causal Predictor under Random Effcts) to existing algorithms for causal 
#'   discovery. This function creates a random DAG with 'p' nodes and each node 
#'   contains on average 'k' edges to other nodes for 'nenv' environments 
#'   corresponding to interventions 'interType'. RICP is compared against 
#'   'random', 'pooled regression', 'GES', 'LinGAM', 'ICP' and 'nonlinearICP'. 
#'   Methods of comparison include 'FWER', 'average Jaccard similarity' and 
#'   'success probability'.
#'
#' @param p integer: The number of nodes the DAG should contain. 
#' @param k integer: The averge number of edges each node in the DAG contains.
#' @param nenv integer: The number of environments corresponding to the number of 
#'    interventions + 1 (observational data). 
#' @param renv vector: Contains two elements: The minimal and maximal number of 
#'    observations per subenvironment. 
#' @param rBeta vector: Contains two elements: The minimal and maximal bounds from 
#'    which to unifromly sample beta from. 
#' @param tau double: The variance parameter of the random effects. 
#' @param interType string: Choose from 'do', 'soft', simultaneous-noise', 
#'    'relaxed-do
#' @param interMean vector: This corresponds to the mean of the intervention. For 
#'    each intervention we sample the mean from 'interMean'. 
#' @param interStrength vector: This corresponds to the strength of the intervention. 
#'    For each intervention we sample the strength from 'interStrength' and 
#'    pre-multiply the intervention with it. 
#' @param nInter string: If "one" (default) then only one predictor variable per
#'    environment is being intervened on. If "multiple" the intervention is on 
#'    multiple predictors in the same environment. The number of interventions in
#'    this case is sampled uniformly from '1:(p-2)'. 
#' @param subenvs boolean: If TRUE then the RICP algorithm for subenvironments is
#'    used. 
#' @param nsubenvs integer/list: Each subenvironment corresponds to a realization 
#'    of the random effect within an environment. If 'nsubenvs' is an integer, then 
#'    each environment contains an equal number of subenvironments. Otherwise, 
#'    each list element contains a vector with the number of subenvironments per
#'    environment. Default is 'NULL' corresponding to no subenvironments. 
#' @param methods vector: The methods for which the causal discovery algorithm 
#'    should be run. Choose from 'random', 'pooled regression', 'GES', 'LinGAM', 
#'    'ICP', 'nonlinearICP' and 'RICP'
#' @param test string: If subenvs is TRUE, then test choices are "LRT-lme4" (default), 
#'   "LRT-nlme" or "Wald-test". The difference mostly lies in either the underlying 
#'   test (LRT vs Wald-test) and the package used to fit the LMM. The 'lme4' package 
#'   does not allow for concrete specification of the structure of the covariance 
#'   matrices for the noise and the random effects whereas 'nlme' does. The 
#'   drawback is that 'nlme' is much slower. 'lme4' fits a generic (unstructured) 
#'   covariance matrix. If subenvs is FALSE, then one can only specify the 
#'   package "lme4" or "nmle". Default is "LRT-lme4".
#' 
#' @return list: Contains the following objects
#'    DAG matrix: p x p adjacency matrix of the random DAG used.
#'    acceptedSets list: Elements are again lists where each entry in each list 
#'       contains the accepted sets for a specific parent node using a causal 
#'       discovery method. 
#'    FWER list: Elements are again lists where each entry in each list 
#'       corresponds to the FWER for the hypothesis test  for a specific parent 
#'       node using a causal discovery method. 
#'    avgJaccard list: Average Jaccard Similarity for each node and for each 
#'       causal discovery method. The average is taken over the nodes in the DAG.
#'    successProbability list: The success probability of identifying the correct 
#'       set of causal parents for all nodes in the DAG. 

runSimRICP <- function(p = 4, k = 2, nenv = 10, renv = c(50, 80), rBeta = c(-5, 5), 
                       tau = 1, alpha = 0.05, interType = "do",
                       interMean = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10), 
                       interStrength = c(0.1, 0.2, 0.5, 1, 2, 5, 10), 
                       nInter = "one", subenvs = FALSE, nsubenvs = NULL, 
                       methods = "RICP", test = "LRT-lme4") {
  # simulate a DAG
  if(subenvs) {
    tmp <- simDAGwsubenvs(p = p, k = k, nenv = nenv, nsubenvs = nsubenvs, 
                          renv = renv, rBeta = rBeta, tau = tau, 
                          interType = interType, interMean = interMean, 
                          interStrength = interStrength, nInter = nInter)
  } else{
    tmp <- simDAG(p = p, k = k, nenv = nenv, renv = renv, rBeta = rBeta, tau = tau, 
                  interType = interType, interMean = interMean, 
                  interStrength = interStrength, nInter = nInter)
  }
  data <- tmp$X
  adjacencyMat <- tmp$Adjacency
  ExpInd <- tmp$ExpInd
  interInd <- tmp$interInd
  acceptedSets <- list()
  
  # GES
  if("GES" %in% methods) {
    GES.accep <- list()
    score <- new("GaussL0penObsScore", data)
    GES.fit <- ges(score)
    cpdag <- GES.fit$essgraph$.in.edges
    for(i in 1:p) {
      plausPa <- cpdag[[i]]
      if(length(plausPa) > 0) {
        GES.accep[[i]] <- plausPa[which(sapply(plausPa, function(j) i %in% cpdag[[j]]) == F)]
      } else{
        GES.accep[[i]] <- integer(0)
      }
    }
    acceptedSets[["GES"]] <- GES.accep
  }
  
  # LinGAM
  if("LinGAM" %in% methods) {
    LinGAM.accep <- list()
    LinGAM.fit <- lingam(data)
    for(i in 1:p) {
      LinGAM.accep[[i]] <- which(LinGAM.fit$Bpruned[i, ] != 0)
    }
    acceptedSets[["LinGAM"]] <- LinGAM.accep
  }
  
  # random, pooled regression, ICP, RICP
  ICP.accep <- nonlinearICP.accep <- RICP.accep <- list()
  for(i in 1:p) {
    X <- data[, -i, drop = FALSE]
    Y <-  data[, i]
    
    # random 
    if("random" %in% methods) {
      random.accep <- list()
      estRand <- sample(0:1, size = 1, prob = c(1-alpha, alpha))
      random.accep[[i]] <- if(estRand == 0) integer(0) else sample((1:p)[-i], size = 1)
    }
    
    # pooled regression
    if("pooled regression" %in% methods) {
      pooledRegr.accep <- list()
      dat <- cbind(Y, X) %>% as.data.frame()
      colnames(dat) <- c("Y", colnames(X))
      form <- as.formula(paste("Y ~ ", paste(colnames(X), collapse= "+")))
      pvals <- summary(lm(form, data = dat))$coefficients[, 4]
      estPooledRegr <- which(pvals <= alpha/ncol(X))
      if("(Intercept)" %in% names(estPooledRegr)) {
        if(length(names(estPooledRegr)) == 1) {
          pooledRegr.accep[[i]] <- integer(0) # only intercept (empty set)
        } else{
          estPooledRegr <- estPooledRegr[-1]
        }
      } 
      if(!exists("retObj")) {
        pooledRegr.accep[[i]] <- gsub("X", "", names(estPooledRegr)) %>% as.integer()
      }
    }
    
    # ICP
    ExpIndICP <- ExpIndRICP <- ExpInd
    interIndEnv <- which(sapply(1:nenv, function(j) i %in% interInd[[j]]))
    if(subenvs) {
      ExpIndList <- lapply(1:nenv, 
                           function(j) {
                             if(j > 1) {
                               indStart <- length(ExpInd[1:(j-1)] %>% unlist())
                               indEnd <- length(ExpInd[1:j] %>% unlist())
                               (indStart + 1):indEnd
                             } else{
                               1:length(ExpInd[[1]])
                             }
                           }
      )
      if(sum(interIndEnv) > 0) {
        ExpIndRICP <- ExpIndRICP[-interIndEnv]
        ExpIndICP <- lapply((1:nenv)[-interIndEnv], 
                            function(j) rep(j, length(ExpInd[[j]]))) %>% unlist()
        interIndEnv <- ExpIndList[interIndEnv] %>% unlist()
      } else{
        ExpIndICP <- lapply(1:nenv, 
                            function(j) rep(j, length(ExpInd[[j]]))) %>% unlist()
      }
    } else{
      if(sum(interIndEnv) > 0) {
        interIndEnv <- which(ExpInd %in% interIndEnv)
        ExpIndICP <- ExpIndRICP <- ExpInd[-interIndEnv]
      }
    }
    if(sum(interIndEnv) > 0){
      X <- X[-interIndEnv, , drop = FALSE]
      Y <- Y[-interIndEnv]
    }
    if("ICP" %in% methods) {
      ICP.fit <- ICP(X, Y, ExpIndICP, alpha = alpha, test = "normal",
                     showAcceptedSets = F, showCompletion = F, stopIfEmpty = T)
      estICP <- Reduce(intersect, ICP.fit$acceptedSets)
      if(is.null(estICP)) {
        ICP.accep[[i]] <- integer(0)
      } else{
        ICP.accep[[i]] <- gsub("X", "", colnames(X)[estICP]) %>% as.integer()
      }
    }
    
    # nonlinearICP
    if("nonlinearICP" %in% methods) {
      nonlinearICP.fit <- nonlinearICP(X, Y, as.factor(ExpIndICP), alpha = alpha)
      estNonlinearICP <- Reduce(intersect, nonlinearICP.fit$acceptedSets)
      if(length(estNonlinearICP) == 0) {
        nonlinearICP.accep[[i]] <- integer(0)
      } else{
        nonlinearICP.accep[[i]] <- gsub("X", "", colnames(X)[estNonlinearICP]) %>% as.integer()
      }
    }
    
    # RICP
    if("RICP" %in% methods) {
      RICP.fit <- RICP(X, Y, ExpIndRICP, alpha = alpha, subenvs = subenvs, 
                       showAcceptedSets = F, showProgress = F, stopIfEmpty = T, 
                       test = test)
      if(is.null(RICP.fit$estimate)) {
        RICP.accep[[i]] <- integer(0)
      } else{
        RICP.accep[[i]] <- gsub("X", "", RICP.fit$estimate) %>% as.integer()
      }
    }
  }
  if("random" %in% methods) {acceptedSets[["random"]] <- random.accep}
  if("pooled regression" %in% methods) {acceptedSets[["pooled regression"]] <- pooledRegr.accep}
  if("ICP" %in% methods) {acceptedSets[["ICP"]] <- ICP.accep}
  if("nonlinearICP" %in% methods) {acceptedSets[["nonlinearICP"]] <- nonlinearICP.accep}
  if("RICP" %in% methods) {acceptedSets[["RICP"]] <- RICP.accep}
  
  # compute metrics
  fwer <- successProb <- avgJaccard <- list()
  for(method in methods) {
    fwer[[method]] <- successProb[[method]] <- avgJaccard[[method]] <- 0
    fwerCount <- 0
    jaccard <- c()
    for(i in 1:p) {
      truePa <- which(adjacencyMat[i, ] == 1) %>% as.vector()
      if(length(acceptedSets[[method]][[i]]) > 0) {
        fwerCount <- fwerCount + sapply(acceptedSets[[method]][[i]], 
                                        function(j) {!(j %in% truePa)}) %>% sum()
        jaccardNum <- length(intersect(truePa, acceptedSets[[method]][[i]]))
        jaccardDenom <- length(union(truePa, acceptedSets[[method]][[i]]))
        jaccard[i] <- jaccardNum/jaccardDenom
        lenInter <- length(intersect(acceptedSets[[method]][[i]], truePa))
        if(lenInter == length(truePa) && lenInter == length(acceptedSets[[method]][[i]])) {
          successProb[[method]] <- successProb[[method]] + 1
        }
      } else{
        fwerCount <- fwerCount
        jaccard[i] <- if(length(truePa) == 0) 0 else 0
        if(length(truePa) == 0) {successProb[[method]] <- successProb[[method]] + 1} 
      }
    }
    fwer[[method]] <- 1/(p*(p-1))*fwerCount
    avgJaccard[[method]] <- mean(jaccard)
    successProb[[method]] <- successProb[[method]]/p
  }
  retObj <- list(DAG = adjacencyMat, acceptedSets = acceptedSets, FWER = fwer, 
                 avgJaccard = avgJaccard, successProbability = successProb)
  return(retObj)
}