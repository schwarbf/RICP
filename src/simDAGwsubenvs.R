#' @title simDAGwsubenvs.R
#'
#' @description Simulates a DAG with 'p' nodes and each node contains on average
#'    'k' edges to other nodes for 'nenv' environments corresponding to 
#'    interventions 'interType'. Each environment contains 'nsubenvs' number of 
#'    subenvironments corresponding to realizations of the random effects. 
#'
#' @param p integer: The number of nodes the DAG should contain. 
#' @param k integer: The averge number of edges each node in the DAG contains.
#' @param nenv integer: The number of environments corresponding to the number of 
#'    interventions + 1 (observational data). 
#' @param nsubenvs integer/list: Each subenvironment corresponds to a realization 
#'    of the random effect within an environment. If 'nsubenvs' is an integer, then 
#'    each environment contains an equal number of subenvironments. Otherwise, 
#'    each list element contains a vector with the number of subenvironments per
#'    environment. 
#' @param renv vector: Contains two elements: The minimal and maximal number of 
#'    observations per subenvironment. 
#' @param rBeta vector: Contains two elements: The minimal and maximal bounds from 
#'    which to unifromly sample beta from. 
#' @param tau double: The variance parameter of the random effects. 
#' @param interType string: Choose from 'do', 'soft', simultaneous-noise', 
#'    'relaxed-do.
#' @param interMean vector: This corresponds to the mean of the intervention. For 
#'    each intervention we sample the mean from 'interMean'. 
#' @param interStrength vector: This corresponds to the strength of the intervention. 
#'    For each intervention we sample the strength from 'interStrength' and 
#'    pre-multiply the intervention with it. 
#' @param nInter string: If "one" (default) then only one predictor variable per
#'    environment is being intervened on. If "multiple" the intervention is on 
#'    multiple predictors in the same environment. The number of interventions in
#'    this case is sampled uniformly from '1:(p-2)'. 
#' 
#' @return list: Contains the following elements: 
#'    Adjacency matrix: A p x p adjacency matrix corresponding to the DAG.
#'    X matrix: A (nenv x subenv) x p matrix containing the observations of all 
#'       nodes .
#'    Beta matrix: A p x p matrix containing the fixed effects. 
#'    interInd list: Each list element is a vector containing the nodes which 
#'       were intervened on in the respective environment. 
#'    ExpInd list: Each element contains an indicator which observation corresponds
#'       to which subenvironemnt. 

simDAGwsubenvs <- function(p = 4, k = 2, nenv = 10, nsubenvs = 10, renv = c(50, 80), 
                           rBeta = c(-5, 5), tau = 1, interType = "do",
                           interMean = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10), 
                           interStrength = c(0.1, 0.2, 0.5, 1, 2, 5, 10), 
                           nInter = "one") {
  seed <- sample(1:1000, 1)
  set.seed(seed)
  
  # sampling DAG
  Adj <- matrix(0, ncol = p, nrow = p)
  Adj[upper.tri(Adj)] <- sample(0:1, p*(p-1)/2, replace = TRUE, 
                                prob = c(1-(k/(p-1)), k/(p-1)))
  colnames(Adj) <- rownames(Adj) <- paste0("X", 1:p)
  
  # sample fixed effects
  Sign <- Beta <- B <- Adj
  ind <- which(Sign == 1, arr.ind = T)
  Sign[ind] <- sample(c(-1, 1), length(Sign[ind]), replace = TRUE)
  Beta[ind] <- runif(length(Sign[ind]), min = rBeta[1], max = rBeta[2])
  
  # sample data within environments 
  n <- lapply(1:nenv, function(i) ceiling(runif(nsubenvs, min = renv[1], max = renv[2])))
  eps <- matrix(0, ncol = p, nrow = sum(unlist(n)))
  X <- matrix(0, ncol = p, nrow = sum(unlist(n))) 
  colnames(X) <- colnames(Adj)
  B <- list()
  interInd <- list()
  i <- 1
  for(env in 1:nenv) {
    if(env != 1) {
      type <- sample(interType, 1)
      if(nInter == "one") {
        interInd[[env]] <- sample(1:p, 1)
      } else{
        if(p - 2 > 0){
          numInter <- sample(c("one", "multiple"), size = 1, prob = c(1/6, 5/6))
          if(numInter == "one") {
            interInd[[env]] <- sample(1:p, 1)
          } else{
            nInter <- sample(1:(p-2), 1)
            interInd[[env]] <- sample(1:p, nInter) 
          }
        } else{
          interInd[[env]] <- sample(1:p, 1)
        }
      }
      interInd[[env]] <- sample(1:p, size = 1)
      interMean <- if(length(interMean) > 1) sample(interMean, size = 1) else interMean
      interStrength <- if(length(interStrength) > 1) sample(interStrength, size = 1) else interStrength
      interSig <- interMean + interStrength*runif(1)
    }
    B[[env]] <- list()
    for(subenv in 1:nsubenvs) {
      if(env == 1 && subenv == 1) {
        rowInd <- 1:n[[env]][1]
      } else{
        rowInd <- (max(rowInd) + 1):(max(rowInd) + n[[env]][subenv])
      }
      sdEps <- sample(c(0.1, 0.2, 0.5, 1, 2, 5, 10), 1)
      eps[rowInd, ] <- matrix(rnorm(p*n[[env]][subenv], sd = sdEps), 
                              nrow = n[[env]][subenv], ncol = p)
      B[[env]][[subenv]] <- Sign
      B[[env]][[subenv]][ind] <- rnorm(length(Sign[ind]), sd = tau)
      for(j in p:1) {
        X[rowInd, j] <- X[rowInd, j:p, drop = FALSE] %*% 
          (Beta[j, j:p] + B[[env]][[subenv]][j, j:p]) + eps[rowInd, j] 
        if(env != 1 && j %in% interInd[[env]]) {
          if(type == "do") {
            X[rowInd, j] <- interSig
          } else if(type == "soft") {
            X[rowInd, j] <- X[rowInd, j] + interSig
          } else if(type == "simultaneous-noise") {
            X[rowInd, j] <- X[rowInd, j:p, drop = FALSE] %*% 
              (Beta[j, j:p] + B[[env]][[subenv]][j, j:p]) + interSig*eps[rowInd, j]
          } else if(type == "relaxed-do") {
            X[rowInd, j] <- interMean + interStrength*runif(length(rowInd))
          } else{
            stop("Intervention type is not supported!")
          }
        }
      }
    }
  }
  ExpInd <- lapply(1:nenv, 
                   function(i) lapply(1:nsubenvs, function(j) {rep(j, n[[i]][j])}) %>% 
                     unlist())
  
  return(list(Adjacency = Adj, X = X, Beta = Beta, interInd = interInd, 
              ExpInd = ExpInd))
}