#' @title RICP.R
#'
#' @description Algorithm to estimate the identifiable causal predictors under
#'    random effects. 
#'
#' @param X data.frame: Contains the fixed-effects predictor.
#' @param Y vector: Contains the response value for all environments.
#' @param ExpInd vector/list: Contains the indicators for the environments.
#' @param alpha float: Determines coverage of the confidence regions.
#' @param subenvs boolean: If TRUE, then the RICP method for subenvironments is 
#'    applied.
#' @param showAcceptedSets boolean: If TRUE, print out information about accepted
#'    sets of variables. 
#' @param showProgress boolean: If TRUE, print out progress of computation.
#' @param stopIfEmpty boolean: If TRUE, the algorithm stops if the intersection 
#'    of accepted sets is empty. 
#' @param test string: Ifsubenvs is TRUE, then a hypothesis test can be specified. 
#'   The choices are "LRT-lme4" (default), "LRT-nlme" or "Wald-test". The 
#'   difference mostly lies in either the underlying test (LRT vs Wald-test) and 
#'   the package used to fit the LMM. The 'lme4' package does not allow for concrete
#'   specification of the structure of the covariance matrices for the noise and 
#'   the random effects whereas 'nlme' does. The drawback is that 'nlme' is much 
#'   slower. 'lme4' fits a generic (unstructured) covariance matrix.
#' 
#' @return list: Contains the following objects: 
#'   estimate vector: Contains the estimate of the predictors that are invariant
#'      under random effects. 
#'   confInt data.frame: Contains the confidence intervals for the fixed-effects 
#'      using the fit on the observational data. 
#'   alpha double: Significance level of chosen hypothesis test. 
#'   colnames vector: Names of all plausible causal predictors under random effects.
#'   dimX vector: First element contains number of rows of input 'X' and second
#'      element contains the number of columns of input 'X'.
#'   coeff vector: Fixed-effects estimates using data from the observational data.
#'   coeffSE data.frame: Contains the standard errors (SE) of the estimated 
#'      fixed-effect coefficients also using the fit on the observational data. 
#'   acceptedSets list: Each element contains the vector of accepted sets. 
#'   testResults list: Contains the test results for each possible set of 
#'      predictors. 
#'   stopIfEmpty boolean: TRUE if the algorithm was stopped if the intersection of 
#'      plausible causal predictors was empty. FALSE otherwise. 
#'   nEnv integer: Number of environments. 
#' @export

RICP <- function(X, Y, ExpInd, alpha = 0.05, subenvs = FALSE, showAcceptedSets = TRUE, 
                 showProgress = TRUE, stopIfEmpty = TRUE, test = "LRT-lme4"){
  # check inputs: data 
  if (is.vector(X) & is.numeric(X)) 
    X <- matrix(X, ncol = 1)
  if (!is.matrix(X) & !is.data.frame(X)) 
    stop("'X' must be a matrix or data frame")
  if (length(ucol <- unique(colnames(X))) < min(3, ncol(X))) 
    colnames(X) <- paste("X", 1:ncol(X), sep = "_")
  if (!is.vector(Y)) 
    stop("'Y' must be a vector")

  # Initialization 
  accepted <- rejected <- candidates <- list()
  intersection <- NULL
  confInt <- matrix(NA, nrow = 2, ncol = ncol(X))
  colnames(confInt) <- colnames(X)
  coeff <- coeffSE <- vector(mode = "list", length = ncol(X))
  testResults <- list()
  continue <- TRUE
  
  # adding intercept to model 
  X <- cbind(1, X) 
  colnames(X) <- c('intercept', colnames(X)[-1])
  
  # create candidates (powerset)
  toBin <- function(n, bitsToKeep) {
    binVec <- rev(as.numeric(intToBits(n)))
    return(binVec[- c(1:(length(binVec) - bitsToKeep))])
  }
  candidates <- list()
  bitsToKeep <- ceiling(log2(2^length(1:ncol(X)) - 1))
  for (j in ((1:2^length(1:(ncol(X)-1))) - 1)) {
    candidates[[j + 1]] <- (colnames(X))[which(toBin(j, bitsToKeep) == 1)]
  }
  candidates <- unique(candidates)
  len <- sapply(candidates, length)
  candidates <- candidates[order(len[1:length(candidates)])]
  
  # Test all candidates (except empty set)
  setIter <- 0
  ncandidates <- length(candidates)
  while(continue && setIter < ncandidates) {
    setIter <- setIter + 1
    candidate <- candidates[[setIter]]
    if(length(candidate) == 0) { # empty set 
      if(!subenvs) {
        tmp <- getpval(X[, 1, drop = FALSE], Y, ExpInd, alpha = alpha)
      } else{
        tmp <- getpvalwsubenvs(X[, 1, drop = FALSE], Y, ExpInd, alpha = alpha, 
                               test = test)
      }
      if (tmp$res == "not rejected") {
        confInt <- matrix(0, nrow = 2, ncol = ncol(X) - 1)
        colnames(confInt) <- colnames(X)[-1]
        testResults <- append(testResults, paste0("Empty set: ", tmp$res))
        intersection <- character(0)
        if(stopIfEmpty)
          continue <- FALSE
        if (showAcceptedSets)
          cat(paste("accepted empty set \n"))
      }
    } else{ # any other set
      notcandidate <- colnames(X)[-c(1, which(colnames(X) %in% candidates[[setIter]]))]
      if(!subenvs) {
        tmp <- getpval(X[, c("intercept", candidate), drop = FALSE], Y, ExpInd, alpha = alpha)
      } else{
        tmp <- getpvalwsubenvs(X[, c("intercept", candidate), drop = FALSE], Y, ExpInd, alpha = alpha, 
                               test = test)
      }
      testResults <- append(testResults, paste0(paste(candidate, collapse = ", "), ": ", tmp$res))
      if (tmp$res == "not rejected") {
        if (is.null(intersection)) {
          intersection <- candidate
        } else {
          intersection <- intersect(intersection, candidate)
        }
        if (length(intersection) == 0 && stopIfEmpty) 
          continue <- FALSE
        accepted[[length(accepted) + 1]] <- candidate
        if (showAcceptedSets){
          cat(paste0("accepted set of variables: ", paste(candidate, collapse = ","), "\n"))
        }
        confInt[1, candidate] <- pmin(confInt[1, candidate, drop = FALSE], 
                                      tmp$confInt[, 1, drop = FALSE], na.rm = TRUE)
        confInt[2, candidate] <- pmax(confInt[2, candidate, drop = FALSE], 
                                      tmp$confInt[, 2, drop = FALSE], na.rm = TRUE)
        for (pred in candidate) {
          coeff[[pred]] <- c(coeff[[pred]], tmp$coeff[1 + which(candidate == pred)])
          coeffSE[[pred]] <- c(coeffSE[[pred]], tmp$coeffSE[1 + which(candidate == pred)])
        }
        if (length(notcandidate) >= 1) {
          confInt[1, notcandidate] <- pmin(confInt[1, notcandidate, drop = FALSE],
                                           0, na.rm = TRUE)
          confInt[2, notcandidate] <- pmax(confInt[2, notcandidate, drop = FALSE],
                                           0, na.rm = TRUE)
        }
      }
    }
    if (showProgress && ncandidates > 0) {
      cat(paste0("*** ", round(100 * setIter/ncandidates), 
                 "% complete: tested ", setIter, " out of ", ncandidates, 
                 " sets of variables \n"))
    }
  }
  retobj <- list(estimate = intersection, confInt = confInt, alpha = alpha, 
                 colnames = colnames(confInt), dimX = dim(X), coeff = coeff, 
                 coeffSE = coeffSE,  acceptedSets = accepted, testResults = testResults, 
                 stopIfEmpty = stopIfEmpty, nEnv = length(ExpInd))
  class(retobj) <- "InvariantCausalPredictionUnderRandomEffects"
  return(retobj)
}