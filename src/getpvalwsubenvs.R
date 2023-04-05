#' @title getpvalwsubenvs.R
#'
#' @description Tests whether a given subset of predictors satisfies the 
#'    hypothesis of 'invariance under random effects' using subenvironments. 
#'
#' @param X data.frame: Contains the fixed-effects predictor.
#' @param Y vector: Contains the response value for all environments and 
#'    subenvironments. 
#' @param ExpInd list: Each element contains a vector which indicates which 
#'    observation belongs to which subenvironment. 
#' @param alpha float: Determines coverage of the confidence regions.
#' @param test string: Choose from the type of test from "F-test", "LRT-nlme" or
#'    "LRT-lme04". For the F-test, a Wald test using the package "nlme" is done.
#'    For the other two, a likelihood ratio test (LRT) is used using different 
#'    packages. The package "nlme" allows for specifying explicit covariance 
#'    structures, whereas the package "lme4" does not. Defaults to "LRT-lme4".
#' 
#' @return list: Contains the following elements: 
#'   res string: Test result; either 'rejected' or 'not rejected'.
#'   coeff vector: Fixed-effects estimates using data from the observational data.
#'   coeffSE data.frame: Contains the standard errors (SE) of the estimated 
#'      fixed-effect coefficients also using the fit on the observational data. 
#'   contInt data.frame: Contains the confidence intervals for the fixed-effects 
#'      using the fit on the observational data. 
#'   pvals vector: Contains the p-values from each hypothesis test. 
 
getpvalwsubenvs <- function (X, Y, ExpInd, alpha = 0.05, test = "LRT-lme4") {
  # set the optimizer control parameters 
  ctrl <- lmerControl(optCtrl = list(x_tol_abs = 1e-08, f_tol_abs = 1e-08))
  
  # extract confidence intervals (fit on observational distribution)
  cols <- c("Y", colnames(X))
  ExpIndObs <- ExpInd[[1]]
  rowIndObs <- 1:length(ExpInd[[1]])
  data <- cbind(Y[rowIndObs], X[rowIndObs, , drop = FALSE]) %>% 
    as.data.frame() %>% 
    `colnames<-` (cols)
  formLMM <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1 +",
                              paste0("(", paste0(cols[-1], collapse= " + ")), 
                              "-1 | ExpIndObs)"))
  lmer.fit <- lmer(formLMM, data = data, REML = TRUE, control = ctrl) %>% 
    suppressMessages()
  coeff <- fixef(lmer.fit)
  coeffSE <- summary(lmer.fit)$coefficients[, 2]
  confInt <- confint(lmer.fit, parm = "beta_", method = "Wald", level = 1 - alpha) %>% 
    suppressMessages() 
  colnames(confInt) <- c("lower", "upper")
  if("intercept" %in% colnames(X) && ncol(X) > 1) {confInt <- confInt[-1, , drop = FALSE]}
  
  # perform hypothesis test for observational data against the rest
  nenv <- length(ExpInd)
  envs <- 2:nenv
  pvalvec <- numeric(nenv - 1)
  notRejected <- TRUE
  i <- 1
  while(notRejected & i <= length(envs)) {
    env <- envs[i]
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
    rowIndTest <- ExpIndList[[env]]
    ExpIndTest <- c(ExpIndObs, ExpInd[[env]] + length(unique(ExpIndObs)))
    Xtest <- bdiag(X[rowIndObs, , drop = FALSE], X[rowIndTest, , drop = FALSE]) %>% 
      as.matrix()
    Ztest <- rbind(X[rowIndObs, , drop = FALSE], X[rowIndTest, , drop = FALSE])
    Ytest <- c(Y[rowIndObs], Y[rowIndTest])
    colsTest <- c("Y", sapply(1:ncol(X), function(i){paste0(colnames(X)[i], "_Obs")}),
                  sapply(1:ncol(X), function(i){paste0(colnames(X)[i], "_Test")}))
    
    if(test == "LRT-lme4") {
      Xtest[(length(rowIndObs) + 1):nrow(Xtest), 1:ncol(X)] <- -X[rowIndTest, , drop = FALSE]
      dataTest <- cbind(Ytest, Xtest, Ztest) %>% 
        as.data.frame() %>% 
        `colnames<-` (c(colsTest, cols[-1]))
      formAll <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1 + ",
                                  paste0("(", paste0(cols[-1], collapse= " + ")), 
                                  "-1 | ExpIndTest)"))
      formNested <- as.formula(paste("Y ~ ", paste0(colsTest[-1], collapse= " + "), 
                                     " -1 +", paste0("(", paste0(cols[-1], collapse= " + ")), 
                                     "-1 | ExpIndTest)"))
      lmer.fitNested <- lmer(formNested, data = dataTest, REML = FALSE, control = ctrl) %>% 
        suppressMessages()
      lmer.fitAll <- lmer(formAll, data = dataTest, REML = FALSE, control = ctrl) %>% 
        suppressMessages()
      pvalvec[i] <- anova(lmer.fitNested, lmer.fitAll)[["Pr(>Chisq)"]][2]
    } else if(test == "LRT-nlme") {
      Xtest[(length(rowIndObs) + 1):nrow(Xtest), 1:ncol(X)] <- -X[rowIndTest, , drop = FALSE]
      dataTest <- cbind(Ytest, Xtest, Ztest) %>% 
        as.data.frame() %>% 
        `colnames<-` (c(colsTest, cols[-1]))
      ExpIndTest <<- ExpIndTest
      formAll <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1"))
      formRdmAll <- as.formula(paste("~", paste(cols[-1], collapse= "+"), " -1"))
      formFxd <- as.formula(paste(" ~ ", paste(cols[-1], collapse= "+"), " -1 | ExpIndTest"))
      varFixed <- varIdent(form = formFxd)
      formNested <- as.formula(paste("Y ~ ", paste(colsTest[-1], collapse= "+"), " -1"))
      formRdmNested <- as.formula(paste("~", paste(cols[-1], collapse= "+"), " -1"))
      lme.fitAll <- lmeFit(formAll, random = list(ExpIndTest = pdIdent(formRdmAll)),
                           weights = varFixed, data = dataTest, nrestarts = 5)
      lme.fitNested <- lmeFit(formNested, random = list(ExpIndTest = pdIdent(formRdmNested)),
                              weights = varFixed, data = dataTest, nrestarts = 5)
      pvalvec[i] <- anova(lme.fitNested, lme.fitAll)[["p-value"]][2]
    } else if(test == "Wald-test") {
      dataTest <- cbind(Ytest, Xtest, Ztest) %>% 
        as.data.frame() %>% 
        `colnames<-` (c(colsTest, cols[-1]))
      ExpIndTest <<- ExpIndTest
      formTest <- as.formula(paste("Y ~ ", paste(colsTest[-1], collapse= "+"), " -1"))
      formFxdTest <- as.formula(paste(" ~ ", paste(colsTest[-1], collapse= "+"), 
                                      " -1 | ExpIndTest"))
      formRdmTest <- as.formula(paste("~", paste(cols[-1], collapse= "+"), " -1"))
      varFixedTest <- varIdent(form = formFxdTest)
      lme.fitTest <- lmeFit(formTest, random = list(ExpIndTest = pdIdent(formRdmTest)),
                            weights = varFixedTest, data = dataTest, nrestarts = 5)
      if(ncol(X) > 1) {
        L = c(0, rep(1, ncol(X) - 1), 0, rep(-1, ncol(X) - 1))
      } else{
        L = c(1, -1)
      }
      pvalvec[i] <- anova(lme.fitTest, type = "sequential", L = L)[["p-value"]]
    } else{
      stop(paste0("Test ", test, " is currently not implemented."))
    }
    
    alphaBonfCor <- alpha/(nenv - 1)
    if(pvalvec[i] <= alphaBonfCor) {notRejected <- FALSE}
    i <- i + 1
  }
  result <- if(notRejected) {"not rejected"} else{"rejected"}

  return(list(res = result, coeff = coeff, coeffSE = coeffSE, confInt = confInt, 
              pvals = pvalvec))
}
