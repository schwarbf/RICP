#' @title getpval.R
#'
#' @description 
#'
#' @param X data.frame: Contains the fixed-effects predictor.
#' @param Y vector: Contains the response value for all environments.
#' @param ExpInd list: Contains the indicators which environment observation 
#'    belongs to which partition. 
#' @param alpha float: Determines coverage of the confidence regions.
#' 
#' @return 
 
getpval <- function (X, Y, ExpInd, alpha = 0.05) {
  # extract confidence intervals (fit on all data)
  cols <- c("Y", colnames(X))
  data <- cbind(Y, X) %>% as.data.frame() %>% `colnames<-` (cols)
  formLMM <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1 +",
                              paste0("(", paste0(cols[-1], collapse= " + ")), "-1 | ExpInd)"))
  lmer.fit <- lmer(formLMM, data = data, REML = TRUE) %>% suppressMessages()
  # ============================================================================
  # PRINTING
  # print(colnames(X))
  # ============================================================================
  coeff <- fixef(lmer.fit)
  coeffSE <- summary(lmer.fit)$coefficients[, 2]
  confInt <- confint(lmer.fit, parm = "beta_", method = "Wald", level = 1 - alpha) %>% 
    suppressMessages() 
  colnames(confInt) <- c("lower", "upper")
  if("intercept" %in% colnames(X) && ncol(X) > 1) {confInt <- confInt[-1, , drop = FALSE]}
  
  
  # perform hypothesis test for one against the rest
  # envs <- unique(ExpInd)
  envs <- 1
  notRejected <- TRUE
  i <- 1
  alphaBonfCor <- alpha/length(envs)
  while(notRejected & i <=length(envs)) {
    env <- envs[i]
    indEnv <- which(ExpInd == env)
    
    # check if matrix data[indEnv, ] has full column rank 
    if(!("intercept" %in% colnames(X))) {
      dataTest <- data[indEnv, , drop = FALSE]
      columnSums <- colSums(dataTest)
      colRankTest <- sapply(2:ncol(data), function(i) columnSums[i]/nrow(dataTest) == dataTest[1, i]) %>% sum()
      if(colRankTest == (ncol(data) - 1)) {
        stop(paste0("The matrix of predictors in environment ", env, " has not full ", 
                    "column rank. Did you apply a do-intervention?"))
      }
    } 

    # estimate 'beta + b' on single environment
    formLM <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1 "))
    lm.fit <- lm(formLM, data = data[indEnv, , drop = FALSE])
    confInt_betaPlusB <- confint(lm.fit, alpha = alphaBonfCor/2)
    colnames(confInt_betaPlusB) <- c("lower", "upper")
    
    # get estimate of tau (package 'lme4')
    # ExpIndTest <- ExpInd[-which(ExpInd == env)]
    # formLMMTest <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1 +",
    #                                 paste0("(", paste0(cols[-1], collapse= " + ")), "-1 | ExpIndTest)"))
    # lmer.fit <- lmer(formLMMTest, data = data[-indEnv, ], REML = TRUE) %>% suppressMessages()
    # confInt_beta <- confint(lmer.fit, parm = "beta_", method = "Wald", level = 1 - alphaBonfCor/2) %>% suppressMessages()
    # colnames(confInt_beta) <- c("lower", "upper")
    # tauHat <- VarCorr(lmer.fit)$ExpIndTest[1, 1] %>% sqrt()
    
    # get estimate of tau (package 'nlme')
    # ==========================================================================
    ExpIndTest <<- ExpInd[-which(ExpInd == env)]
    formLMMTest <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1"))
    formFxdTest <- as.formula(paste(" ~ ", paste(cols[-1], collapse= "+"), " -1 | ExpIndTest"))
    formRdmTest <- as.formula(paste("~", paste(cols[-1], collapse= "+"), " -1"))
    varFixedTest <- varIdent(form = formFxdTest)
    # lme.fit <- lmeFit(formLMMTest, random = list(ExpIndTest = pdIdent(formRdmTest)),
    #                   weights = varFixedTest, data = data[-indEnv, ], nrestarts = 5)
    lme.fit <- lmeFit(formLMMTest, random = list(ExpIndTest = pdIdent(formRdmTest)),
                      weights = NULL, data = data[-indEnv, ], nrestarts = 5)
    confInt_beta <- intervals(lme.fit, which = "fixed", level = 1- alphaBonfCor/2)$fixed %>%
      as.data.frame() %>%
      select(c("lower", "upper"))
    tauHat <- VarCorr(lme.fit)[1, "StdDev", drop = FALSE] %>% as.double()
    rm(ExpIndTest, envir = .GlobalEnv)
    # ==========================================================================
    
    # construct confidence intervals for '(beta + b) - b'
    confInt_betaPlusBMinusB <- confInt_beta
    confInt_betaPlusBMinusB[, "lower"] <- confInt_betaPlusB[, "lower"] - qnorm(1-alphaBonfCor/4)*tauHat
    confInt_betaPlusBMinusB[, "upper"] <- confInt_betaPlusB[, "upper"] + qnorm(1-alphaBonfCor/4)*tauHat
    
    # check if 'confInt_beta' and 'confInt_betaPlusBMinusB' overlaps
    checkOverlap <- function(i) {
      confInt_betaPlusBMinusB[i, "upper"] < confInt_beta[i, "lower"] || 
        confInt_beta[i, "upper"] < confInt_betaPlusBMinusB[i, "lower"]
    }
    checkOverlapOut <- sapply(1:nrow(confInt_beta), checkOverlap) %>% sum()
    if(is.na(checkOverlapOut)) {browser()}
    if(checkOverlapOut > 0) {notRejected <- FALSE}
    print(tauHat)
    print(confInt_beta)
    print(confInt_betaPlusB)
    print(confInt_betaPlusBMinusB)
    browser()
    # ==========================================================================
    # PRINTING
    # if(checkOverlapOut > 0) {
    #   print(tauHat)
    #   print(confInt_beta)
    #   print(confInt_betaPlusBMinusB)
    # }
    # ==========================================================================
    # i <- i + 1
  }
  result <- if(notRejected) {"not rejected"} else{"rejected"}

  return(list(out = result, coeff = coeff, coeffSE = coeffSE, confInt = confInt))
}
