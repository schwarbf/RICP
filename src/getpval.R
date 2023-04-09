#' @title getpval.R
#'
#' @description Tests whether a given subset of predictors satisfies the 
#'    hypothesis of 'invariance under random effects' without using subenvironments.
#'
#' @param X data.frame: Contains the fixed-effects predictor.
#' @param Y vector: Contains the response value for all environments.
#' @param ExpInd list: Contains the indicators which environment observation 
#'    belongs to which partition. 
#' @param alpha float: Determines coverage of the confidence regions.
#' @param test string: Choose from the packages "lme4" or "nmle". The package "nlme" 
#'    allows for specifying explicit covariance structures, whereas the package 
#'    "lme4" does not. Defaults to "lme4". (Due to computational reasons, the 
#'    package "nlme" fits in this current version a diagonal noise covariance matrix
#'    with same parameter for all environments. Uncomment the lines 86/87 to fit
#'    covariance matrices with different parameters for each environment.)
#' 
#' @return list: Contains the following elements: 
#'   res string: Test result; either 'rejected' or 'not rejected'.
#'   coeff vector: Fixed-effects estimates using data from all environments.
#'   coeffSE data.frame: Contains the standard errors (SE) of the estimated 
#'      fixed-effect coefficients also using the fit on all environments.
#'   contInt data.frame: Contains the confidence intervals for the fixed-effects 
#'      using the fit on all environments. 
 
getpval <- function (X, Y, ExpInd, alpha = 0.05, test = "lme4") {
  # extract confidence intervals (fit on all data)
  cols <- c("Y", colnames(X))
  data <- cbind(Y, X) %>% as.data.frame() %>% `colnames<-` (cols)
  formLMM <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1 +",
                              paste0("(", paste0(cols[-1], collapse= " + ")), "-1 | ExpInd)"))
  lmer.fit <- lmer(formLMM, data = data, REML = TRUE) %>% suppressMessages()
  coeff <- fixef(lmer.fit)
  coeffSE <- summary(lmer.fit)$coefficients[, 2]
  confInt <- confint(lmer.fit, parm = "beta_", method = "Wald", level = 1 - alpha) %>% 
    suppressMessages() 
  colnames(confInt) <- c("lower", "upper")
  if("intercept" %in% colnames(X) && ncol(X) > 1) {confInt <- confInt[-1, , drop = FALSE]}
  
  
  # perform hypothesis test for one against the rest
  envs <- 1 # change to 'envs <- unique(ExpInd)' to do all against all
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
      colRankTest <- sapply(2:ncol(data), 
                            function(i) columnSums[i]/nrow(dataTest) == dataTest[1, i]) %>% sum()
      if(colRankTest == (ncol(data) - 1)) {
        stop(paste0("The matrix of predictors in environment ", env, " has not full ", 
                    "column rank. Did you apply a do-intervention?"))
      }
    } 

    # estimate 'beta + b' on single environment
    formLM <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1 "))
    lm.fit <- lm(formLM, data = data[indEnv, , drop = FALSE])
    # confInt_betaPlusB <- confint(lm.fit, alpha = alphaBonfCor/2)
    # colnames(confInt_betaPlusB) <- c("lower", "upper")
    confInt_beta_lm <- confint(lm.fit, alpha = alphaBonfCor/2)
    colnames(confInt_beta_lm) <- c("lower", "upper")
    
    # estimate tau on rest of environments
    if(test == "lme4"){
      ExpIndTest <- ExpInd[-which(ExpInd == env)]
      formLMMTest <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), 
                                      " -1 +", paste0("(", paste0(cols[-1], collapse= " + ")), 
                                      "-1 | ExpIndTest)"))
      lmer.fit <- lmer(formLMMTest, data = data[-indEnv, ], REML = TRUE) %>% 
        suppressMessages()
      # confInt_beta <- confint(lmer.fit, parm = "beta_", method = "Wald",
      #                         level = 1 - alphaBonfCor/2) %>% suppressMessages()
      confInt_beta_lmm <- confint(lmer.fit, parm = "beta_", method = "Wald",
                                  level = 1 - alphaBonfCor/2) %>% suppressMessages()
      colnames(confInt_beta_lmm) <- c("lower", "upper")
      # tauHat <- VarCorr(lmer.fit)$ExpIndTest[1, 1] %>% sqrt()
    } else if(test == "nmle") {
      ExpIndTest <<- ExpInd[-which(ExpInd == env)]
      formLMMTest <- as.formula(paste("Y ~ ", paste(cols[-1], collapse= "+"), " -1"))
      formFxdTest <- as.formula(paste(" ~ ", paste(cols[-1], collapse= "+"), " -1 | ExpIndTest"))
      formRdmTest <- as.formula(paste("~", paste(cols[-1], collapse= "+"), " -1"))
      varFixedTest <- varIdent(form = formFxdTest)
      # lme.fit <- lmeFit(formLMMTest, random = list(ExpIndTest = pdIdent(formRdmTest)),
      #                   weights = varFixedTest, data = data[-indEnv, ], nrestarts = 5)
      lme.fit <- lmeFit(formLMMTest, random = list(ExpIndTest = pdIdent(formRdmTest)),
                        weights = NULL, data = data[-indEnv, ], nrestarts = 5)
      # confInt_beta <- intervals(lme.fit, which = "fixed", level = 1- alphaBonfCor/2)$fixed %>%
      #   as.data.frame() %>%
      #   select(c("lower", "upper"))
      confInt_beta_lmm <- intervals(lme.fit, which = "fixed", level = 1- alphaBonfCor/2)$fixed %>%
        as.data.frame() %>%
        select(c("lower", "upper"))
      # tauHat <- VarCorr(lme.fit)[1, "StdDev", drop = FALSE] %>% as.double()
      rm(ExpIndTest, envir = .GlobalEnv)
    } else{
      stop(paste0("Test ", test, " is currently not implemented."))
    }
    
    # construct confidence intervals for '(beta + b) - b'
    # confInt_betaPlusBMinusB <- confInt_beta
    # confInt_betaPlusBMinusB[, "lower"] <- confInt_betaPlusB[, "lower"] # - qnorm(1-alphaBonfCor/4)*tauHat
    # confInt_betaPlusBMinusB[, "upper"] <- confInt_betaPlusB[, "upper"] # + qnorm(1-alphaBonfCor/4)*tauHat
    
    # check if 'confInt_beta' and 'confInt_betaPlusBMinusB' overlaps
    # checkOverlap <- function(i) {
    #   confInt_betaPlusBMinusB[i, "upper"] < confInt_beta[i, "lower"] || 
    #     confInt_beta[i, "upper"] < confInt_betaPlusBMinusB[i, "lower"]
    # }
    checkOverlap <- function(i) {
      confInt_beta_lm[i, "upper"] < confInt_beta_lmm[i, "lower"] ||
        confInt_beta_lmm[i, "upper"] < confInt_beta_lm[i, "lower"]
    }
    # checkOverlapOut <- sapply(1:nrow(confInt_beta), checkOverlap) %>% sum()
    checkOverlapOut <- sapply(1:nrow(confInt_beta_lm), checkOverlap) %>% sum()
    if(checkOverlapOut > 0) {notRejected <- FALSE}
    i <- i + 1
  }
  result <- if(notRejected) {"not rejected"} else{"rejected"}

  return(list(res = result, coeff = coeff, coeffSE = coeffSE, confInt = confInt))
}
