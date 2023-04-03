#' @title lmeFit.R
#'
#' @description This function helps to deal with errors when fitting a linear
#'    mixed-effects model using the package 'nlme'. During testing a lot of fatal
#'    errors occurred when fitting a lme-object
#'
#' @param form formula: Specifies the formula for the fixed-effects. 
#' @param random list: Each element corresponds to a hierarchical level for LMMs. 
#' @param weights varFunc: A variance function for the noise term as in glms.
#' @param data data.frame: Contains any variables that are specified in 'form'
#'    or in 'random'. 
#' 
#' @return lme.fit lme: A fitted lme-object. 

lmeFit <- function(form, random, weights = NULL, data, nrestarts = 5) {
  ctrl <- lmeControl(maxIter = 100, msMaxIter = 100,
                     niterEM = 25, tolerance = 1e-6,
                     msTol = 1e-7, opt = "opt", optimMethod = "BFGS")
  withRestarts(
    tryCatch(
      {
        if(!is.null(weights)) {
          lme.fit <- lme(form, random = random, weights = weights, data = data, 
                         method = 'ML', control = ctrl) %>%
            suppressWarnings()
        } else{
          lme.fit <- lme(form, random = random, data = data, 
                         method = 'ML', control = ctrl) %>%
            suppressWarnings()
        }
      }, error = function(e) {
        if(!exists("errorLog")) {
          errorLog <<- paste0(e, "\n")
        } else{
          errorLog <<- c(errorLog, paste0(e, "\n"))
        }
        nrestarts <<- nrestarts - 1
        invokeRestart("rerun")
      }
    ), 
    rerun = function() {
      message("trying again")
      if(nrestarts > 0) {
        lmeFit(form, random, weights, data, nrestarts = nrestarts)
      } else{
        error <- errorLog
        rm("errorLog", envir = .GlobalEnv)
        stop(error)
      }
    }
  )
  return(lme.fit)
}
