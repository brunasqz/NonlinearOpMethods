#' Check parameter (Options.list)
#'
#' \code{checkparameters} is a function that checks if the parameter Options.list is complete
#' to use in some functions in this package.
#'
#' @section About the parameter Options.list:
#' This list contains all parameters that is not fundamental for the function. In other words,
#' Options.list receives instances that have default values in all functions that are used,
#' but the user can have a greater control of the operation of the code, if desired.
#' \itemize{
#' \item use_bracketing = TRUE or FALSE, used in function \code{linesearch}
#' \item step_bracketing, (number), used in function \code{linesearch}
#' \item A and B (numbers), used in function \code{linesearch}
#' \item Method ("quadraticinterpolation" or "goldensection"), used in function \code{linesearch}
#' \item rhoBacktracking e cBacktracking, used in function \code{linesearch}
#' \item Eps, epsilon, error tolerance, default = 1e-6
#' \item Eps_f for the functions \code{stopping_conditions}, \code{quasiNewton} and \code{conjugateGradient}
#' \item Eps_df for the functions \code{stopping_conditions}, \code{quasiNewton} and \code{conjugateGradient}
#' \item maxNI maximum number iterations for the functions \code{stopping_conditions}, \code{quasiNewton} and
#' \code{conjugateGradient}
#' }
#'
#' @param Options.list A list
#' @return Returns a complete Options.list with default values.



# Options.list <- list( A = NULL, B = NULL, maxNI, Eps, Eps_df, Eps_f, Method, RhoBacktracking, cBacktracing, use_bracketing, step_bracketing)



checkparameters <- function(obj.list, x.list, Options.list) {
  #------Checking Obj.list-----------------#

  if(!("gradientFunc") %in% names(obj.list)) {

   dFun <- gradient

    #Complete the list
   obj.list$gradientObj <- dFun
  } else {

    dFun <- obj.list$gradientObj
  }

  #------Checking x.list-----------------#

  if(!("dF(x)") %in% x.list) {

    dfx_k1 <- dFun(x.list$X, obj.list$functionObj)

    #Complete the list
    x.list$`dF(x)` <- dfx_k1
  }

  # ----- Checking parameters of Options.list ----- #

  #maximum number of iterations (Quasi-Newton, Conjugate Gradient)
  if(!("maxNI") %in% names(Options.list)) {

    Options.list$maxNI <- 200
  } else {
    #verify if maxNI is correct
    if(is.numeric(Options.list$maxNI) == FALSE) {
      message("Parameter maxNI in Options.list isn't numeric.")
    }

  }


  if(!("RhoBacktracking") %in% names(Options.list)) {

    Options.list$RhoBacktracking <- 0.5
  } else {
    if(is.numeric(Options.list$RhoBacktracking) == FALSE) {
      message("Parameter RhoBacktracking in Options.list isn't numeric.")
    }
  }

  if(!("cBacktracking") %in% names(Options.list)) {

    Options.list$cBacktracking <- 1e-4
  } else {
    if(is.numeric(Options.list$cBacktracking) == FALSE) {
      message("Parameter cBacktracking in Options.list isn't numeric.")
    }
  }

  if(!("Eps") %in% names(Options.list)) {

    Options.list$Eps <- 1e-6
  }

  if(!("Eps_f") %in% names(Options.list)) {

    Options.list$Eps_f <- 1e-6
  }

  if(!("Eps_df") %in% names(Options.list)) {

    Options.list$Eps_df <- 1e-8
  }

  return(list(obj.list, x.list, Options.list))

}


