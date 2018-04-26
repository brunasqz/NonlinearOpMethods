#' Line search
#'
#' \code{linesearch} is a function that englobes the methods: backtracking, quasi-newton and goldensection.
#' The extra parameters in this functions for control  and use in methos must be passed via Options.list.
#'
#'
#' @param obj.list A list that must have the names, functionObj (for the objective function to be optimized),
#' gradientObj (for the corresponding gradient function), as the example \cr
#' \code{ objFG <- list(functionObj = f, gradientObj = df)} \cr \code{f} and \code{df} returns values.
#' @param X.list A list that must have the names, x (for the point), fx (for the value of the function in x),
#' dfx (for the value of the gradiente in x), as the example \cr
#' \code{x_example <- list(x = c(1,1), fx = 12, dfx = c(-1, 2))}. \cr
#' @param searchD A search direction for the method.
#' @param Options.list A list as the example \cr
#' \code{options <- list(A = , B = , Eps = ,  Method = , RhoBacktracking = , cBacktracing = , use_bracketing = ,
#' step_bracketing =)} \cr
#'
#' @return Returns a list that contains a \code{X.list} in the point optimum and the step length \code{alpha}
#' for the univariate function
#'
#' @section Options.list:
#' In this function, Options.list may contain the following names:
#' \itemize{
#' \item use_bracketing = TRUE or FALSE, for determine whether bracketing should be used when choosing the
#' interval for search.
#' \item step_bracketing, for determine the step length in bracketing; if use_bracketing = TRUE, step_bracketing
#' must have a value.
#' \item A e B (numbers), user-determined search interval; default A = 0 and B = 1.
#' \item Method (a character string), in format "quadraticinterpolation" or "goldensection", if its NULL or not
#' passed, by default, the method used is backtracking.
#' \item rhoBacktracking and cBacktracking (numbers), control variables in backtracking.
#' \item Eps, epsilon = error tolerance, default = 1e-6
#' }
#'
#' @seealso The documentations of fuctions \code{backtracking}, \code{bracketing}, \code{quadraticinterpolation}
#' and \code{goldensection} in this package.
#'
#' @examples
#' options1 <-  list(Method = "goldensection", step_bracketing = 0.5, use_bracketing = TRUE)
#' options2 <- list(use_bracketing = FALSE, rhoBacktrackting = 0.6, cBacktrackting = 1e-6, Eps = 1e-7)
#' options3 <- list(Method = "quadraticinterpolation", A = 3, B = 5, use_bracketing = FALSE)
#'
#' f <- function{return(fx)}
#' dfun <- functuib{return(dfx)}
#' x <- list(x = c(1,2), fx = f(c(1,2)), dfx = dfun(c(1,2)))
#' objfunctions <- list(functionObj = f, gradientObj = dfun)
#'
#' linesearch(objfunctions, x, c(-3,2), options1)
#' linesearch(objfunctions, x, c(-3,2), options2)
#' linesearch(objfunctions, x, c(-3,2), options4)
#'
#'
#' @export


linesearch <- function (obj.list,
                        X.list,
                        searchD,
                        Options.list)
{
  #Checking the parameters
  checkparameters(obj.list, X.list, Options.list)

  #Copy of values
  obj <- obj.list$functionObj
  x <- X.list[[1]]
  eps <- Options.list$Eps

  #
  phi_fk <- univariate_f(obj, x, searchD)

  if((!("Method") %in% names(Options.list)) | is.null(Options.list$Method)) {

    rhoD <- Options.list$RhoBacktracking
    cB <- Options.list$cBacktracking

    #Backtracking - linear function default
    alphak <- backtracking(obj, X.list, searchD, rhoD, cB)

  } else {

    if((!("use_bracketing") %in% names(Options.list)) & (!("step_bracketing") %in% names(Options.list))) {
      a <- 0; b <- 1
    } else {

      if(identical(Options.list$use_bracketing, FALSE)) {
        a <- 0; b <- 1
      } else {
        step_b <- Options.list$step_bracketing
        results <- bracketing(phi_fk, step_b)

        a <- results[[1]]
        b <- results[[2]]
      }
    }

    # ----- Use the parameters entered by the to calculate alpha -----#

    if(Options.list$Method == "goldensection") {

      alphak <- goldensection(phi_fk, a, b, eps)

    } else if(Options.list$Method == "quadraticinterpolation") {

      alphak <- quadraticinterpolation(phi_fk, a, b, eps)

    }
  }

  #finds x
  x <- x + as.numeric(alphak*searchD)
  fx <- obj(x)
  dfx <- gradient(x, obj, fx)
  x_optimum <- list(x, fx, dfx)

  return(list(x_optimum, alphak))
}

