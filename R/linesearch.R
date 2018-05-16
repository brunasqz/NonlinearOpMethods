#' Line search
#'
#' \code{linesearch} is a function that englobes the methods: backtracking, quasi-newton and goldensection.
#' The extra parameters in this functions for control  and use in methos must be passed via Options.list.
#'
#'
#' @param obj.list A list that must have the names, functionObj (for the objective function to be optimized),
#' gradientObj (for the corresponding gradient function), as the example \cr
#' \code{ objFG <- list(functionObj = f, gradientObj = df)} \cr \code{f} and \code{df} returns values.
#' @param x.list A list with the current solution. It must have the names \cr
#'   x: a vector with its value in the search space \cr
#'   fx: a scalar with its objective value \cr
#'   dfx: a vector with its gradient value \cr
#' @param searchD A search direction for the method.
#' @param method A character string, in format "quadraticinterpolation" or "goldensection", if its NULL,
#' by default, the method used is backtracking.
#' @param rhoD A number, control variable in backtracking.
#' @param cB A number, control variable in backtracking.
#' @param use.bracketing A boolean, for determine whether bracketing should be used when choosing the
#' interval for search.
#' @param step.bracketing  A numver for determine the step length in bracketing; if use_b = TRUE, step_b
#' must have a value.
#' @param eps, epsilon = error tolerance, default = 1e-6
#'
#' @return Returns a list that contains a \code{x.list} in the point optimum and the step length \code{alpha}
#' for the univariate function
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
#' linesearch(objfunctions, x, c(-3,2), method = "goldensection", step_b = 0.5, use_b= TRUE)
#' linesearch(objfunctions, x, c(-3,2), use_b = FALSE, rhoB = 0.6, cB= 1e-6, eps = 1e-7)
#' linesearch(objfunctions, x, c(-3,2), method = "quadraticinterpolation", use_b = FALSE)
#'
#' @section To do:
#' Modify the documentation, remove comments
#'
#' @export


linesearch <- function (obj.list, x.list, searchD, method = NULL, rhoD = 0.5, cB = 1e-4, use.bracketing = FALSE, step.bracketing = 0.05, eps = 1e-6)
{
  #Checking the parameters
  outcheck <- checkparameters(obj.list, x.list)
  obj.list <- outcheck[[1]]
  x.list <- outcheck[[2]]

  #Copy of values
  obj <- obj.list$functionObj
  x <- x.list$x
  eps <- Options.list$Eps

  #
  phi.fk <- univariate_f(obj, x, searchD)

  if(is.null(method)) {

    #Backtracking - linear function default
    alphak <- backtracking(obj, x.list, searchD, rho = rhoD, c = cB)

  } else {

      if(identical(use.bracketing, FALSE)) {
        a <- 0; b <- 1
      } else {

        results <- bracketing(phi = phi.fk, sstep = step.bracketing)

        a <- results[[1]]
        b <- results[[2]]
      }

    # ----- Use the parameters entered by the to calculate alpha -----#

    if(identical(method,"goldensection")) {

      alphak <- goldensection(phi.fk, a, b, eps)

    } else if(identical(method, "quadraticinterpolation")) {

      alphak <- quadraticinterpolation(phi.fk, a, b, eps)

    }
  }

  #Finds x
  x <- x + as.numeric(alphak*searchD)
  fx <- obj(x)
  dfx <- gradient(x, obj, fx)
  x.optimum <- list(x = x, fx = fx, dfx = dfx)

  return(list(x.optimum, alphak))
}

