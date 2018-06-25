#' Unconstrained optimization with Quasi-Newton
#'
#' \code{quasiNewton} solves a multivariate function using a Quasi-Newton with
#' BFGS udpate.
#'
#' @param obj.list Either an objective function or a list with the following
#' names \cr
#'   f: objective function
#'   df: gradient of f (it not provided, defaults to the numerical version)
#' @param x.list Either a vector with an initial solution or a list with the
#' following names \cr
#'   x: a vector with its value in the search space \cr
#'   fx: a scalar with its objective value \cr
#'   dfx: a vector with its gradient value \cr
#' @param maxNI maximum number of iterations
#' @param eps.df tolerance in the norm of the gradient
#' @param eps.x tolerance in the norm of the difference between two consecutive
#' solutions
#' @param alpha0 Initial step size in the backtracking
#' @param rho A constant to reduce alpha in backtracking
#' @param c A small constant, control parameter.
#'
#' @return Returns a list with the (approximate) optimum.
#'
#' @section Todo:
#' Maybe include alpha0, rho and c into ... and get the appropriate parameters
#' in each function?
#'
#' @examples
#' # The popular Rosenbrock function
#' f <- function(x)
#' {
#'   x1 <- x[1]
#'   x2 <- x[2]
#'   return ( 100*(x2 - x1^2)^2 + (1 - x1)^2 )
#' }
#' x0 <- c(-1.2,1) #usual starting point
#' # Run the Quasi-Newton with default parameters (very close to c(1,1))
#' x.list <- quasiNewton(f, x0) #x.list$x = c(0.9996996 0.9993989))
#' # Notice the result can be improved by using "CFD" in the numerical gradient
#' # instead of "FFD":
#' x.list <- quasiNewton(f, x0, method = "cfd") #x = c(1,1)
#'
#' @references
#' \enumerate{
#' \item Nocedal, Jorge; Wright, Stephen J.; \emph{Numerical Optimization}, 2nd ed., page 37.
#' \item Wikipedia, \emph{Quasi-newton method} \url{https://en.wikipedia.org/wiki/Quasi-Newton_method}.
#' }
#' @export


quasiNewton <- function(obj.list, x.list,
                        maxNI = 50,
                        eps.df = 1e-6,
                        eps.x = 1e-6,
                        alpha0 = 1,
                        rho = 0.5,
                        c = 1e-4,
                        ...)
{
  ## Initial considerations
  # Check input
  out <- check_parameters(obj.list, x.list, ...)
  obj.list <- out[[1]]
  x.list <- out[[2]]
  hesI <- diag(length(x.list$x)) #initial approximation of the Hessian inverse

  for (k in 1:maxNI)
  {
    searchD <- as.numeric(-hesI %*% x.list$dfx)
    alpha <- backtracking(obj.list, x.list, searchD, alpha0, rho, c)
    xnew <- x.list$x + searchD*alpha
    xnew.list <- evaluate(obj.list, xnew)
    hesI <- bfgs(x.list, xnew.list, hesI)

    ## Stopping condition (maybe convert into another function?)
    test.x = norm( xnew.list$x - x.list$x, type = "2" ) <= eps.x
    test.df = norm( xnew.list$dfx, type = "2") <= eps.df
    if (test.x | test.df)
      break

    x.list <- xnew.list #update x
  }

  return(xnew.list)
}
