#' Constrained optimization with Gradient Projection
#'
#' \code{gradientproj} solves a multivariate function with constraints using the
#' Gradient Projection.
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
#' @param constraint A list, with the following names \cr
#' xmin: lower restriction \cr
#' xmax: upper restraint \cr
#' @param maxNI maximum number of iterations
#' @param eps tolerance for stop codition
#' @param alpha0 Initial step size in the backtracking
#' @param rho A constant to reduce alpha in backtracking
#' @param c A small constant, control parameter.
#' @param ... parameters used in the check_parameters function
#'
#' @return Returns a list with the (approximate) optimum.
#'
#' @examples
#' f <- function(x) {
#'  f1 <- x[1]^2 + 2*x[2]^2 + 3*x[3]^2 - 2*x[1] -4*x[2] -6*x[3] + 6
#'  return(f1)
#' }
#' x0 <- list(x = c(2.8,3,3.5), fx = f(c(2.8,3,3.5)))
#' const <- list(xmin = c(2,2,2), xmax = c(4,4,4))
#' gradientproj(f, x0, const)
#'
#' #can also be used to list:
#' obj <- list(functionObj = f)
#' gradientproj(obj, x0, const)
#'
#' @references
#' \enumerate{
#' \item Nocedal, C.T.; \emph{Iterative Methods for optimization}.
#' }
#' @export

#Reference: Kelley - Iterative Methods for optimization
gradientproj <- function(obj.list, x.list, constraint, maxNI = 50, eps = 1e-4, alpha0 = 1, c = 1e-4, rho = 0.5, ...) {
  ## Initial considerations
  # Check input
  out <- check_parameters(obj.list, x.list, ...)
  obj.list <- out[[1]]
  x.list <- out[[2]]

  if(any(x.list$x > constraint$xmax) | any(x.list$x < constraint$xmin))
    stop("Error: x point out of the restrictions")

  x.k1 <- x.list
  dfx.k1 <- x.list$dfx
  searchD <- alpha0*dfx.k1

  for(i in 1:maxNI) {
    x.k <- x.k1

    #Backtracking
    alpha <- backtracking.gp(obj.list, x.k, constraint, alpha0, rho, c)
    searchD <- -dfx.k1

    x.k1$x <- projection(x.k$x + alpha*searchD, constraint)
    x.k1 <- evaluate(obj.list, x.k1$x)

    #Gradient in the new point
    dfx.k1 <- x.k1$dfx

    if (identical(stopping_condition(x.k1, x.k), TRUE)) #default eps.df or eps.f?
    {
      break
    }
  }
  return(x.k1)

}

projection <- function(x, constraint) {
  u <- constraint$xmax
  l <- constraint$xmin
  out1 <- pmin(x,u)
  out <- pmax(out1, l)

  return(out)
}

backtracking.gp <- function(obj.list, x.list, constraint, alpha0 = 1, rho = 0.5, c = 1e-4)
{
  x <- x.list$x
  fx <- x.list$fx
  dfx <- x.list$dfx

  alpha <- alpha0 #initial alpha
  sdfx <- sum(-dfx * dfx) #dot product between dfx and searchD
  x.proj <- projection(x - alpha*dfx, constraint)

  while( fx - obj.list$f(x.proj) < c*alpha*sdfx )
  {
    alpha <- alpha*rho
    x.proj <- projection(x - alpha*dfx, constraint)
  }

  return(alpha)
}


