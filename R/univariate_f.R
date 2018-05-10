#' Univariate function
#'
#' In iterative optimization methods, given a point x and a search direction s,
#' a line search consists in finding the (approximate) minimum of a function f
#' in this direction, i.e., we wish to find \code{alpha} such that
#' \eqn{f(x + \alpha*s)} is minimum. Under these conditions, x and s are fixed,
#' and the problem becomes a function \eqn{\phi(\alpha) := f(x + s*\alpha)} of
#' only the step \code{alpha}.
#' \code{univariate_f} returns this unidimensional function for a given f, x and
#' s.
#'
#' @param f The objective function.
#' @param x A number or vector with length n representing the current point.
#' @param searchD A number or vector, representing the search direction.
#' @return Returns an univariate version of the original function
#' \eqn{\phi(\alpha) := f(x + s*\alpha)} of the step length \code{alpha}.
#'
#' @examples
#' ## A sum of translated quadratics
#' f <- function (x)
#' {
#'   i <- seq(1, length(x))
#'   return (sum( (x - i)^2 ))
#' }
#' x <- c(0,0,0,0) #current point
#' searchD <- c(1,2,3,4) #search direction
#' phi <- univariate_f(f, x, searchD)
#' # Compare phi() with f()
#' phi(0) #equivalent to f(x)
#' phi(1) #equivalent to f(x + 1*searchD) = f(c(1,2,3,4)), optimum here
#' @export

univariate_f <- function(f, x, searchD) {
  phi <- f(x)
  function(alpha) {
    phi <- f(x + alpha * searchD)

    return(phi)
  }
}
