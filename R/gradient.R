#' Gradient by the numerical approximation
#'
#' \code{gradient} is a function based on the algorithm of numerical differentiation for estimating  the gradient
#' using finite difference approximations.
#'
#' @param x A point (number or vector).
#' @param f A objective function.
#' @param fx A number, the value of function in x.
#' @param eps A number, for calculating the gradient.
#' @return g the numerical approximation of the gradient.
#'
#' @examples
#' f <- functios(x) {... return x}
#' x1 <- c(2,2)
#' gradient(x1, f)
#' gradiente(x1, f, fx = f(x1), eps = 1e-3)
#'
#' @references
#' \enumerate{
#' \item Wikipedia, \emph{Numerical differentiation}, \url{https://en.wikipedia.org/wiki/Numerical_differentiation}.
#' }

gradient <- function(x, f, fx = NULL, eps = 1e-4)
{
  g = c()
  #ncf <- 0

  if (is.null(fx))
  {
    fx <- f(x)
  }

  for (i in 1:length(x)){
    e <- rep(0,length(x))
    e[i] <- 1
    g[i] <- (f(x + (eps * e)) - fx) / eps

    #ncf <- ncf + 1
  }

  return (g)
}
