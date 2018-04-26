#' Univariate function
#'
#' \code{univariate_f} is a function in therms of \code{alpha} and used in optimization methods.
#' Usually represented by the greek letter "phi" has the following relation with the objective
#' function \eqn{\phi(\alpha) = f(x + S*\alpha)}, where S represents the search direction.
#'
#' @param fx A objective function.
#' @param x A point.
#' @param searchD A search direction.
#' @return Returns the step length \code{alpha} for the univariate function.
#'
#' @examples
#' f <- function(x) {... return(x)}
#' phi <- function(f, (1,2), (-3,2))
#' phi(alpha)
#' @export



univariate_f <- function(fx, x, searchD) {
  f <- fx(x)
  function(alpha) {
    f <- fx(x + alpha * searchD)

    return(f)
  }
}
