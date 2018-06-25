#' Shortcut to evaluate a new point at an objective function and its gradient
#'
#' @param obj.list A list with the problem functions
#' @param x A vector solution to be evaluated
#'
#' @return A list with x evaluated at \code{f} and \code{df}

evaluate <- function(obj.list, x)
{
  # Evaluate a point x at f and df and return an appropriate list
  x.list <- list("x" = x)
  x.list$fx <- obj.list$f(x)
  x.list$dfx <- obj.list$df(x, x.list$fx)

  return (x.list)
}
