#' Check input parameters on quasi-Newton
#'
#' The quasiNewton algorithm can be run by providing different types of input,
#' but they all need to be accomodated to a common format. This is performed
#' in this function
#'
#' @param obj.list Either a function or a list with the following names \cr
#'   f: the objective function \cr
#'   df: its gradient (optional; where the numerical version is used by default)
#'   \cr
#' @param x.list Either a vector with initial solution or a list with the
#' following names \cr
#'   x: a vector with its value in the search space \cr
#'   fx: a scalar with its objective value \cr
#'   dfx: a vector with its gradient value \cr
#' @param ... Optional parameters passed to the \code{gradient} function
#'
#' @return Returns two lists: \cr
#'   \code{obj.list}: a list with names \code{f} and \code{df} with the
#'   objective function and its gradient, respectively \cr
#'   \code{x.list}: a list with the names \code{x}, \code{fx} and \code{dfx}
#'   with the initial point properly evaluated.


check_parameters <- function(obj.list, x.list, ...)
{
  ## Verify obj.list
  # Check if obj.list is a list of functions or simply a function
  if (!(is.list(obj.list)))
  {
    obj.list <- list("f" = obj.list)
  }

  # Check if gradient provided
  if(!("df") %in% names(obj.list))
  {
    obj.list$df <- function (x, fx) gradient(obj.list$f, x, fx = fx, ...)
  }

  ## Verify x.list
  # If x.list is simply a vector, convert it into a list
  if (!(is.list(x.list)))
  {
    x.list <- list("x" = x.list)
  }

  # Evaluate x at the objective function if required
  if(!("fx") %in% names(x.list))
  {
    x.list$fx <- obj.list$f(x.list$x)
  }

  # Evaluate x at the gradient if required
  if(!("dfx") %in% names(x.list))
  {
    x.list$dfx <- obj.list$df(x.list$x, x.list$fx)
  }

  return (list(obj.list, x.list))
}

