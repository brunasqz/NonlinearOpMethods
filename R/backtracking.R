#' Backtracking line search.
#'
#' \code{backtracking} is a line search method to determine a step length that minimize the function
#' in a given search direction. The method starts with a \code{alpha = 1} and iteratively modify it.
#' Based on the book "Numerical Optimization", Jorge Nocedal, for ensure that the algorithm makes reasonable progress along
#' in the given search direction, is used the control parameters with a default value, \code{rho = 0.5}
#' and \code{c = 1e-4}.
#'
#' @param obj A objective function.
#' @param X.list A list that must have the names, X (for the point), F(x) (for the value of the function in x),
#' dF(x) (for the value of the gradiente in x), as the example \cr
#' \code{x_example <- list(X = c(1,1), `F(x)` = 12, `dF(x)` = c(-1, 2))}. \cr
#' @param searchD A search direction for the method.
#' @param rho A number, control parameter.
#' @param c A number, control parameter.
#' @return Returns the step length \code{alpha} for the univariate function.
#'
#'
#' @examples
#' f <- function(x) {... return(x)}
#' backtracking(f, x_example, c(1,-2))
#' backtracking(f, x_example, searchD = c(1,-2), rho = 0.1, c = 1e-6)
#'
#' @references
#' \enumerate{
#' \item Nocedal, Jorge; Wright, Stephen J.; \emph{Numerical Optimization}, 2nd ed., page 37.
#' \item Wikipedia, \emph{Backtracking line search} \url{https://en.wikipedia.org/wiki/Backtracking_line_search}.
#' }
#' @export



backtracking <- function(obj, X.list, searchD, rho = 0.5, c = 1e-4)
{
  alpha <- 1
  x <- X.list$X
  fx <- X.list$`F(x)`
  dfx <- X.list$`dF(x)`


  k <- alpha * searchD

  while( obj(x + as.numeric(k)) > fx + (c*alpha)%*%(t(dfx))%*%searchD)
  {

    alpha <- alpha*rho
    k <- alpha*searchD #alpha*searchD

    #rho <- rho*(0.1)
    #if(rho < 1e-6)
    #{
    #  rho <- rhoD
    #}
  }
  message("Method: Backtracking")
  return(alpha)

}
