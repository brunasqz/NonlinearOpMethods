#' Hooke and Jeeves algorithm (Pattern search)
#'
#' \code{hjalgorithm} is a function based on the algorithm of Hooke and Jeeves
#' for numerical optimization of functions without the use of the gradient.
#'
#' @param obj A objective function.
#' @param x A number or vector with length n, indicating the current point.
#' @param h A number, search step size
#' @param eps A number, tolerance for h.
#' @param constraint A list, with the following names \cr
#' xmin: lower restriction \cr
#' xmax: upper restraint \cr
#' @param c1
#' @param c2
#'
#' @return Returns a list with the (approximate) optimum and the number of objective
#' function evaluations.
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
#' hjalgorithm(x0, f)
#' #With constraint
#' const <- list(xmin = c(2, 2), xmax = c(4,4))
#' x0 <- c(3,3)
#' hjalgorithm(rosenblock, c(3,3), constraint = const)
#'
#' @references
#' \enumerate{
#' \item Wikipedia, \emph{Pattern search (optimization)}, \url{https://en.wikipedia.org/wiki/Pattern_search_(optimization)}.
#' }
#' @export

hjalgorithm <- function(obj, x, h = 0.25, eps = 1e-6, constraint = NULL, c1 = 1.1, c2 = 0.5) {
  # Check input
  if(is.null(constraint)) {
    constraint$xmax <- Inf
    constraint$xmin <- -Inf
  }

  #Initialization of variables
  nfe <- 0
  searchD <- rbind(diag(1, nrow = length(x)), diag(-1, nrow = length(x)))

  while(h > eps) {
    x.k <- x

    fx.k <- obj(x.k)
    nfe <- nfe + 1
    for(i in 1:nrow(searchD)) {
      x.k1 <- x
      x.k1 <- x.k1 + h * searchD[i,]
      x.k1 <- projection(x.k1, constraint)

      fx.k1 <- obj(x.k1)
      nfe <- nfe + 1

      if(fx.k1 < fx.k) {
        x <- x.k1
        fx.k <- fx.k1
        h <- h * c1
      }
    }

    if(identical(x.k, x)) {
      h <- h * c2
    }
  }
  out <- list(xopt = x, fx = obj(x), NumberObjEvaluations = nfe)

  return(out)
}

projection <- function(x, constraint) {
  u <- constraint$xmax
  l <- constraint$xmin
  out1 <- pmin(x,u)
  out <- pmax(out1, l)

  return(out)
}
