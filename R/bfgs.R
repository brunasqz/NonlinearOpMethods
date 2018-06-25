#' Broyden Fletcher Goldfarb Shanno update method.
#'
#' \code{bfgs} is an update scheme to get a new Hessian inverse approximation in
#' a Quasi-Newton method.
#'
#' @param x.list A list with the current solution
#' @param xnew.list A list with the new solution in the same iteration
#' @param H Current Hessian inverse approximation
#'
#' @seealso The documentation of the Quasi-newton method \code{quasiNewton} in
#' this package.
#'
#' @references
#' \enumerate{
#' \item Wikipedia, \emph{Quasi-newton method} \url{https://en.wikipedia.org/wiki/Quasi-Newton_method}.
#' \item Ramirez, Jaime A.; Campelo, Felipe; Guimaraes, Frederico G.; Batista, Lucas S.; Takahashi, Ricardo H. C.;
#' \emph{Notas de aula de Otimizacao}, pages 26:28.
#' }

bfgs <- function(x.list, xnew.list, H)
{
  dx <- xnew.list$x - x.list$x
  y <- xnew.list$dfx - x.list$dfx
  I = diag(length(dx))
  # Split the terms for easier visualization
  term1 <- I - (dx %*% t(y)) / as.numeric( (t(y) %*% dx) )
  term2 <- I - (y %*% t(dx)) / as.numeric( (t(y) %*% dx) )
  term3 <- (dx %*% t(dx)) / as.numeric( (t(y) %*% dx) )
  return ( term1 %*% H %*% term2 + term3 )
}

# bfgs <- function(r, v, hess) {
#
#   I <- diag(length(v))
#   v <- matrix(v, ncol = 1) #xk - xk1
#   r <- matrix(r, ncol = 1) #dfxk - dfxk1
#   tv <- t(v)
#   tr <- t(r)
#   const1 <- v %*% tr / as.numeric(tr %*% v)
#   const2 <- r %*% tv / as.numeric(tr %*% v)
#   C <- (I - const1) %*% hess %*% (I - const2) + v %*% tv / as.numeric(tr %*% v)
#
#
#   return(C)
# }
