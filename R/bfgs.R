#' Broyden Fletcher Goldfarb Shanno method.
#'
#' \code{bfgs} is a iterative method for update the approximate hessian (or its inverse) in Quasi-Newton method.
#'
#' @param r A vector, the difference of the gradient values at the points (dfx) of the current and previous iteration.
#' @param v A vector, the difference of the points (x) of the current and previous iteration.
#' @param hess A matrix, hessian.
#' @return C, a approximation of the inverse hessian.
#'
#' @seealso The documentation of the Quasi-newton method \code{quasiNewton} in this package.
#'
#' @examples
#' r <- df_k - df_k1
#' v <- x_k - x_k1
#' hess <- diag(length(x_k))
#' hess <- hess + bfgs(r, v, hess)
#'
#' @references
#' \enumerate{
#' \item Wikipedia, \emph{Quasi-newton method} \url{https://en.wikipedia.org/wiki/Quasi-Newton_method}.
#' \item Ramirez, Jaime A.; Campelo, Felipe; Guimaraes, Frederico G.; Batista, Lucas S.; Takahashi, Ricardo H. C.;
#' \emph{Notas de aula de Otimizacao}, pages 26:28.
#' }

bfgs <- function(r, v, hess) {

  I <- diag(length(v))
  v <- matrix(v, ncol = 1) #xk - xk1
  r <- matrix(r, ncol = 1) #dfxk - dfxk1
  tv <- t(v)
  tr <- t(r)
  const1 <- v %*% tr / as.numeric(tr %*% v)
  const2 <- r %*% tv / as.numeric(tr %*% v)
  C <- (I - const1) %*% hess %*% (I - const2) + v %*% tv / as.numeric(tr %*% v)


  return(C)
}
