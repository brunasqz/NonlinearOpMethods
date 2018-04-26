#' Bracketing.
#'
#' \code{bracketing} is a function that searches for some restricted range that can have the
#' minimum point between the two values of a unimodal function. For this implementation,
#' the initial point have a default value (\code{alpha = 0}) and a fixed step size for move
#' in a favorable direction. In addition, the algorithm was based on the search with accelerated
#' step size, that consists of doubling the step size as long as the move results in an improvement
#' of the objective function.
#'
#'
#' @param phif A univariate function.
#' @param sstep A number, the sstep defines whether the search will be in a positive or negative direction.
#' @param alpha1 A number, the initial point.
#' @return Returns a range in the format of a list.
#'
#' @section About the parameter phif:
#' \code{phif} should be a univariate function that returns the value of the function
#' at a given point.
#'
#' @seealso The documentation of the function \code{univariate_f} in this package.
#'
#' @examples
#' phi <- function(alpha) {...}
#' bracketing(phi)
#' bracketing(phi, -0.5)
#' bracketing(phi, sstep = -0.5, alpha1 = 3)
#'
#' @references
#' \enumerate{
#' \item Rao, Singiresu S.; \emph{Engineering Optimization Theory and and Practice}, 4th ed., pages 273:279.
#'
#' }
#' @export


bracketing <- function (phi_f, sstep = 0.05, alpha1 = 0)
{
  alpha2 <- alpha1 + sstep
  alphai <- alpha2

  phi_f1 <- phi_f(0)
  phi_f2 <- phi_f(alpha2)


  while(phi_f2 < phi_f1) {
    sstep <- 2 * sstep
    alphai <- alpha2
    alpha2 <- alpha1 + sstep

    phi_f1 <- phi_f(alphai)
    phi_f2 <- phi_f(alpha2)
  }

  return(list(alphai, alpha2))

}
