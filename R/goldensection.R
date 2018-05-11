#' Golden Section Method.
#'
#' \code{goldensection} is a line search method that computes the optimum of a
#' univariate function over an interval \code{[a,b]} where it is known to be
#' unimodal.
#'
#' @param phi A univariate function.
#' @param a,b Numbers, initial interval where the function is unimodal.
#' @param eps A number, a tolerance to stop the optimization when the current
#' interval becomes smaller \code{eps}.
#' @return Step length \code{alpha} where phi attains its minimum in
#' \code{[a,b]} (assuming unimodality).
#'
#' @section About the parameter phi:
#' \code{phi} should be either an originally univariate function
#' \code{phi(alpha)}, or the result of a multi-variate one over a fixed interval
#' \code{phi(alpha) := f(x + alpha*searchD)}, which should come from the
#' \code{univariate_f} function.
#'
#' @seealso The documentation of the function \code{univariate_f} in this package.
#'
#' @examples
#' # A univariate function
#' phi <- function(alpha)
#' {
#'   return (alpha^5 - 5*alpha^3 - 20*alpha + 5)
#' }
#' # Plot the function for convenience over the interval [-5,5]
#' plot(phi, -5, 5)
#' # Compute the optimum in the interval [0,5] where it is unimodal
#' alpha.opt1 <- goldensection(phi, a = 0, b = 5) #approx 2
#' # If it is changed to [-5,0] where it is not unimodal, an extreme is returned
#' alpha.opt2 <- goldensection(phi, a = -5, b = 0) #return -5
#'
#' # Considering a multi-variate function, assume given the current point x and
#' # the search direction searchD
#' f <- function (x)
#' {
#'   i <- seq(1, length(x))
#'   return (sum( (x - i)^2 ))
#' }
#' x <- c(0,0,0,0) #current point
#' searchD <- c(1,2,3,4) #search direction
#' phi <- univariate_f(f, x, searchD) #return univariate version
#' # Compute optimum over the line x + alpha*searchD in [0,2]
#' alpha.opt <- goldensection(phi, a = 0, b = 2) #approx 1
#'
#' @references
#' \enumerate{
#' \item Rao, Singiresu S.; \emph{Engineering Optimization Theory and and Practice}, 4th ed., page 267.
#' \item Ramirez, Jaime A.; Campelo, Felipe; Guimaraes, Frederico G.; Batista, Lucas S.; Takahashi, Ricardo H. C.;
#' \emph{Notas de aula de Otimizacao}, pages 10:13.
#' }
#' @export

goldensection <- function(phi, a = 0, b = 1, eps.alpha = 1e-6)
{
  gr <- (sqrt(5) - 1)/2 #golden ratio

  # Internal values of the interval
  alpha.a <- b - gr * (b - a)
  alpha.b <- a + gr * (b - a)
  phi.a <- phi(alpha.a)
  phi.b <- phi(alpha.b)
  ncf <- 2 #number of function evaluations

  while ((b - a) > eps.alpha)
  {
    if(phi.b < phi.a){
      a <- alpha.a
      alpha.a <- alpha.b
      alpha.b <- a + gr * (b - a)
      phi.a <- phi.b
      phi.b <- phi(alpha.b)
      ncf <- ncf + 1

    }else{
      b <- alpha.b
      alpha.b <- alpha.a
      alpha.a <- b - gr * (b - a)
      phi.b <- phi.a
      phi.a <- phi(alpha.a)
      ncf <- ncf + 1
    }
  }
  alpha.opt <-(a + b)/2

  message("Method: Golden Section. Number of calls of the objective function: ", ncf)
  return(alpha.opt)
}

