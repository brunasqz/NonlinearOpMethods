#' Golden Section Method.
#'
#' \code{goldensection} is a line search method that reduces the initial interval
#' to get a optimum with the highest precision needed. In this implementation the reduction
#' of the interval, when it happens, occurs by dividing by half.
#'
#' @param phif A univariate function.
#' @param a,b Numbers, initial interval.
#' @param eps A number.
#' @return Returns the step length \code{alpha}.
#'
#' @section About the parameter phif:
#' \code{phif} should be a univariate function that returns the value of the function
#' at a given point.
#'
#' @seealso The documentation of the function \code{univariate_f} in this package.
#'
#' @examples
#' phi <- function(alpha) {...}
#' goldensection(phi, a = 0, b = 10, eps = 1e-4)
#' goldensection(phi, 3, 7, 1e-5)
#' @references
#' \enumerate{
#' \item Rao, Singiresu S.; \emph{Engineering Optimization Theory and and Practice}, 4th ed., page 267.
#' \item Ramirez, Jaime A.; Campelo, Felipe; Guimaraes, Frederico G.; Batista, Lucas S.; Takahashi, Ricardo H. C.;
#' \emph{Notas de aula de Otimizacao}, pages 10:13.
#' }
#' @export

goldensection <- function(phif, a = 0, b = 1, eps = 1e-6)
{
  #phif -> phi function, depends alpha
  #a, b -> intervals

  phi <- (sqrt(5) - 1)/2

  alphaa <- b - phi * (b - a)
  alphab <- a + phi * (b - a)

  #implementation of the theta functions
  thetaA <- phif(alphaa)
  thetaB <- phif(alphab)

  ncf <- 2

  while ((b - a) > eps)
  {
    #print(c(a, xa, xb, b, fx(a), thetaA, thetaB, fx(b)))

    if(thetaB < thetaA){
      a <- alphaa
      alphaa <- alphab
      alphab <- a + phi * (b - a)
      thetaA <- thetaB
      thetaB <- phif(alphab)

      ncf <- ncf + 1
    }else{
      b <-alphab
      alphab <- alphaa
      alphaa <- b - phi * (b - a)
      thetaB <- thetaA;
      thetaA <- phif(alphaa)

      ncf <- ncf + 1
    }
  }
  alpha <-(a + b)/2

  message("Method: Golden Section. Number of calls of the objective function: ", ncf)
  return(alpha)
}

