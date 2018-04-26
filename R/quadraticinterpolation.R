#' Quadratic Interpolation Method.
#'
#' \code{quadraticinterpolation} is a one-dimensional minimization method that finds the minimizing
#' step length \code{alpha} through the approximation by a quadratic function. If the \code{alpha}
#' is not sot sufficiently to the minimum a refit of the function is used.
#'
#' @section Reduction of points, interval:
#' To reduce the interval between the three points, first evaluate the value of the function
#' in alpha and in the poinst: lower (\code{l}), upper (\code{u}) and middle (\code{m}).
#' After this evaluation, it is analyzed (through comparisons) at which two points the
#' alpha is between and are chosen the next values for the refite
#'
#' @param phif A univariate function.
#' @param l,u Numbers, initial interval (or points).
#' @param eps A number, used in the stop condition.
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
#' quadraticinterpolation(phi, a = 0, b = 10, eps = 1e-4)
#' quadraticinterpolation(phi, 3, 7, 1e-5)
#' @references
#' \enumerate{
#' \item Rao, Singiresu S.; \emph{Engineering Optimization Theory and and Practice}, 4th ed., pages 273:279.
#'
#' }
#' @export


quadraticinterpolation <- function(phi_f, l = 0, u = 1, eps = 1e-6)
{
  #polynomials evaluated at 3 points(l, m, u) for interpolation: phi_fl, phi_fm, phi_fu
  if(l == u)
  {
    l <- 0
  }
  m <- (l + u)/2
  phi_fl <- phi_f(l)
  phi_fm <- phi_f(m)
  phi_fu <- phi_f(u)


  #a, b and c are coefficients of the polynomial (hx)
  y <- ((l - m)*(m - u)*(u - l))
  a <- (phi_fl*m*u*(u - m) + phi_fm*u*l*(l - u) + phi_fu*l*m*(m - l))/y
  b <- (phi_fl*(m^2 - u^2) + phi_fm*(u^2 - l^2) + phi_fu*(l^2 - m^2))/y
  c <- -(phi_fl*(m - u) + phi_fm*(u - l) + phi_fu*(l - m))/y

  alpha <- -b/(2*c)
  falpha <- phi_f(alpha)

  ncf <- 4

  hx <- function(x)  c*x^2 + b*x + a #environment global

  while((abs((hx(alpha) - falpha)) > eps))
  {
    #choice of l, m and C new values
    if(alpha <= m)
    {
      if(phi_fm >= falpha){
        u <- m; phi_fu <- phi_fm
        m <- alpha; phi_fm <- falpha
      }else{
        l <- alpha; phi_fl <- falpha
      }
    }else{
      if(phi_fm >= falpha){
        l <- m; phi_fl <- phi_fm
        m <- alpha; phi_fm <- falpha
      }else{
        u <- alpha; phi_fu <- falpha
      }

    }


    #refitting the function
    y <- ((l - m)*(m - u)*(u - l))
    a <- (phi_fl*m*u*(u - m) + phi_fm*u*l*(l - u) + phi_fu*l*m*(m - l))/y
    b <- (phi_fl*(m^2 - u^2) + phi_fm*(u^2 - l^2) + phi_fu*(l^2 - m^2))/y
    c <- -(phi_fl*(m - u) + phi_fm*(u - l) + phi_fu*(l - m))/y

    #quadratic interpolation

    alpha <- -b/(2*c)
    falpha <- phi_f(alpha)
    ncf <- ncf + 1

  }
  message("Method: Quadratic Interpolation. Number of calls of the objective function: ", ncf)
  return(alpha)
}
