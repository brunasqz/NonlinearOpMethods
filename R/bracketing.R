#' Bracketing.
#'
#' \code{bracketing} is a function that determines a range (bracket)
#' \code{[a, b]} of values where a univariate function has possibly a local
#' minimum. This function is able to implement an accelerated search.
#'
#' @param phi A univariate function.
#' @param alpha0 A number, starting point.
#' @param phi0 A number, objective value of \code{alpha0}.
#' @param sstep A number with the step to increase alpha0 until a bracket is
#' determined. It can accept positive or negative values, see below.
#' @param t Multiplicator factor to increase sstep in each iteration. If
#' \code{t > 1}, an accelerated search is performed. Use \code{t = 1} for the
#' fixed step version.
#' @param kmax Maximum number of iterations. It can be used to determine
#' degenerate problems with no finite minimum, when the bracket does not exist.
#'
#' @return Returns a bracket in the format of a list. These intervals can be
#' used in other interval halving methods, such as the goldensection.
#'
#' @details
#' Given a starting point \code{alpha0}, the bracketing process follows: \cr
#' 1. Set \code{alpha1 <- alpha0} and \code{alpha2 <- alpha1 + sstep} \cr
#' 2. While the function is decreasing in this direction, keep increasing alpha1
#' and alpha2, that is, Repeat the following until
#' \code{phi(alpha2) > phi(alpha1)}: \cr
#'  2.1 Set sstep <- t*sstep \cr
#'  2.2 Set alpha0 <- alpha1 \cr
#'  2.3 Set alpha1 <- alpha2 \cr
#'  2.3 Set alpha2 <- alpha1 + sstep \cr
#' In the end, either \code{[alpha1, alpha2]} or \code{[alpha0, alpha2]} can be
#' chosen as bracket. The latter was preferred here for safety purposes. Also,
#' notice that the method here can employ an accelerated search by setting
#' \code{t > 1}, in which the step of alpha increases in each iteration.
#'
#' @section About the parameter sstep:
#' If \code{phi} comes from \code{phi(alpha) := f(x + alpha*s)} for a given
#' point \code{x} and a \emph{descent} direction \code{s}, then \code{sstep > 0}
#' will correctly yield a descent direction in the univariate version. For
#' general single-variable functions, the option \code{sstep < 0} is left to
#' handle specific situations, as shown in the examples.
#'
#' @seealso The documentation of the functions \code{univariate_f} and
#' \code{goldensection} in this package.
#'
#' @examples
#' # A univariate function
#' phi <- function(alpha)
#' {
#'   return (alpha^5 - 5*alpha^3 - 20*alpha + 5)
#' }
#' # Plot the function for convenience over the interval [-5,5]
#' plot(phi, -5, 5)
#' # In the inverval [0,5], it has a minimum in alpha = 2. Starting from
#' # alpha0 = 0, for instance, the following bracket is obtained
#' bracket <- bracketing(phi, alpha0 = 0) #(0.75, 3.15), containing 2
#' # If phi0 = phi(alpha0) = 5 is provided, the results are the same, but the
#' # number of evaluations is reduced by one
#' bracket <- bracketing(phi, alpha0 = 0, phi0 = 2)
#' # Setting t = 1 (fixed step size) generates a smaller interval, but increases
#' # the number of evaluations considerably
#' bracket <- bracketing(phi, alpha0 = 0, t = 1) #(1.95, 2.05)
#' # Starting from alpha0 = 4, only one iteration is executed because phi is
#' # increasing in this direction (the optimum is actually at 4)
#' bracket <- bracketing(phi, alpha0 = 4, sstep = 0.05) #(4, 4.05)
#' # However, changing the sign of sstep yields a correct bracketing of the
#' # true minimum.
#' bracket <- bracketing(phi, alpha0 = 4, sstep = -0.05) #(0.85, 3.25)
#' # Finally, from -4 to the left this function has no finite minimum, and thus
#' # the algorithm goes until the total iterations are run. This is an indicator
#' # of a problem with no minimum.
#' bracket <- bracketing(phi, alpha0 = -4, sstep = -0.05)
#'
#' @references
#' \enumerate{
#' \item Rao, Singiresu S.; \emph{Engineering Optimization Theory and and Practice}, 4th ed., pages 273:279.
#'
#' }
#' @export


bracketing <- function(phi, alpha0 = 0, phi0 = NULL, sstep = 0.05, t = 2, maxNI = 50)
{
  ncf <- 0 #number of function evaluations

  # Handle starting point. If objective value not provided, compute it here
  alpha1 <- alpha0
  phi1 <- phi0
  if (is.null(phi1))
  {
    phi1 <- phi(alpha1)
    ncf <- ncf + 1
  }

  alpha2 <- alpha1 + sstep
  phi2 <- phi(alpha2)
  ncf <- ncf + 1

  # Use a for loop to limit the number of evaluations, specially if the function
  # an infinite minimum
  for (k in 1:maxNI)
  {
    if (phi2 > phi1)
    {
      break
    }

    sstep <- t*sstep #increase the step in alpha

    alpha0 <- alpha1
    phi0 <- phi1

    alpha1 <- alpha2
    phi1 <- phi2

    alpha2 <- alpha1 + sstep
    phi2 <- phi(alpha2)
    ncf <- ncf + 1
  }

  # In the end, return (alpha0, alpha2) as bracketing. Notice that the correct
  # order is (alpha0, alpha2) if sstep is positive, otherwise, it should return
  # (alpha2, alpha0)

  if (sstep > 0)
  {
    bracket <- list(alpha0, alpha2)
  }
  else
  {
    bracket <- list(alpha2, alpha0)
  }

  message("Bracketing. Number of function evaluations: ", ncf)
  return(bracket)

}
