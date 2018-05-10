#' Gradient by the numerical approximation
#'
#' \code{gradient} is a function based on the algorithm of numerical differentiation for estimating the gradient
#' using finite difference approximations.
#'
#' @param f A function, representing the objective function.
#' @param x A number or vector with length n, indicating the current point.
#' @param fx A number, the objective value calculated at x.
#' @param dx A number, the small perturbation in x.
#' @param method character string, specifying the discretization method.
#' "ffd": Forward finite differences. Requires n+1 evaluations (or n if fx is provided)
#' "cfd": Central (or Symmetrical) finite differences. Requires 2n evaluations (fx is not used)
#'
#' @return df the numerical approximation of the gradient.
#'
#' @examples
#' f <- function(x) {return ( sum(x^2) )}
#' df_analytical <- function (x) {return (2*x)}
#' x <- c(1,1) #current point
#' gradient(f, x) #uses ffd by default
#' gradient(f, x, method = "cfd")
#' df_analytical(x)
#'
#' x <- seq(1,100)
#' dfx_analytical <- df_analytical(x)
#' dfx_ffd <- gradient(f, x, method = "ffd")
#' dfx_cfd <- gradient(f, x, method = "cfd")
#' norm(dfx_analytical - dfx_ffd, type = "2") #error in the FFD approximation
#' norm(dfx_analytical - dfx_cfd, type = "2") #error in the CFD approximation
#'
#' @references
#' \enumerate{
#' \item Wikipedia, \emph{Numerical differentiation}, \url{https://en.wikipedia.org/wiki/Numerical_differentiation}.
#' }

gradient <- function(f, x, fx = NULL, dx = 1e-6, method = c("ffd", "cfd"))
{
  method <- match.arg(method)
  if (method == "ffd")
  {
    df <- ffd(x, f, fx, dx)
  }
  else
  {
    df <- cfd(x, f, dx)
  }

  return (df)
}

## Sub-functions
# Finite forward differences
ffd <- function(x, f, fx = NULL, dx = 1e-6)
{
  n = length(x) #number of variables
  df = numeric(n) #initialize the gradient

  if (is.null(fx))
  {
    fx <- f(x)
  }

  for (i in 1:n){
    e <- numeric(n)
    e[i] <- 1
    df[i] <- (f(x + (dx * e)) - fx) / dx
  }

  return (df)
}

# Central finite differences
cfd <- function(x, f, dx = 1e-6)
{
  n = length(x) #number of variables
  df = numeric(n) #initialize the gradient

  for (i in 1:n){
    e <- numeric(n)
    e[i] <- 1
    df[i] <- (f(x + dx*e) - f(x - dx*e)) / (2*dx)
  }

  return (df)
}
