% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hjalgorithm.R
\name{hjalgorithm}
\alias{hjalgorithm}
\title{Hooke and Jeeves algorithm (Pattern search)}
\usage{
hjalgorithm(obj, x, h = 0.25, eps = 1e-06, constraint = NULL, c1 = 1.1,
  c2 = 0.5)
}
\arguments{
\item{obj}{A objective function.}

\item{x}{A number or vector with length n, indicating the current point.}

\item{h}{A number, search step size}

\item{eps}{A number, tolerance for h.}

\item{constraint}{A list, with the following names \cr
xmin: lower restriction \cr
xmax: upper restraint \cr}

\item{c1}{A constant to update the search step size}

\item{c2}{A constant to update the search step size}
}
\value{
Returns a list with the (approximate) optimum and the number of objective
function evaluations.
}
\description{
\code{hjalgorithm} is a function based on the algorithm of Hooke and Jeeves
for numerical optimization of functions without the use of the gradient.
}
\examples{
# The popular Rosenbrock function
f <- function(x)
{
  x1 <- x[1]
  x2 <- x[2]
  return ( 100*(x2 - x1^2)^2 + (1 - x1)^2 )
}
x0 <- c(-1.2,1) #usual starting point
hjalgorithm(f, x0)
#With constraint
const <- list(xmin = c(2, 2), xmax = c(4,4))
x0 <- c(3,3)
hjalgorithm(f, c(3,3), constraint = const)

}
\references{
\enumerate{
\item Wikipedia, \emph{Pattern search (optimization)}, \url{https://en.wikipedia.org/wiki/Pattern_search_(optimization)}.
}
}
