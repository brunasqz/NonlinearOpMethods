#Reference: Kelley - Iterative Methods for optimization
gradientproj <- function(obj.list, x.list, constraint, maxNI = 50, eps = 1e-4, alpha0 = 1, c = 1e-4, rho = 0.5, ...) {
  ## Initial considerations
  # Check input
  out <- check_parameters(obj.list, x.list, ...)
  obj.list <- out[[1]]
  x.list <- out[[2]]

  x.k1 <- x.list
  dfx.k1 <- x.list$dfx

  for(i in 1:maxNI) {
    dfx.k <- dfx.k1
    x.k <- x.k1

     searchD <- projection(x.k$x - dfx.k, constraint) - x.k$x
    # x.proj <- projection(x.k$x, constraint)
    # x.listproj <- evaluate(obj.list, x.proj)
    # searchD <- -x.listproj$dfx

    #Backtracking
    alpha <- backtracking.gp(obj.list, x.k, constraint, searchD, alpha0, rho, c)
    x.k1$x <- x.k$x + alpha*searchD
    x.k1 <- evaluate(obj.list, x.k1$x)

    #Gradient in the new point
    dfx.k1 <- x.k1$dfx

    if (identical(stopping_condition(x.k1, x.k), TRUE)) #default eps.df or eps.f?
    {
      break
    }
  }
  return(x.k1)

}

projection <- function(x, constraint) {
  u <- constraint$xmax
  l <- constraint$xmin
  out1 <- pmin(x,u)
  out <- pmax(out1, l)

  return(out)
}

backtracking.gp <- function(obj.list, x.list, constraint, searchD, alpha0 = 1, rho = 0.5, c = 1e-4)
{
  x <- x.list$x
  fx <- x.list$fx
  dfx <- x.list$dfx

  alpha <- alpha0 #initial alpha
  sdfx <- sum(dfx * searchD) #dot product between dfx and searchD
  x.proj <- projection(x - dfx, constraint)

  while( obj.list$f(x) - obj.list$f(x.proj) > fx + c*alpha*sdfx )
  {
    alpha <- alpha*rho
    x.proj <- projection(x - dfx, constraint)
  }

  return(alpha)
}

