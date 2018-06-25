#Reference: Kelley - Iterative Methods for optimization
gradientproj <- function(obj.list, x.list, constraint, maxNI = 50, eps = 1e-4, methodgradient = c("ffd", "cfd")) {
  #Copy of values
  obj <- obj.list$f

  projection <- function(constraint, obj) {
    out <- pmax(x,constraint$l)
    out <- pmin(out, constraint$u)


    return(obj(out))
  }

  if(!("df") %in% names(obj.list)) {
    dFun <- gradient
  } else {
    dFun <- obj.list$df
  }
  if(!("dfx") %in% names(x.list)) {
    dfx.k1 <- dFun(obj, x.list$x, method = methodgradient)
  } else {
    dfx.k1 <- x.list$dfx
  }

  x.k1 <- x.list

  for(i in 1:maxNI) {
    dfx.k <- dfx.k1
    x.k <- x.k1

    searchD <- dfx.k

    #Backtracking
    alpha <- backtracking.gp(projection, x.k, searchD, alpha0, rho, c)
    xp <- x.k$x + as.numeric(searchD*alpha)
    x.k1 <- list(x = xp,
                 fx = obj(xp),
                 dfx = dFun(obj, xp, method = methodgradient))

    #Gradient in the new point
    dfx.k1 <- x.k1$dfx

    if (identical(stopping_condition(x.k1, x.k), TRUE)) #default eps.df or eps.f?
    {
      break
    }
  }


}

backtracking.gp <- function(f, x.list, searchD, alpha0 = 1, rho = 0.5, c = 1e-4)
{
  x <- x.list$x
  fx <- x.list$fx
  dfx <- x.list$dfx

  alpha <- alpha0 #initial alpha
  sdfx <- sum(dfx * searchD) #dot product between dfx and searchD

  while( f(x + alpha*searchD) > fx + c*alpha*sdfx )
  {
    alpha <- alpha*rho
  }

  return(alpha)
}


projection <- function(x, l, u) {
  out <- pmax(x,l)
  out <- pmin(out, u)
  return(out)
}
