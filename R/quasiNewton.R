quasiNewton <- function(obj.list, x.list, maxNI = 50, eps.df = 1-6, eps = 1e-4, alpha0 = 1, rho = 0.5, c = 1e-4, methodgradient = c("ffd", "cfd")){

  #Copy of values
  obj <- obj.list$f

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

  #--------Quasi-Newton Method--------------#

  hesI <- diag(length(x.list$x))

  for (k in 1:maxNI){

    dfx.k <- dfx.k1
    x.k <- x.k1

    searchD <- as.numeric(-hesI %*% dfx.k)

    #Backtracking
    alpha <- backtracking(obj, x.k, searchD, alpha0, rho, c)
    xp <- x.k$x + as.numeric(searchD*alpha)
    x.k1 <- list(x = xp,
                 fx = obj(xp),
                 dfx = dFun(obj, xp, method = methodgradient))

    #Gradient in the new point
    dfx.k1 <- x.k1$dfx

    #BFGS constants
    v <- x.k$x - x.k1$x
    r <- dfx.k - dfx.k1

    C <- bfgs(r, v, hesI)

    #Hessian new value
    hesI <- hesI + C

    if (identical(stopping_condition(x.k1, x.k), TRUE)) #default eps.df or eps.f?
    {
      break
    }
  }

  return(x.k1)
}
