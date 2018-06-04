SUMT <- function(x.list, obj.list, rho0 = 1, c = 10, eps.f = 1e-4, eps.df = 1e-6, maxNI = 100, ...) {

  out <- checkparameters(obj.list, x.list)
  x.list <- out[[2]]
  obj.list <- out[[1]]

  x.k1 <- x.list
  for(k in 1:maxNI) {

    obj.list$functionObj <- p(obj.list, rho0)
    x.k <- x.k1

    x.k1 <- quasiNewton(obj.list = obj.list, x.list = x.k, maxNI = maxNI, eps = eps.f, ...)

    if (identical(stopping_condition(x.k1, x.k, eps.f, eps.df), TRUE)) #default eps.df or eps.f?
    {
      break
    }

    rho0 <- c*rho0
    message(rho0)
  }

  return(x.k1)

}

p <- function(obj.list, rho)
{
  function (x) {
    Fx <- obj.list$f(x)
    #Validation
    if(!("g") %in% names(obj.list)) {
      Gx <- 0
    } else {
      Gx <- obj.list$g(x)
    }

    if(("hf") %in% names(obj.list)) {
      Hfx <- obj.list$hf(x)
    } else {
      Hfx <- 0
    }

    Gx <- as.matrix(Gx)
    Hfx <- as.matrix(Hfx)

    g <- 0
    for(i in length(Gx)) {
      g <- g + rho*max(Gx[i], 0)^2
    }

    hf <- 0
    for(i in length(Hfx)) {
      hf <- hf + rho*(Hfx[i]^2)
    }

    p <- Fx + g + hf

    #p <- Fx + rho*sum(Hx^2) + rho*sum(max(Gx, 0)^2)
    return(p)
  }
}



