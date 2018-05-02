#' Quasi-Newton Method.
#'
#' \code{quasiNewton} is a optimization method that uses the search direction as the result of the multiplication
#' between the inverse of the hessian and the gradient evaluated in the point. For update the approximate hessian
#' is used the bfgs method. If the user does not pass in \code{obj.list} the desired gradient function, as default
#' the numerical calculation by the finite differences is used to approximate the gradient.
#'
#' @param obj.list A list that must have the names, functionObj (for the objective function to be optimized),
#' gradientObj (for the corresponding gradient function), as the example \cr
#' \code{ objFG <- list(functionObj = f, gradientObj = df)} \cr \code{f} and \code{df} returns values.
#' @param X.list A list that must have the names, x (for the point), F(x) (for the value of the function in x),
#' dF(x) (for the value of the gradiente in x), as the example \cr
#' \code{x_example <- list(X = c(1,1), `F(x)` = 12, `dF(x)` = c(-1, 2))}. \cr
#' @param Options.list A list as the example \cr
#' \code{options <- list(A = , B = , Eps = ,  Method = , RhoBacktracking = , cBacktracing = , use_bracketing = ,
#' step_bracketing =, maxNI =)} \cr
#' @return Returns a X.list in the point optimum.
#'
#' @section Options.list:
#' In this function, Options.list may contain the following names:
#' \itemize{
#' \item use_bracketing = TRUE or FALSE, in linesearch for determine whether bracketing should be used when
#' choosing the interval for search.
#' \item step_bracketing, in linesearch for determine the step length in bracketing; if use_bracketing = TRUE,
#' step_bracketing must have a value.
#' \item A e B (numbers), in linesearch, user-determined search interval; default A = 0 and B = 1.
#' \item Method (character), in format "quadraticinterpolation" or "goldensection", if its NULL or not
#' passed, by default, the method used is backtracking (for linesearch function).
#' \item rhoBacktracking and cBacktracking (numbers), control variables in backtracking (for linesearch
#' function).
#' \item Eps, epsilon = error tolerance, default = 1e-6
#' \item Eps_f, for minimum variation of the objective function in the iterations.
#' \item Eps_df, for minimum variation of the gradient in the iterations.
#' \item maxNI, maximum number for iterations.
#' }
#'
#' @examples
#' options1 <-  list(Method = "goldensection", step_bracketing = 0.5, use_bracketing = TRUE, Eps_f = 1e-5, Eps_df = 1e-8)
#' options2 <- list(use_bracketing = FALSE, rhoBacktrackting = 0.6, cBacktrackting = 1e-6, Eps = 1e-7, maxNI = 300)
#' options3 <- list(Method = "quadraticinterpolation", A = 3, B = 5, use_bracketing = FALSE, maxNI = 400, Eps_f = 1e-5, Eps_df = 1e-8)
#'
#' f <- function{return(fx)}
#' dfun <- functuib{return(dfx)}
#' x <- list(X = c(1,2), `F(x)` = f(c(1,2)), `dF(x)` = dfun(c(1,2)))
#' objfunctions <- list(functionObj = f, gradientObj = dfun)
#'
#' quasiNewton(objfunctions, x, c(-3,2), options1)
#' quasiNewton(objfunctions, x, c(-3,2), options2)
#' quasiNewton(objfunctions, x, c(-3,2), options4)
#'
#' @references
#' \enumerate{
#' \item Nocedal, Jorge; Wright, Stephen J.; \emph{Numerical Optimization}, 2nd ed., pages 46:51.
#' \item Rao, Singiresu S.; \emph{Engineering Optimization Theory and and Practice}, 4th ed., page 288:291.
#' \item Ramirez, Jaime A.; Campelo, Felipe; Guimaraes, Frederico G.; Batista, Lucas S.; Takahashi, Ricardo H. C.;
#' \emph{Notas de aula de Otimizacao}, pages 26:28.
#' }
#' @export

quasiNewton <- function(obj.list, x.list, Options.list){
  #Checking the parameters
  outcheck <- checkparameters(obj.list, x.list, Options.list)
  obj.list <- outcheck[[1]]
  x.list <- outcheck[[2]]
  Options.list <- outcheck[[3]]

  #Copy of values
  obj <- obj.list$functionObj
  dFun <- obj.list$gradientObj

  x_k1 <- x.list
  dfx_k1 <- x.list$`dF(x)`

  eps <- Options.list$Eps

  #--------Quasi-Newton Method--------------#

  hesI <- diag(length(x.list$X))

  for (k in 1:(Options.list$maxNI)){

    dfx_k <- dfx_k1
    x_k <- x_k1

    searchD <- as.numeric(-hesI %*% dfx_k)

    #Linesearch
    out <- linesearch(obj.list, x_k, searchD, Options.list)
    x_k1 <- out[[1]] #list
    alpha <- out[[2]] #numeric

    #Gradient in the new point
    dfx_k1 <- dFun(x_k1$X, obj)

    #BFGS constants
    v <- x_k$X - x_k1$X
    r <- dfx_k - dfx_k1

    C <- bfgs(r, v, hesI)

    #Hessian new value
    hesI <- hesI + C

    if (identical(stopping_condition(x_k1, x_k, Options.list), TRUE))
    {
      break
    }
  }
  x_k1

}
