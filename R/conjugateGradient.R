#' Conjugate Gradient Method.
#'
#' \code{conjugateGradient} is a method that uses the search direction as opposed to the gradient at a point x
#' of the iteration. If the user does not pass in \code{obj.list} the desired gradient function, as default the
#' numerical calculation by the finite differences is used to approximate the gradient. It is important to remind
#' that calculating the gradient will be accurate only for linear functions.
#'
#' @param obj.list A list that must have the names, functionObj (for the objective function to be optimized),
#' gradientObj (for the corresponding gradient function), as the example \cr
#' \code{ objFG <- list(functionObj = f, gradientObj = df)} \cr \code{f} and \code{df} returns values.
#' @param X.list A list that must have the names, X (for the point), F(x) (for the value of the function in x),
#' dF(x) (for the value of the gradiente in x), as the example \cr
#' \code{x_example <- list(X = c(1,1), `F(x)` = 12, `dF(x)` = c(-1, 2))}. \cr
#' @param Options.list A list as the example \cr
#' \code{options <- list(A = , B = , Eps = ,  Method = , RhoBacktracking = , cBacktracing = , use_bracketing = ,
#' step_bracketing =, maxNI =)} \cr
#' @param eps A number.
#' @return Returns a X.list in the point optimum
#'
#' @section Options.list:
#' In this function, Options.list may contain the following names:
#' \itemize{
#' \item use_bracketing = TRUE or FALSE, for determine whether bracketing should be used when choosing the
#' interval for search (for linesearch).
#' \item step_bracketing, for determine the step length in bracketing; if use_bracketing = TRUE, step_bracketing
#' must have a value (for linesearch).
#' \item A e B (numbers), user-determined search interval; default A = 0 and B = 1 (for linesearch).
#' \item Method (character), in format "quadraticinterpolation" or "goldensection", if its NULL or not
#' passed, by default, the method used is backtracking (for linesearch function).
#' \item rhoBacktracking and cBacktracking (numbers), control variables in backtracking (for linesearch
#' function).
#' \item Eps, epsilon = error tolerance, default = 1e-6
#' \item maxNI, maximum number for iterations.
#' \item Eps_f, for minimum variation of the objective function in the iterations.
#' \item Eps_df, for minimum variation of the gradient in the iterations.
#' }
#'
#' @seealso The documentation of the functions \code{checkparameters} and \code{linesearch}.
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
#' conjugateGradient(objfunctions, x, c(-3,2), options1)
#' conjugateGradient(objfunctions, x, c(-3,2), options2)
#' conjugateGradient(objfunctions, x, c(-3,2), options4)
#'
#'
#' @references
#' \enumerate{
#' \item Ramirez, Jaime A.; Campelo, Felipe; Guimaraes, Frederico G.; Batista, Lucas S.; Takahashi, Ricardo H. C.;
#' \emph{Notas de aula de Otimizacao}, pages 6:9.
#' \item Wikipedia, \emph{Conjugate Gradient Method}, \url{https://en.wikipedia.org/wiki/Conjugate_gradient_method}.
#' }
#' @export

conjugateGradient <- function(obj.list, x.list, Options.list, eps = 1e-4) {
  #Checking parameters
  outcheck <- checkparameters(obj.list, x.list, Options.list)
  obj.list <- outcheck[[1]]
  x.list <- outcheck[[2]]
  Options.list <- outcheck[[3]]

  #Copy of values
  obj <- obj.list$functionObj
  dFun <- obj.list$gradientObj

  dfx_k1 <- x.list$`dF(x)`
  x_k1 <- x.list

  #---------Conjugate Gradient----------------#
  searchD_k1 <- -dfx_k1


  for(k in 1:(Options.list$maxNI)) {
    dfx_k <- dfx_k1
    x_k <- x_k1
    searchD_k <- searchD_k1

    #Linesearch
    out <- linesearch(obj.list, x_k, searchD_k, Options.list)
    x_k1 <- out[[1]]
    alpha <- out[[2]]

    dfx_k1 <- dFun(x_k1$X, obj)
    searchD_k1 <- -dfx_k1 + as.numeric(abs(searchD_k1) ^2 / abs(searchD_k)^2) * searchD_k

    if (identical(stopping_condition(x_k1, x_k, Options.list), TRUE))
    {
      break
    }
  }
  return(x_k1)
}
