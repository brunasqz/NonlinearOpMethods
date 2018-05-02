#' Stopping conditions in iterative methods
#'
#' \code{stopping_conditions} is a function used in iterative methods for stopping the process
#' when the conditions are fulfilled. These conditions are stabilization of the function
#' value, reduction of the gradient norm to a value less than epsilon and maximum number of
#' evaluation of objective funtion.
#'
#' @param x_k1 A list that contains the information of the point in next iteration and must
#' have the following names, x (for the point), fx (for the value of the function in x), dfx
#' (for the value of the gradiente in x), as the example \cr
#' \code{x_example <- list(x = c(1,1), fx = 12, dfx = c(-1, 2))}. \cr
#' @param x_k A list that contains the information of the point in current iterantion, in the
#' same pattern as \code{x_k1}.
#' @param Options.list A list as the example \cr
#' \code{options <- list(eps_f = , eps_df =)} \cr
#'
#' @return Returns a results of boolean algebra.
#'
#' @section Options.list:
#' In this function, Options.list may contain the following names:
#' \itemize{
#' \item Eps_f, for minimum variation of the objective function in the iterations.
#' \item Eps_df, for minimum variation of the gradient in the iterations.
#' }
#'
#' @seealso The documentation of function \code{checkparameters} for the parameter Options.list.
#'
#' @examples
#' options <- list(Eps_f = 1e-9, Eps_df = 1e-8)
#' x1 <- list(x = c(1,1), fx = 12, dfx = c(-1,2))
#' xk <- list(x = c(0.008, 1.01), fx = 11.00009, dfx = c(-0.9999, 2.9999))
#' stopping_condition(x1, xk, options)


stopping_condition <- function(x_k1, x_k, Options.list){ #eps_f eps_df in Options.list

  cond1 <- (abs(x_k1$`F(x)` - x_k$`F(x)`)) <= (Options.list$Eps_f)
  cond2 <- (normE(x_k1$`dF(x)`)) <= (Options.list$Eps_df)
  # cond3 <- n_f > nfmax

  return (cond1 | cond2 )
}
