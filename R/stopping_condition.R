#' Stopping conditions in iterative methods
#'
#' \code{stopping_conditions} is a function used in iterative methods for stopping the process
#' when the conditions are fulfilled. These conditions are stabilization of the function
#' value, reduction of the gradient norm to a value less than epsilon and maximum number of
#' evaluation of objective funtion.
#'
#' @param x.k1 A list that contains the information of the point in current iterantion.
#'  It must have the names \cr
#'   x: a vector with its value in the search space \cr
#'   fx: a scalar with its objective value \cr
#'   dfx: a vector with its gradient value \cr
#' @param x.k A list that contains the information of the point in current iterantion, in the
#' same pattern as \code{x_k1}.
#' @param eps.f A number for minimum variation of the objective function in the iterations.
#' @param eps.df A number for minimum variation of the gradient in the iterations.
#' @return Returns a results of boolean algebra.
#'
#' @examples
#' x1 <- list(x = c(1,1), fx = 12, dfx = c(-1,2))
#' xk <- list(x = c(0.008, 1.01), fx = 11.00009, dfx = c(-0.9999, 2.9999))
#' stopping_condition(x1, xk, eps_f = 1e-9, eps_df = 1e-8)


stopping_condition <- function(x.k1, x.k, eps.f = 1e-9, eps.df = 1e-8){ #eps_f eps_df in Options.list

  cond1 <- (abs(x.k1$fx - x.k$fx)) <= (eps.f)
  cond2 <- (normE(x.k1)) <= (eps.df)
  # cond3 <- n_f > nfmax

  return (cond1 | cond2 )
}
