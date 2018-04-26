#' Euclidean norm.
#'
#' @param vec A vector.
#' @return Returns the Euclidean norm of \code{vec}.
#' @examples
#' normE(c(-1,2))
#' normE(c(10, 25))

normE <- function(vec){
  normV <- sqrt(sum(vec^2));
  return (normV);
}
