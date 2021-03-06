#' Euclidean norm.
#'
#' @param vec A vector.
#' @return Returns the Euclidean norm of \code{vec}.

normE <- function(vec){
  normV <- sqrt(sum(vec^2));
  return (normV);
}
