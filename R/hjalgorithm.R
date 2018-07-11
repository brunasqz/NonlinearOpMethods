hjalgorithm <- function(x, obj, h, eps) {
  searchD <- rbind(diag(1, nrow = length(x)), diag(-1, nrow = length(x)))
  while(h > eps) {
    x.k <- x

    for(i in 1:nrow(searchD)) {
      x.k1 <- x
      x.k1 <- x.k1 + h * searchD[i,]

      if(obj(x.k1) < obj(x)) {
        x <- x.k1
      }
    }

    if(identical(x.k, x)) {
      h <- h / 2
    }
  }
  x
}





#x base point
#h pattern size -> h_k
#search directions (vj) is the jth column of a direction matrix V
#Kelley -> let V = I be the matrix of coordinate directions
# hjalgorithm <- function(x, obj, h, searchD) {
#   for(k in 1:length(h))
#   {
#     xb <- x
#     xc <- x
#     flag <- 1 #used to signal failure and trigger a shrink step
#     out <- hjexplore(x, xc, obj, h[k], searchD)
#     xb <- out[[1]]
#     xc <- out[[2]]
#     flag <- out[[3]]
#
#     while(identical(flag, 1)) {
#       d <- x - xb
#       xb <- x
#       xc <- x + d
#
#       out <- hjexplore(x, xc, obj, h[k], searchD)
#       xb <- out[[1]]
#       xc <- out[[2]]
#       flag <- out[[3]]
#
#       if(identical(flag, 0)) {
#         xc <- x
#         out <- hjexplore(x, xc, obj, h[k], searchD)
#         xb <- out[[1]]
#         xc <- out[[2]]
#         flag <- out[[3]]
#       }
#     }
#   }
# }
#
# hjexplore <- function(xb, xc, obj, h, searchD, flag) {
#   fb <- obj(xb)
#   d <- 0
#   flag <- 0
#   xcb <- xb
#   fcb <- fb
#   xt <- xc
#
#   for(i in 1:nrow(searchD)) {
#     p <- xt + h * searchD[i,]
#
#     fp <- obj(p)
#     if(fp > fb) {
#       p <- xt - h * searchD[i,]
#     }else if (fp < fb) {
#       xt <- xcb
#       p <- xt
#       fb <- f(xcb)
#     }
#
#     if(!identical(xcb,xb)) {
#       flag = 1
#       xb <- xcb
#     }
#   }
#
#   return(list(xb, xc, flag))
#
# }
