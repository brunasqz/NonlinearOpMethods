# #Função 1 - pag 61 livro Jan Snyman
#
g1 <- function(x)
{
  gx <- c()
  gx[1] <- 1 - x[1]
  gx[2] <- -x[2]

  return(gx)
}

f1 <- function(x) {
  return((1/3)*(x[1] + 1)^3 + x[2])
}

obj1 <- list(functionObj = f1,
             f = f1,
             g = g1)

x1 <- list(x = c(2, 1),
           fx = f1(c(2, 1)))


#Exemplo 2 -> tp campelo 2017/2
f2 <- function(x) {
  return((x[1] - 5)^2 + (x[2] - 3)^2)
}

g2 <- function(x) {
  g2x <- c()

  g2x[1] <- x[1] + x[2] -3
  g2x[2] <- -x[1] + 2*x[2] -4

  return(g2x)
}



obj2 <- list(functionObj = f2,
             f = f2,
             g = g2)

x2 <- list(x = c(0, 0),
           fx = f2(c(0, 0)))


#Exemplo 3 -> tp campelo 2017/2

f3 <- function(x){
  fx <- x[1]^4 - 2*(x[1]^2)*x[2] + x[1]^2 +x[1]*(x[2]^2) - 2*x[1] +4
  return (fx)
}
h3 <- function(x){

  hfx <- x[1]^2 + x[2]^2 -2

  return(hfx);
}
g3 <- function(x){
  g <- 0.25*x[1]^2 + 0.75*x[2]^2 - 1;

  return(g);
}

obj3 <- list(functionObj = f3,
             f = f3,
             g = g3,
             hf = h3)

x3 <- list(x = c(2, -4),
           fx = f2(c(2, -4)))

#Functions Test for quasinewton and conjugateGradient
#Function 1 x* = [1,1,1] f(x*) = 0
# f1 <- function(x) {
#   f1 <- x[1]^2 + 2*x[2]^2 + 3*x[3]^2 - 2*x[1] -4*x[2] -6*x[3] + 6
#
#   return(f1)
# }
#
# x1.list <- list(x = c(3,3,3), fx = f1(c(3,3,3)))
#
# obj1 <- list(functionObj = f1)
#
# #Function 2 x* = [1,1], f1(x*) = 0
# f2 <- function(x) {
#   f2 <- x[1]^4 - 2*(x[1]^2)*x[2] + x[1]^2 + x[2]^2 - 2*x[1] + 1
#
#   return(f2)
# }
#
# x2.list <- list(x = c(3,3), fx = f2(c(3,3)))
# obj2 <- list(functionObj = f2)
#
# #Function3 x* = [2,1], f(x*) = 0
# f3 <- function(x) {
#   f3 <- x[1]^4 - 8*x[1]^3 + 25*(x[1]^2) + 4*(x[2]^2) -4*x[1]*x[2] -32*x[1] + 16
#   return(f3)
# }
#
# x3.list <- list(x = c(3,3), fx = f3(c(3,3)))
# obj3 <- list(functionObj = f3)
#
# #Function 4 x* = [1,1] f(x*) = 0
# f4 <- function(x) {
#   f4 <- 100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
#
#   return(f4)
# }
#
# x4.list <- list(x = c(1,1), fx = f4(c(1,1)))
# obj4 <- list(functionObj = f4)
#
# #Function 5 - x* = [0.57085, -0.93955591, 0.7681755]
# f5 <- function(x) {
#   f5 <- x[1]^4 + x[1]^3 - x[1] + x[2]^4 - x[2]^2 + x[2] + x[3]^2 - x[3] + x[1]*x[2]*x[3]
#
#   return(f5)
#
# }
#
# x5.list <- list(x = c(1,-1,1), fx = f5(c(1,-1,1)))
#
# obj5 <- list(functionObj = f5)
