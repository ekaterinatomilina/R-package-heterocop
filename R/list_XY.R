# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

gaussCopule <- function(R, n){
  d = dim(R)[1]
  A = as.matrix(eigen(R)$vectors%*%diag(sqrt(eigen(R)$values),2,2))
  data = matrix(0, n, d)
  for (i in c(1:n)){
    z = as.matrix(rnorm(d, mean = 0, sd = 1))
    x = pnorm(A%*%z)
    data[i,] = t(x)}
  data = data.frame(data)
  return(data)
}


#' Data generator
#'
#' @param R Matrix, 2x2 correlation matrix
#' @param n Integer, length of the simulated data
#' @param list_qlaws List, list of two quasi-inverse function ( must be chosen among (qbeta, qbinom, qchisq, qubif, qexp, qf, qgamma, qgeom, qhyper, qnbinom, qnorm, qpois, qt, qweibull))
#' @param list_param List, list with lists of parameters
#'
#' @return list of two vectors of length n such that the joint distribution function of those vectors is given  by F(.,.) = C(F_1(.), F_2(.)) where C is the gaussian Copula related to R and F_1, F_2 are respectively the distribution functions of the margins.
#' @export
#'
#' @examples list_XY(matrix(c(1,0.5,0.5,1),2,2), 100, list(loi1 = qexp, loi2 = qnorm), list(3,c(0,1)))
list_XY <- function(R, n, list_qlaws, list_param){
  XY = gaussCopule(R, n)
  dd = 2
  liste_Xs = c()
  for (i in 1:dd){
    if (length(list_param[[i]]) == 2){
      liste_Xs[[i]] = list_qlaws[[i]](XY[,i], list_param[[i]][1], list_param[[i]][2])
    }
    else if (length(list_param[[i]]) == 3){
      liste_Xs[[i]] = list_qlaws[[i]](XY[,i], list_param[[i]][1], list_param[[i]][2], list_param[[i]][3])
    }
    else {liste_Xs[[i]] = list_qlaws[[i]](XY[,i], list_param[[i]])}
  }
  return(liste_Xs)
}


