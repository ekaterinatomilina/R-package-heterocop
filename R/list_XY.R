# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' gaussCopule
#'
#' @description gaussian copula data simulator based on the generalized inverse method
#'
#' @param R correlation matrix of the copula
#' @param n vnumber of observations
#' @return matrix of d vectors (in columns) of length n such that the joint distribution function of those vectors is given  by F(.,.) = C(F_1(.),..., F_d(.)) where C is the gaussian Copula related to a block diagonal matrix and F_1,.., F_d are respectively the distribution functions of the margins.
#' @examples list_XY_d(100, c(3,1,4), c(0.2,0.1,0.5), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(5,4,4), random = TRUE)
#' @export
#'

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






