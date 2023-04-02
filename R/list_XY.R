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






