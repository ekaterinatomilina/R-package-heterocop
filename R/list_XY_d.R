VecEstim <- function(Corr){
  rho = c()
  dim = as.numeric(dim(Corr)[1])
  for (i in 1:(dim-1)){
    rho = c(rho, Corr[i, (i+1):dim])
  }
  return(rho)
}

gaussCopule_ <- function(vec,d,n){
  norm.cop <- copula::normalCopula(vec, dim = d, dispstr = "un")
  XY_g = data.frame(copula::rCopula(n, norm.cop))
  return(as.matrix(XY_g))
}

matrice_diag_blocs <- function(blocs, coeff){
  d = sum(blocs)
  n = length(blocs)
  R = matrix(0,d,d)
  R[1:blocs[1],1:blocs[1]] = coeff[1]
  compteur = 2
  for (j in 2:n){
    if (blocs[j]==1){R[j,j]=1}
    else{
      R[(sum(blocs[1:j-1])+1):sum(blocs[1:j]), (sum(blocs[1:j-1])+1):sum(blocs[1:j])] = coeff[compteur]
      compteur = compteur+1
    }
  }
  for (i in 1:d){R[i,i] = 1}
  return(R)
}


#' CopulaSim
#'
#' @description multidimensional data simulator based on the generalized inverse method
#'
#' @param n Integer, length of the simulated data
#' @param blocks vector, vector with the size of the blocks on the diagonal of the matrix
#' @param coeff vector, coefficient associated to each block (there is no need to inform the coefficient for one dimensional block (see examples))
#' @param list_qlaws  List, list of  quasi-inverse functions (must be chosen among (qbeta, qbinom, qchisq, qubif, qexp, qf, qgamma, qgeom, qhyper, qnbinom, qnorm, qpois, qt, qweibull)) with their specified parameters
#' @param repetition vector, number of repetition of each laws
#' @param random if TRUE the distribution of each laws will be randomly distributed among the blocks. if FALSE the order enter as a parameter will be preserved. (by default : random = TRUE)
#'
#' @return matrix of d vectors (in columns) of length n such that the joint distribution function of those vectors is given  by F(.,.) = C(F_1(.),..., F_d(.)) where C is the gaussian Copula related to a block diagonal matrix and F_1,.., F_d are respectively the distribution functions of the margins.
#' @examples list_XY_d(100, c(3,1,4), c(0.2,0.1,0.5), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(5,4,4), random = TRUE)
#' @export
#'
CopulaSim <- function(n, blocks, coeff, list_qlaws, repetition, random = TRUE){

  R = matrice_diag_blocs(blocks, coeff)
  d = sum(blocks)
  XY = gaussCopule(R,n)
  
  nb = c()
  for (l in 1:length(repetition))
    nb = append(nb, rep(l,repetition[l]))
  if (random == TRUE){
    nb = sample(nb)
  }
  
  liste_Xs = matrix(0,n,d)
  for (j in 1:d){
    i = nb[j]
    liste_Xs[,j] = list_qlaws[[i]](XY[, i])
  }
  return(liste_Xs)
}

#'



