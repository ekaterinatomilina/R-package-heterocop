# First, we are going to define the function GCopula that simulates a Gaussian copula by the inversion method

#' GCopula
#'
#' @description gaussian copula data simulator based on the generalized inverse method
#'
#' @param R correlation matrix of the copula
#' @param n number of observations
#' 
#' @return matrix of d vectors (in columns) of length n such that the joint distribution function of those vectors is given  by F(.,.) = C(F_1(.),..., F_d(.)) where C is the gaussian Copula related to a block diagonal matrix and F_1,.., F_d are respectively the distribution functions of the margins.
#' 
#' @examples GCopula(100, c(3,1,4), c(0.2,0.1,0.5), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(5,4,4), random = TRUE)
#' 
#' @export


GCopula <- function(R, n){
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

#an intermediate function to generate a block-wise diagonal matrix
matrice_diag_blocs <- function(blocs, coeff){
  step=c(0,cumsum(blocs))
  d = sum(blocs)
  R = matrix(0,d,d)
  for(i in 1:(length(step)-1)){
    R[((step[i]+1):step[i+1]),((step[i]+1):step[i+1])]=coeff[i]
  }
  for (i in 1:d){R[i,i] = 1}
  return(R)
}

#Now, we are going to use it in order to simulate d vectors linked by a Gaussian copula

#' CopulaSim
#'
#' @description multidimensional gaussian copula data simulator based on the generalized inverse method
#'
#' @param n Integer, length of the simulated data (=number of observations)
#' @param blocks vector, vector with the size of the blocks on the diagonal of the matrix
#' @param coeff vector, coefficient associated to each block (there is no need to inform the coefficient for one dimensional block (see examples))
#' @param list_qlaws  List, list of  quasi-inverse functions (must be chosen among (qbeta, qbinom, qchisq, qubif, qexp, qf, qgamma, qgeom, qhyper, qnbinom, qnorm, qpois, qt, qweibull)) with their specified parameters
#' @param repetition vector, number of repetition of each laws
#' @param random if TRUE the distribution of each laws will be randomly distributed among the blocks. if FALSE the order enter as a parameter will be preserved. (by default : random = TRUE)
#'
#' @return matrix of d vectors (in columns) of length n such that the joint distribution function of those vectors is given  by F(.,.) = C(F_1(.),..., F_d(.)) where C is the gaussian Copula related to a block diagonal matrix and F_1,.., F_d are respectively the distribution functions of the margins.
#' 
#' @examples CopulaSim(100, c(3,1,4), c(0.2,0.1,0.5), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(5,4,4), random = TRUE)
#' 
#' @export

CopulaSim <- function(n, blocks, coeff, list_qlaws, repetition, random = TRUE){
  
  R = matrice_diag_blocs(blocks, coeff)
  d = sum(blocks)
  XY = GCopula(R,n)
  
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

# A function for estimating the correlation matrix


# First, we are going to create several auxiliary functions

#Empirical CDF
nombre_groupe <- function(y_var, list_Y)+return(length(which(list_Y==y_var)))
list_fdr <- function(Z, arg ){
  if(arg == "m"){F <- as.numeric((rank(Z, ties.method = "max")-sapply(Z,nombre_groupe,list_Y = Z))/(length(Z)+1))}
  else {F <- as.numeric(rank(Z, ties.method = "max")/(length(Z)+1))}
}

#This function estimates the CDF values for each variable
list_fdr_d <- function(X, Type){
  F <- array(0,dim=c(dim(X)[1],dim(X)[2],2))
  for (i in 1:dim(X)[2]){
    if (Type[i] == "C"){F[,i,1] = list_fdr(X[,i],1)
    F[,i,2]=c(rep(0,dim(X)[1]))}
    else{
      F[,i,1]=list_fdr(X[,i], "p")
      F[,i,2] = list_fdr(X[,i], "m")
    }
  }
  return(F)
}

### Copule gaussienne ###

C_R_2D <- function(u, v, R){
  return(mvtnorm::pmvnorm(lower = c(-Inf, -Inf), c(qnorm(u,0,1), qnorm(v,0, 1)), corr = R, sigma=NULL, algorithm = mvtnorm::GenzBretz(), keepAttr=FALSE))
}


### densité associée à la copule gaussienne ###

c_R_2D <- function(x1, x2, rho) + return(exp(-0.5*((rho**2)*(qnorm(x1,0,1)**2+qnorm(x2,0,1)**2)-2*rho*qnorm(x1,0,1)*qnorm(x2,0,1))/(1-rho**2))/sqrt(1-(rho**2)))


### Log-vraisemblance dans le cas continu/continu ###

L_n_CC <- function(theta, F1, F2){
  L = 0
  delta = 10**-9
  
  for (i in 1:length(F1)){
    L <- L + log(max(delta,c_R_2D(F1[i], F2[i], theta)))
    
  }
  return(-L)
}

### Copule gaussienne ###


### Log-Vraisemblance dans le cas discret/discret ###

L_n_DD <- function(theta, F1m, F1p, F2m, F2p){
  R = matrix(c(1, theta, theta, 1), 2, 2)
  delta = 10**-9
  L = 0
  for (j in 1:length(F1m)){
    S = C_R_2D(F1p[j], F2p[j], R) + C_R_2D(F1m[j], F2m[j], R) - C_R_2D(F1p[j], F2m[j], R) - C_R_2D(F1m[j], F2p[j], R)
    L <- L + log(max(S, delta))
  }
  return(-L)
}

### Log-Vraisemblance dans le cas continu/discret ###

L_n_CD <- function(theta, F1, F2m, F2p){
  L = 0
  delta = 10**-9
  for (i in (1:length(F1))){
    integr <- integrate(c_R_2D, lower = F2m[i], upper = F2p[i], x1 = F1[i], rho = theta)$value
    L <- L + log(max(delta,integr))
  }
  return(-L)
}

#### type des listes

continu_discret <- function(listeX){
  list_type = c()
  for (i in 1:dim(listeX)[2]){
    X_i = listeX[,i]
    if (all(X_i%%1 == rep(0, length(X_i)))){
      list_type = c(list_type, "D")
    }
    else {list_type = c(list_type, "C")}
  }
  return(list_type)
}

#' rho_estim
#'
#' @description given d random vectors, this function estimates the correlation matrix associated to the gaussian copula C such that the joint distribution function of the random vectors is C(F1(.), ..., Fd(.)) with F1, ..., Fd the distribution functions of the margins.
#'
#' @param list_var list, list of d random vectors with the same length
#' @param Shape The shape returned by the function ("vector" or "matrix")
#'
#' @return A correlation matrix or a vector with the estimates of the correlation coefficients
#' @import mvtnorm
#'
#' @examples list_X = CopulaSim(100, c(3,1,4), c(0.5,0.8,0.3), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(5,4,4), random = TRUE)
#'           rho_estim(list_X, "matrix")
#'
#' @export

rho_estim <- function(list_var, Shape){
  Type = continu_discret(list_var)
  F = list_fdr_d(list_var, Type)
  rho = c()
  M_rho = diag(length(Type))
  
  for (i in 1:(length(Type)-1)){
    for (j in (i+1):length(Type)){
      if (Type[i] == "C" & Type[j] == "C"){
        rho_ij <- optimize(L_n_CC, c(-1,1), F1 = F[,i,1], F2 = F[,j,1], maximum = FALSE)$minimum
      }
      if(Type[i] == "C" & Type[j] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = F[,i,1], F2m = F[,j,1], F2p = F[,j,2], maximum = FALSE)$minimum
      }
      if(Type[j] == "C" & Type[i] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = F[,j,1], F2m = F[,i,1], F2p = F[,i,2], maximum = FALSE)$minimum
      }
      if(Type[j] == "D" & Type[i] == "D"){
        rho_ij <- optimize(L_n_DD, c(-1,1), F1m = F[,i,1], F1p = F[,i,2], F2m = F[,j,1], F2p = F[,j,2], maximum = FALSE)$minimum
      }
      rho = c(rho, rho_ij)
      M_rho[i,j] = rho_ij
      M_rho[j,i] = rho_ij
    }
  }
  if (Shape == "matrix"){return(M_rho)}
  else {return(rho)}
}

#thresholding functions
matrix_cor_01 <- function(M_est, TS){
  M_ = M_est
  M_[which(M_ <= TS)] = 0
  M_[which(M_ != 0)] = 1
  return(M_)
}

matrix_cor_w <- function(M_est, TS){
  M_ = M_est
  M_[which(M_ <= TS)] = 0
  return(M_)
}


#Plotting the graph

#' cor_network_graph
#' @description Create a correlation networks for heterogeneous data
#'
#' @param list_var list, list of d random vectors with the same length
#' @param TS scalar, defined the threshold at which the random variables are considered as being correlated
#' @param vertex_col a color for the vertices of the graph (default color : light blue)
#' @param Weighted boolean, if TRUE, edges will be weighted (the smaller the estimated coefficient the longer the edge). (Default to FALSE)
#'
#' @return a graph where the vertices are correlations and the edges are random variables.
#' @examples
#' list_X = list_XY_d(100, c(3,1,4,5), c(0.5,0.8,0.3), list(qexp, qnorm, qbinom), list(3, c(0,1), c(20,0.5)), c(5,4,4), random = TRUE)
#' cor_network_graph(list_X, 0.3, TRUE, "orange")
#' 
#' @export

cor_network_graph <- function(list_var, TS, Weighted = FALSE, vertex_col = "lightsteelblue2"){
  
  M = rho_estim(list_var, "matrix")
  if (Weighted == FALSE){
    network_ref <- igraph::graph_from_adjacency_matrix(matrix_cor_01(M, TS), mode="undirected", diag=F)
  }
  else{network_ref <- igraph::graph_from_adjacency_matrix(matrix_cor_w(M, TS), mode="undirected", weighted = TRUE, diag=F)}
  
  par(bg="grey", mar=c(1,1,1,1))
  plot(network_ref,
       vertex.size=12,
       vertex.label = names(list_var),
       vertex.color=vertex_col,
       vertex.label.cex=0.8,
       vertex.label.color="black",
       vertex.frame.color="black",
       edge.color = "black")
  
  text(1,1, stringr::str_glue("Réseau de corrélations (seuil = ",TS,")") ,col="black", cex=1.5)
}

