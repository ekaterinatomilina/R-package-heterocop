
### fonction de répartition empirique
nombre_groupe <- function(y_var, list_Y)+return(length(which(list_Y==y_var)))

list_fdr <- function(Z, arg ){
  if(arg == "m"){F <- as.numeric((rank(Z, ties.method = "max")-sapply(Z,nombre_groupe,list_Y = Z))/(length(Z)+1))}
  else {F <- as.numeric(rank(Z, ties.method = "max")/(length(Z)+1))}
}

list_fdr_2 <- function(listeX, Type){
  list_F = c()
  for (i in 1:2){
    X_i = listeX[[i]]
    if (Type[[i]] == "C"){list_F[[i]] = list(list_fdr(X_i,1))}
    else{
      Fp = list_fdr(X_i, "p")
      Fm = list_fdr(X_i, "m")
      list_F[[i]] = list(Fm,Fp)
    }
  }
  return(list_F)
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
  for (i in 1:length(listeX)){
    X_i = listeX[[i]]
    if (all(X_i%%1 == rep(0, length(X_i)))){
      list_type = c(list_type, "D")
    }
    else {list_type = c(list_type, "C")}
  }
  return(list_type)
}

#### Fonction qui permet d'estimer le rho associé

#' rho_estim_2
#'
#' @description  given two random vectors, this function estimates the correlation coefficient associated to the gaussian copula C such that the joint distribution function of the random vectors is C(F1(.), F2(.)) with  F1 and F2 the margins distribution function.
#'
#' @param X list, list of iid random variables
#' @param Y list, list of iid random variables (with the same length as X)
#'
#' @return an estimation of the correlation matrix parameter
#'
#' @examples XY = list_XY(matrix(c(1,0.5,0.5,1),2,2), 100, list(loi1 = qexp, loi2 = qnorm), list(3,c(0,1)))
#'           rho_estim_2(XY[[1]], XY[[2]])

#' @export
rho_estim_2 <- function(X, Y){

  list_var = list(X,Y)
  Type = continu_discret(list_var)
  listeF = list_fdr_2(list_var, Type)

  rho = c()
  if (Type[[1]] == "C" & Type[[2]] == "C"){
    rho_12 <- optimize(L_n_CC, c(-1,1), F1 = unlist(listeF[[1]]), F2 = unlist(listeF[[2]]), maximum = FALSE)$minimum
  }
  if(Type[[1]] == "C" & Type[[2]] == "D"){
    rho_12 <- optimize(L_n_CD, c(-1,1), F1 = unlist(listeF[[1]]), F2m = unlist(listeF[[2]][[1]]), F2p = unlist(listeF[[2]][[2]]), maximum = FALSE)$minimum
  }
  if(Type[[2]] == "C" & Type[[1]] == "D"){
    rho_12 <- optimize(L_n_CD, c(-1,1), F1 = unlist(listeF[[2]]), F2m = unlist(listeF[[1]][[1]]), F2p = unlist(listeF[[1]][[2]]), maximum = FALSE)$minimum
  }
  if(Type[[2]] == "D" & Type[[1]] == "D"){
    rho_12 <- optimize(L_n_DD, c(-1,1), F1m = unlist(listeF[[1]][[1]]), F1p = unlist(listeF[[1]][[2]]), F2m = unlist(listeF[[2]][[1]]), F2p = unlist(listeF[[2]][[2]]), maximum = FALSE)$minimum
  }
  rho = c(rho, rho_12)
  return(rho)
}



