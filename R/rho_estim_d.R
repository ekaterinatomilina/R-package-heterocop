

list_fdr_d <- function(listeX, Type){
  list_F = c()
  for (i in 1:length(listeX)){
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



#' rho_estim_d
#'
#' @param list_var list, list of d random vectors with the same length
#' @param Shape The shape returned by the function ("vector" or "matrix")
#'
#' @return A correlation matrix or a vector with the estimates of the correlation coefficients
#' @import mvtnorm
#'
#' @examples rho_estim_d(list(X,Y,Z), "matrix")
rho_estim_d <- function(list_var, Shape){
  Type = continu_discret(list_var)
  listeF = list_fdr_d(list_var, Type)
  rho = c()
  M_rho = diag(length(Type))

  for (i in 1:(length(Type)-1)){
    for (j in (i+1):length(Type)){
      if (Type[[i]] == "C" & Type[[j]] == "C"){
        rho_ij <- optimize(L_n_CC, c(-1,1), F1 = unlist(listeF[[i]]), F2 = unlist(listeF[[j]]), maximum = FALSE)$minimum
      }
      if(Type[[i]] == "C" & Type[[j]] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = unlist(listeF[[i]]), F2m = unlist(listeF[[j]][[1]]), F2p = unlist(listeF[[j]][[2]]), maximum = FALSE)$minimum
      }
      if(Type[[j]] == "C" & Type[[i]] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = unlist(listeF[[j]]), F2m = unlist(listeF[[i]][[1]]), F2p = unlist(listeF[[i]][[2]]), maximum = FALSE)$minimum
      }
      if(Type[[j]] == "D" & Type[[i]] == "D"){
        rho_ij <- optimize(L_n_DD, c(-1,1), F1m = unlist(listeF[[i]][[1]]), F1p = unlist(listeF[[i]][[2]]), F2m = unlist(listeF[[j]][[1]]), F2p = unlist(listeF[[j]][[2]]), maximum = FALSE)$minimum
      }
      rho = c(rho, rho_ij)
      M_rho[i,j] = rho_ij
      M_rho[j,i] = rho_ij
    }
  }
  if (Shape == "matrix"){return(M_rho)}
  else {return(rho)}
  return(M_rho)
}
