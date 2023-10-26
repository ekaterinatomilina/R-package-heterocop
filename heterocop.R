#' GCopula
#'
#' @description This function enables the user to simulate a dataset of n observations by d variables which cumulative distribution function corresponds to a multivariate Gaussian cumulative distribution function of correlation matrix R. In the context of our Gaussian copula model, it corresponds to the simulation of our marginals F_1(X_1), ..., F_d(X_d).
#'
#' @param R the correlation matrix of the copula. It has to be a dxd block matrix.
#' @param n the number of observations
#' 
#' @return The dataset, a data frame containing d vectors (in columns) of length n containing the values F_1(X_1), ..., F_d(X_d) such that the joint distribution function of those vectors is given by F(.,.) = C(F_1(.),..., F_d(.)) where C is the Gaussian Copula of correlation matrix R. F_1,.., F_d denote the marginal distribution functions of the variables X_1, ..., X_d.
#' 
#' @examples M <- diag_block_matrix(c(3,1,4), c(0.2,0.1,0.5))
#' @examples GCopula(M, 100)
#' 
#' @export


GCopula <- function(R, n){
  d = dim(R)[1]
  A = as.matrix(eigen(R)$vectors%*%diag(sqrt(eigen(R)$values)))
  z = matrix(rnorm(d*n, mean = 0, sd = 1), ncol=n)
  data <- data.frame(t(pnorm(A%*%z)))
  return(data)
}

#an intermediate function to generate a block-wise diagonal matrix
diag_block_matrix <- function(blocs, coeff){
  step=c(0,cumsum(blocs))
  d = sum(blocs)
  R = matrix(0,d,d)
  for(i in 1:(length(step)-1)){
    R[((step[i]+1):step[i+1]),((step[i]+1):step[i+1])]=coeff[i]
  }
  diag(R)=1
  return(R)
}

#Now, we are going to use it in order to simulate d vectors linked by a Gaussian copula

#' CopulaSim
#'
#' @description This function enables the user to simulate a dataset of n observations by d variables which cumulative distribution function corresponds to a Gaussian copula. The simulation is done via the generalized inverse method.
#'
#' @param n Integer corresponding to the number of observations.
#' @param blocks vector which contains the size of the blocks on the diagonal of the matrix
#' @param coeff vector of coefficients associated to each block (there is no need to mention the "1" coefficient for one dimensional blocks)
#' @param list_qlaws  List, list of  quasi-inverse functions (must be chosen among (qbeta, qbinom, qchisq, qubif, qexp, qf, qgamma, qgeom, qhyper, qnbinom, qnorm, qpois, qt, qweibull)) with their specified parameters
#' @param repetition vector, number of variables following each law
#' @param random if TRUE the laws will be randomly distributed among the variables, if FALSE their order will be preserved. (by default : random = TRUE)
#'
#' @return matrix of d vectors (in columns) of length n such that the joint distribution function of those vectors is given  by F(.,.) = C(F_1(.),..., F_d(.)) where C is the gaussian Copula related to a block diagonal matrix and F_1,.., F_d are respectively the distribution functions of the margins.
#' 
#' @examples CopulaSim(100, c(3,2,4), c(0.2,0.7,0.5), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(3,4,2), random = TRUE)
#' 
#' @export

CopulaSim <- function(n, blocks, coeff, list_qlaws, repetition, random = TRUE){
  
  R = diag_block_matrix(blocks, coeff)
  d = sum(blocks)
  XY = GCopula(R,n)
  
  nb = c()
  for (l in 1:length(repetition)){
    nb = append(nb, rep(l,repetition[l]))
    }
  if (random == TRUE){
    nb = sample(nb)
  }
  
  vars = matrix(0,n,d)
  for (j in 1:d){
    i = nb[j]
    vars[,j] = list_qlaws[[i]](XY[, j])
  }
  return(data.frame(vars))
}


#This function estimates the CDF values for each variable
fdr_d <- function(X, Type){ ###### WARNING: X MUST BE A DATA FRAME!
  F <- array(0,dim=c(dim(X)[1],dim(X)[2],2))
  F[,,2] <- sapply(X,rank,ties.method="max")/(dim(X)[1]+1) 
  F_s <- rbind(0,apply(F[,,2],2,sort))
  for (i in 1:dim(X)[2]){
    if (Type[i] == "D"){
      for(j in 1:dim(X)[1]){
        F[j,i,1]=F_s[min(which(F_s[,i]==F[j,i,2]))-1,i]
      }
    }
  }
  return(F) # OK; SEE ALSO F[,1,2][rank(F[,1,2], ties.method="max")-1] (WARNING: TAKE CARE OF RANK=0)
}

#Gaussian copula expression

c_R_2D <- function(x1, x2, rho) + return(exp(-0.5*((rho**2)*(qnorm(x1,0,1)**2+qnorm(x2,0,1)**2)-2*rho*qnorm(x1,0,1)*qnorm(x2,0,1))/(1-rho**2))/sqrt(1-(rho**2)))

C_R_2D<- function(u1, u2, l1, l2, R){
  return(mvtnorm::pmvnorm(upper=c(qnorm(u1,0,1), qnorm(u2,0,1)),mean=c(0,0),lower=c(qnorm(l1,0,1), qnorm(l2,0,1)),corr = R, sigma=NULL, algorithm = mvtnorm::GenzBretz(), keepAttr=FALSE))
}

#Calculating log-likelihood

L_n_CC <- function(theta, F1, F2){
  delta = 10**-9
  mysummands <- c_R_2D(F1, F2, theta)
  L <- sum(log(mapply(max,mysummands,MoreArgs=list(delta))))
  return(-L)
}

L_n_DD <- function(theta, F1m, F1p, F2m, F2p){
  R = matrix(c(1, theta, theta, 1), 2, 2)
  delta = 10**-9
  mysummands <- mapply(C_R_2D,F1p, F2p, F1m, F2m, MoreArgs=list(R))
  L <- sum(log(mapply(max,mysummands,MoreArgs=list(delta))))
  return(-L)
}

### Log-Vraisemblance dans le cas continu/discret ###
## F1: vector of the CDF values for the continuous variable
## F2p: vector of the CDF values for the discrete variables
## F2m: F2m[i] is the value in F2p that precedes F2p[i]
L_n_CD <- function(theta, F1, F2m, F2p){
  delta=10**-9
  mysummands <- pnorm(qnorm(F2p),mean=theta*qnorm(F1),sd=sqrt(1-theta**2)) - pnorm(qnorm(F2m),mean=theta*qnorm(F1),sd=sqrt(1-theta**2))
  L <-sum(log(mapply(max,mysummands,MoreArgs=list(delta))))
  return(-L)
}


#' rho_estim
#'
#' @description given d random vectors X_1, .., X_d, this function estimates the correlation matrix associated to the gaussian copula C such that the joint distribution function of the random vectors is C(F1(.), ..., Fd(.)) with F1, ..., Fd the distribution functions of the margins.
#'
#' @param data list, list of d random vectors with the same length
#' @param Type a vector of characters corresponding to the type of each function, "C" for continuous and "D" for discrete.
#'
#' @return A correlation matrix with the estimates of the correlation coefficients of the copula.
#' @import mvtnorm
#'
#' @examples data = CopulaSim(100, c(3,2,4), c(0.5,0.8,0.3), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(3,4,2), random = TRUE)
#'           rho_estim(data, rep("C",9))
#'
#' @export

rho_estim <- function(data,Type){
  F = fdr_d(data, Type)
  M_rho = diag(length(Type))
  for (i in 1:(length(Type)-1)){
    for (j in (i+1):length(Type)){
      print(c(i,j))
      if (Type[i] == "C" & Type[j] == "C"){
        rho_ij <- optimize(L_n_CC, c(-1,1), F1 = F[,i,2], F2 = F[,j,2], maximum = FALSE)$minimum
      }
      if(Type[i] == "C" & Type[j] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = F[,i,2], F2m = F[,j,1], F2p = F[,j,2],maximum = FALSE)$minimum
      }
      if(Type[j] == "C" & Type[i] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = F[,j,2], F2m = F[,i,1], F2p = F[,i,2],maximum = FALSE)$minimum
      }
      if(Type[j] == "D" & Type[i] == "D"){
        rho_ij <- optimize(L_n_DD, c(-1,1), F1m = F[,i,1], F1p = F[,i,2], F2m = F[,j,1], F2p = F[,j,2],maximum = FALSE)$minimum
      }
      M_rho[i,j] = rho_ij
      M_rho[j,i] = rho_ij
    }
  }
  return(M_rho)
}

rho_estim_par <- function(data,Type){
  F = fdr_d(data, Type)
  M_rho = diag(length(Type))
  
  Ncpus <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(Ncpus)
  doParallel::registerDoParallel(cl)
  
M_rho <- foreach(i=1:(length(Type)-1), .combine='rbind',.export=c("c_R_2D", "C_R_2D","L_n_CC", "L_n_CD","L_n_DD"))%dopar%{
    rho_i <- c(rep(0,i-1),1)
    for (j in (i+1):length(Type)){
      if (Type[i] == "C" & Type[j] == "C"){
        rho_ij <- optimize(L_n_CC, c(-1,1), F1 = F[,i,2], F2 = F[,j,2], maximum = FALSE)$minimum
      }
      if(Type[i] == "C" & Type[j] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = F[,i,2], F2m = F[,j,1], F2p = F[,j,2],maximum = FALSE)$minimum
      }
      if(Type[j] == "C" & Type[i] == "D"){
        rho_ij <- optimize(L_n_CD, c(-1,1), F1 = F[,j,2], F2m = F[,i,1], F2p = F[,i,2],maximum = FALSE)$minimum
      }
      if(Type[j] == "D" & Type[i] == "D"){
        rho_ij <- optimize(L_n_DD, c(-1,1), F1m = F[,i,1], F1p = F[,i,2], F2m = F[,j,1], F2p = F[,j,2],maximum = FALSE)$minimum
      }
      rho_i <- c(rho_i, rho_ij)
    }
      return(rho_i)
}
    M_rho <- rbind(M_rho, c(rep(0,length(Type)-1),1))
    parallel::stopCluster(cl)
    return(M_rho + t(M_rho)-diag(1,length(Type),length(Type)))
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
#' @description Create a correlation network for heterogeneous data
#'
#' @param list_var list, list of d random vectors with the same length
#' @param TS scalar, defined the threshold at which the random variables are considered as being correlated
#' @param vertex_col a color for the vertices of the graph (default color : light blue)
#' @param Weighted boolean, if TRUE, edges will be weighted (the smaller the estimated coefficient the longer the edge). (Default to FALSE)
#'
#' @return a graph where the vertices are correlations and the edges are random variables.
#' @examples
#' data = CopulaSim(100, c(3,2,4), c(0.5,0.8,0.3), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(5,2,2), random = TRUE)
#' cor_network_graph(data, rep("C",9),0.3, TRUE, "orange")
#' 
#' @export

cor_network_graph <- function(data, Type, TS, Weighted = FALSE){
  
  M = rho_estim(data, Type)
  if (Weighted == FALSE){
    network_ref <- igraph::graph_from_adjacency_matrix(matrix_cor_01(M, TS), mode="undirected", diag=F)
  }
  else{network_ref <- igraph::graph_from_adjacency_matrix(matrix_cor_w(M, TS), mode="undirected", weighted = TRUE, diag=F)}
  
  par(bg="white", mar=c(1,1,1,1))
  plot(network_ref,
       vertex.size=18,
       vertex.label = c(1:30),
       vertex.color=c(rep("orange1",3), rep("cornflowerblue",5), rep("orange1",4),rep("cornflowerblue",8),rep("mediumseagreen",8),rep("orange1",2)),
       vertex.label.cex=1.2,
       vertex.label.color="black",
       vertex.frame.color="black",
       edge.color = "black")
  legend(1.2,0.8,c("N(0.5,2)", "P(4)","B(0.7)"),fill=c("orange1","cornflowerblue","mediumseagreen"))
  text(1.1,1, stringr::str_glue("threshold = ",TS) ,col="black", cex=1.5)
}

