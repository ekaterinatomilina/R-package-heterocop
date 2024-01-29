#' matrix_gen
#' 
#' @description This function enables the user to generate a positive semi-definite sparse matrix of size d x d and with initial sparsity parameter gamma. 
#' 
#' @param d the number of columns/rows of the matrix
#' @param gamma initial sparsity parameter
#' 
#' @return a list containing the generated matrix and its new sparsity parameter (ie the proportion of zeros)
#' @examples matrix_gen(30,0.2)
#' @export

matrix_gen <- function(d,gamma){
  if(gamma >= 0 & gamma <=1){
L <- matrix(0,d,d)
params <- runif(d*(d-1)/2,0.3,1)
ind <- sample(1:(d*(d-1)/2),floor(gamma*d*(d-1)/2))
params[ind] <- 0
L[lower.tri(L)] <- params
diag(L) <- 1
R <- L%*%t(L)
C <- round(cov2cor(R),3)
gamma_f = sum(C==0)/(d*d)
return(list(as.matrix(C),gamma_f))}
  else{
    stop("gamma must be between 0 and 1")
  }
}

#' diag_block_matrix
#' 
#' @description This function enables the user to generate a diagonal matrix of size d x d 
#' 
#' @param blocks a vector containing the size of each block
#' @param coeff a vector containing the coefficient corresponding to each block
#' 
#' @return the generated block matrix
#' @examples M <- diag_block_matrix(c(3,1,4), c(0.2,0.1,0.5))
#' @export

diag_block_matrix <- function(blocks, coeff){
  if(!is.numeric(coeff)){
    stop("The coefficients must be numeric.")
  }else{
  step=c(0,cumsum(blocs))
  d = sum(blocs)
  R = matrix(0,d,d)
  for(i in 1:(length(step)-1)){
    R[((step[i]+1):step[i+1]),((step[i]+1):step[i+1])]=coeff[i]
  }
  diag(R)=1
  return(R)
  }
}

#' GCopula
#'
#' @description This function enables the user to simulate a dataset of n observations by d variables which cumulative distribution function corresponds to a multivariate Gaussian cumulative distribution function of correlation matrix R. In the context of our Gaussian copula model, it corresponds to the simulation of our marginals F_1(X_1), ..., F_d(X_d).
#'
#' @param R the correlation matrix of the copula. It has to be a dxd positive semi-definite matrix.
#' @param n the number of observations
#' 
#' @return The dataset, a data frame containing d vectors (in columns) of length n containing the values F_1(X_1), ..., F_d(X_d) such that the joint distribution function of those vectors is given by F(.,.) = C(F_1(.),..., F_d(.)) where C is the Gaussian Copula of correlation matrix R. F_1,.., F_d denote the marginal distribution functions of the variables X_1, ..., X_d.
#' 
#' @examples M <- diag_block_matrix(c(3,1,4), c(0.2,0.1,0.5))
#'  GCopula(M, 100)
#' @import matrixcalc
#' @export


GCopula <- function(R, n){
  if(sum(R>1 | R< -1)==0){ #checking that all coefficients are between -1 and 1
  if(matrixcalc::is.positive.semi.definite(R)){ # checking that the matrix is semi-positive definite
  d = dim(R)[1]
  A = as.matrix(eigen(R)$vectors%*%diag(sqrt(eigen(R)$values)))
  z = matrix(rnorm(d*n, mean = 0, sd = 1), ncol=n)
  data <- data.frame(t(pnorm(A%*%z)))
  return(data)
  }
  else{
    stop("The matrix is not positive semi-definite")
  }
  }
  else{
    stop("The coefficients in the matrix must be between -1 and 1")
  }
}

#' CopulaSim
#'
#' @description This function enables the user to simulate a dataset of n observations by d variables which cumulative distribution function corresponds to a Gaussian copula. The simulation is done via the generalized inverse method.
#'
#' @param n Integer corresponding to the number of observations.
#' @param blocks vector which contains the size of the blocks on the diagonal of the matrix
#' @param coeff vector of coefficients associated to each block (there is no need to mention the "1" coefficient for one dimensional blocks)
#' @param list_qlaws  List, list of  quasi-inverse functions (must be chosen among (qbeta, qbinom, qchisq, qubif, qexp, qf, qgamma, qgeom, qhyper, qnbinom, qnorm, qpois, qt, qweibull)) with their specified parameters
#' @param repetition vector, number of variables following each law
#' @param random if TRUE the probability distributions will be randomly distributed among the variables, if FALSE their order will be preserved. (by default : random = TRUE)
#' @param vartype a vector containing the letters "C" (continuous) and "D"(discrete) denoting the type of each variable, in order
#'
#' @return a list containing a matrix of d vectors (in columns) of length n such that the joint distribution function of those vectors is given  by F(.,.) = C(F_1(.),..., F_d(.)) where C is the gaussian Copula related to a block diagonal matrix and F_1,.., F_d are respectively the distribution functions of the margins, and a vector containing the re-ordered variable types
#' 
#' @examples M <- diag_block_matrix(c(3,2,4), c(0.2,0.1,0.5))
#' CopulaSim(100, M, list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(3,4,2), random = TRUE)
#' 
#' @export

CopulaSim <- function(n, R, list_qlaws, repetition, random = FALSE,vartype){
  if(dim(R)[1]==sum(repetition)){
  d = sum(repetition)
  XY = GCopula(R,n)
  names <- colnames(XY)
  
  nb = c()
  for (l in 1:length(repetition)){
    nb = append(nb, rep(l,repetition[l]))
  }
  if (random == TRUE){
    idx = sample(1:length(nb))
    nb <- nb[idx]
    names <- names[idx]
    vartype <- vartype[idx]
  }
  
  vars = matrix(0,n,d)
  for (j in 1:d){
    i = nb[j]
    vars[,j] = list_qlaws[[i]](XY[, j])
  }
  colnames(vars) <- names
  
  return(list(data.frame(vars),vartype))
  }else{
    stop("The total number of distribution repetitions must be equal to the size of the matrix")
  }
}

#' fdr_d
#'
#' @description This function enables the user to estimate the empirical cumulative distribution function of each variable
#' 
#' @param X data frame of d variables and n observations
#' @param Type a vector of length d containing the type of each variable, "C" denoting continuous, "D" denoting discrete
#'
#' @return an array of size 2 x n x d where the first dimension contains the n x d values for F- (0 for continuous variables), the second one contains the n x d values for F
#' 
#' @examples M <- diag_block_matrix(c(3,2,4), c(0.2,0.1,0.5))
#'  data <- CopulaSim(100, M, list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(3,4,2), random = TRUE)
#'  fdr_d(data[[1]],data[[2]])
#' 
#' @export

fdr_d <- function(X, Type){
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
  return(F)
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

### Log-likelihood in the mixed case ###

L_n_CD <- function(theta, F1, F2m, F2p){
  delta=10**-9
  mysummands <- pnorm(qnorm(F2p),mean=theta*qnorm(F1),sd=sqrt(1-theta**2)) - pnorm(qnorm(F2m),mean=theta*qnorm(F1),sd=sqrt(1-theta**2))
  L <-sum(log(mapply(max,mysummands,MoreArgs=list(delta))))
  return(-L)
}


#' rho_estim
#'
#' @description given d random vectors X_1, .., X_d, this function estimates the correlation matrix associated to the gaussian copula C such that the joint distribution function of the random vectors is C(F1(.), ..., Fd(.)) with F1, ..., Fd the cumulative distribution functions of the margins.
#'
#' @param data data frame containing the random vectors in columns
#' @param Type a vector of characters corresponding to the type of each function, "C" for continuous and "D" for discrete.
#' @param parallel TRUE if you want the code to be parallelized, FALSE if not
#'
#' @return A correlation matrix with the estimates of the correlation coefficients of the copula.
#' @import mvtnorm
#' @import doParallel
#'
#' @examples data = CopulaSim(100, c(3,2,4), c(0.5,0.8,0.3), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(3,4,2), random = TRUE,rep("C",9))
#'           rho_estim(data[[1]], data[[2]])
#'
#' @export

rho_estim <- function(data,Type,parallel=TRUE){
  
  F = fdr_d(data, Type)
  M_rho = diag(length(Type))
  
  if(parallel==FALSE){
  
  for (i in 1:(length(Type)-1)){
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
      M_rho[i,j] = rho_ij
      M_rho[j,i] = rho_ij
    }
  }
  return(M_rho)
}else{
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



#' cor_network_graph
#' @description Create a correlation network for heterogeneous data
#'
#' @param data the dataset containing d variables (columns) on n observations
#' @param Type a vector containing the type of each variable or the group it belongs to
#' @param TS scalar, defined the threshold at which the random variables are considered as being correlated
#' @param Weighted boolean, if TRUE, edges will be weighted (the smaller the estimated coefficient the longer the edge). (Default to FALSE)
#' @param color a vector containing colors for the vertices
#'@param legend a vector containing the legend
#'
#' @return a graph where the vertices are correlations and the edges are random variables.
#' @examples
#' data = CopulaSim(100, c(3,2,4), c(0.5,0.8,0.3), list(function(p){qexp(p=p, rate=1)}, function(p){qnorm(p=p, mean=0, sd=1)}, function(p){qnorm(p=p, mean=0.76, sd=2)}), c(5,2,2), random = TRUE,vartype=rep("C",9))
#' cor_network_graph(data[[1]],data[[2]],0.3, TRUE)
#' 
#' @export

cor_network_graph <- function(data, Type, TS, parallel=TRUE,Weighted = FALSE,color,legend){
  
  M = rho_estim(data, Type,parallel)
  if (Weighted == FALSE){
    network_ref <- igraph::graph_from_adjacency_matrix(matrix_cor_01(M, TS), mode="undirected", diag=F)
  }
  else{network_ref <- igraph::graph_from_adjacency_matrix(matrix_cor_w(M, TS), mode="undirected", weighted = TRUE, diag=F)}
  
  par(bg="white", mar=c(1,1,1,1))
  plot(network_ref,
       vertex.size=18,
       vertex.label = names(data),
       vertex.color= as.factor(color),
       vertex.label.cex=0.8,
       vertex.label.color="black",
       vertex.frame.color="black",
       edge.color = "black")
  legend(1.2,0.8,legend,fill=c(as.factor(color)))
  text(1.1,1, stringr::str_glue("threshold = ",TS) ,col="black", cex=1.5)
}

