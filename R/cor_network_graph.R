
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
#'
#' @param list_var list, list of d random vectors with the same length
#' @param TS scalar, defined the threshold at which the random variables are considered as being correlated
#' @param vertex_col a color for the vertices of the graph (default color : ligth blue)
#' @param Weighted boolean, if TRUE, edges will be weighted (the smaller the estimated coefficient the longer the edge). (Default to FALSE)
#'
#' @return a graph where the vertices are correlations and the edges are random variables.
#' @export
#'
#'
#' @examples cor_network_graph_ (list(X,Y,Z), 0.3, TRUE, "orange")
#'
cor_network_graph <- function(list_var, TS, Weighted = FALSE, vertex_col = "lightsteelblue2"){

  M = rho_estim_d(list_var, "matrix")
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
