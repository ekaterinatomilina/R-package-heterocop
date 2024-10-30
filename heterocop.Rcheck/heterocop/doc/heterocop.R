## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(heterocop)
require(dplyr)
require(kableExtra)
require(knitr)

## ----warning=FALSE,message=FALSE----------------------------------------------
data(icgc_data)
R <- rho_estim(icgc_data,c(rep("C",10),rep("D",5)),parallel=FALSE)

## ----echo=FALSE,warning=FALSE,message=FALSE-----------------------------------
knitr::kable(head(R[,1:6]), digits = 2) 

## -----------------------------------------------------------------------------
cor_network_graph(R,TS=0.3,legend=c(rep("RNAseq",5),rep("Proteins",5),rep("Mutations",5)))

## -----------------------------------------------------------------------------
R <- diag_block_matrix(c(3,2),c(0.4,0.8))

## ----echo=FALSE,warning=FALSE,message=FALSE-----------------------------------
knitr::kable(R, digits = 2, col.names = NULL)

## -----------------------------------------------------------------------------
R <- matrix_gen(5,0.81)

## ----echo=FALSE,warning=FALSE,message=FALSE-----------------------------------
knitr::kable(R, digits = 2, col.names = NULL)

## -----------------------------------------------------------------------------
R <- diag_block_matrix(c(3,5,2),c(0.7,0.3,0.5))
CopulaSim(5,R,c(rep("qnorm(0,1)",5),rep("qexp(0.5)",3),rep("qbinom(4,0.8)",2)),random=TRUE)

## -----------------------------------------------------------------------------
latent_data <- gauss_gen(R,10)

## ----echo=FALSE,warning=FALSE,message=FALSE-----------------------------------
knitr::kable(latent_data, digits = 2)%>%
kableExtra::kable_styling(font_size = 8, full_width = FALSE)

