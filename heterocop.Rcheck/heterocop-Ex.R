pkgname <- "heterocop"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "heterocop-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('heterocop')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CopulaSim")
### * CopulaSim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CopulaSim
### Title: CopulaSim
### Aliases: CopulaSim

### ** Examples

M <- diag_block_matrix(c(3,4,5),c(0.7,0.8,0.2))
CopulaSim(20,M,c(rep("qnorm(0,1)",6),rep("qexp(0.5)",4),rep("qbinom(4,0.8)",2)),random=TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CopulaSim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cor_network_graph")
### * cor_network_graph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cor_network_graph
### Title: cor_network_graph
### Aliases: cor_network_graph

### ** Examples

R <- diag_block_matrix(c(3,4,5),c(0.7,0.8,0.2))
data <- CopulaSim(20,R,c(rep("qnorm(0,1)",6),rep("qexp(0.5)",4),
rep("qbinom(4,0.8)",2)),random=FALSE)[[1]]
cor_network_graph(R,TS=0.3,binary=TRUE,legend=c(rep("Normal",6),
rep("Exponential",4),rep("Binomial",2)))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cor_network_graph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("diag_block_matrix")
### * diag_block_matrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: diag_block_matrix
### Title: diag_block_matrix
### Aliases: diag_block_matrix

### ** Examples

diag_block_matrix(c(3,4,5),c(0.3,0.4,0.8))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("diag_block_matrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("gauss_gen")
### * gauss_gen

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: gauss_gen
### Title: gauss_gen
### Aliases: gauss_gen

### ** Examples

M <- diag_block_matrix(c(3,4,5),c(0.7,0.8,0.2))
gauss_gen(M,20)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("gauss_gen", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("matrix_cor_ts")
### * matrix_cor_ts

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: matrix_cor_ts
### Title: matrix_cor_ts
### Aliases: matrix_cor_ts

### ** Examples

M <- diag_block_matrix(c(3,4,5),c(0.7,0.8,0.2))
matrix_cor_ts(M,0.5)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("matrix_cor_ts", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("matrix_gen")
### * matrix_gen

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: matrix_gen
### Title: matrix_gen
### Aliases: matrix_gen

### ** Examples

matrix_gen(15,0.81)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("matrix_gen", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rho_estim")
### * rho_estim

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rho_estim
### Title: rho_estim
### Aliases: rho_estim

### ** Examples

M <- diag_block_matrix(c(3,4,5),c(0.7,0.8,0.2))
data <- CopulaSim(20,M,c(rep("qnorm(0,1)",6),rep("qexp(0.5)",4),
rep("qbinom(4,0.8)",2)),random=FALSE)[[1]]
rho_estim(data,c(rep("C",10),rep("D",2)),parallel=FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rho_estim", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
