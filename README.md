heterocop: an R package for Gaussian copula semi-parametric inference
for mixed data
================
Ekaterina Tomilina
2024-10-18

This package enables the user to quantify dependencies between mixed
(continuous, discrete, binary) variables in the framework of a Gaussian
copula model by estimating the correlation matrix of the copula
(Tomilina, Mazo, and Jaffrézic (2024)).

## Context

When working with $d$ mixed variables $X_1, ..., X_d$, it can be
complicated to infer a correlation network as well-known statistical
methods do not work in case of mixed data. Indeed, most classical
network inference methods such as gLasso or Gaussian graphical models
rely on the assumption that the data follow a Gaussian distribution.
Recently, the Non-paranormal distribution was introduced for network
inference for continuous, non-Gaussian data (Liu (2009)). It consists in
a transformation of the cumulative distribution functions via a Gaussian
copula and provides results for non-Gaussian but continuous variables.
We propose an extension of this model to the case of mixed variables.

## The Model

Let $X_1, ..., X_d$ be mixed variables (continuous or discrete). Let
$F_1, ..., F_d$ denote their marginal CDFs, $\Phi^{-1}$ the inverse of
the standard normal CDF and $\Phi_\Sigma$ the Gaussian CDF of
correlation matrix $\Sigma$. We assume that the multivariate CDF of the
vector $(X_1,\dots,X_d)$ is given by:
$$F(X_1, ..., X_d)=C_\Sigma(F_1(X_1), ..., F_d(X_d)):=\Phi_\Sigma(\Phi^{-1}(F_1(X_1)), ..., \Phi^{-1}(F_d(X_d))).$$

# Estimation

In order to estimate the correlation matrix of the copula, the rho_estim
function uses a semiparametric pairwise maximum likelihood estimator. It
returns the estimated correlation matrix of the copula and takes as
arguments the data set and the variable types in a vector. In the
example below, we have used a subset of the ICGC data set (Zhang,
Bajari, and D. (2019)) which contains 5 RNA-seq, 5 protein and 5
mutation variables. We have specified the variable types, where a “C”
stands for “continuous” and a “D” for “discrete”.

``` r
data(icgc_data)
R <- rho_estim(icgc_data,c(rep("C",10),rep("D",5)))
```

A $6\times6$ subset of the obtained copula correlation matrix is
represented below.

|         | ACACA | AKT1S1 |  ANLN | ANXA1 |    AR | ACACA_P |
|:--------|------:|-------:|------:|------:|------:|--------:|
| ACACA   |  1.00 |  -0.13 |  0.09 | -0.32 |  0.36 |    0.58 |
| AKT1S1  | -0.13 |   1.00 | -0.13 | -0.16 | -0.19 |    0.11 |
| ANLN    |  0.09 |  -0.13 |  1.00 |  0.13 | -0.36 |   -0.07 |
| ANXA1   | -0.32 |  -0.16 |  0.13 |  1.00 | -0.36 |   -0.22 |
| AR      |  0.36 |  -0.19 | -0.36 | -0.36 |  1.00 |    0.14 |
| ACACA_P |  0.58 |   0.11 | -0.07 | -0.22 |  0.14 |    1.00 |

# Graphical representation

The cor_network_graph function enables to visualize the obtained
network. It takes as arguments the dataset, the correlation matrix, the
threshold and a legend.

``` r
cor_network_graph(R,TS=0.3,legend=c(rep("RNAseq",5),rep("Proteins",5),rep("Mutations",5)))
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

## Simulation

Our package is also able to simulate data distributed according to the
Gaussian copula model. Two functions enable us to generate two types of
correlation matrices: block-wise and sparse. The diag_block_matrix
function enables the user to get a block-wise correlation matrix. It
takes as arguments a vector containing the block sizes and a vector
containing the coefficients of each block. An example is shown below.

``` r
R <- diag_block_matrix(c(3,2),c(0.4,0.8))
```

|     |     |     |     |     |
|----:|----:|----:|----:|----:|
| 1.0 | 0.4 | 0.4 | 0.0 | 0.0 |
| 0.4 | 1.0 | 0.4 | 0.0 | 0.0 |
| 0.4 | 0.4 | 1.0 | 0.0 | 0.0 |
| 0.0 | 0.0 | 0.0 | 1.0 | 0.8 |
| 0.0 | 0.0 | 0.0 | 0.8 | 1.0 |

The matrix_gen function enables the user to generate a sparse
correlation matrix of initial sparsity parameter $\gamma$, which has to
be specified. It is based on the Cholesky decomposition where a lower
triangular matrix of sparsity $\gamma$ is generated before being
multiplied by its transpose in order to obtain the final matrix. Note
that the initial parameter is not equal to the final parameter, which is
also returned by the function. In the example below, the first element
of the list is the resulting matrix, and the second element of the list
is the final sparsity parameter.

``` r
R <- matrix_gen(5,0.81)
```

<table class="kable_wrapper">
<tbody>
<tr>
<td>

|      |      |     |      |      |
|-----:|-----:|----:|-----:|-----:|
| 1.00 | 0.00 |   0 | 0.66 | 0.00 |
| 0.00 | 1.00 |   0 | 0.00 | 0.49 |
| 0.00 | 0.00 |   1 | 0.00 | 0.00 |
| 0.66 | 0.00 |   0 | 1.00 | 0.00 |
| 0.00 | 0.49 |   0 | 0.00 | 1.00 |

</td>
<td>

|                 |
|:----------------|
| sparsity = 0.64 |

</td>
</tr>
</tbody>
</table>

The CopulaSim function enables the user to generate a data set which CDF
can be expressed as a Gaussian copula of correlation matrix R (to be
specified). In the example below, we first generate a block diagonal
correlation matrix R and then generate the data set. Then, CopulaSim
takes as arguments the number of observations, the correlation matrix of
the copula, a vector containing the probability distributions and their
parameters, the number of repetitions of each distribution, and enables
the user to randomize their order. It returns a list of three elements:
the data frame containing the generated data, the correlation matrix,
and the permutation realized on the rows and columns of R order after
randomization.

``` r
R <- diag_block_matrix(c(3,5,2),c(0.7,0.3,0.5))
CopulaSim(5,R,c(rep("qnorm(0,1)",5),rep("qexp(0.5)",3),rep("qbinom(4,0.8)",2)),random=TRUE)
#> [[1]]
#>             X1         X2         X3          X4         X5        X6        X7
#> 1 -1.089134252  0.5232069  0.3876855  1.36404502 2.22151961 1.3530549 3.3331100
#> 2 -0.002139848  0.2751454 -0.7825396 -0.19900487 0.59975407 4.1010099 0.3227235
#> 3  2.359786314 -1.1978236  2.6050133  0.44768577 0.39552340 7.2441257 3.2497848
#> 4  1.187899455  0.4613101 -0.8278819  0.59213314 2.06272144 0.8695877 0.7899809
#> 5  0.114799322  1.7587785 -1.7438769  0.08446285 0.00618523 0.9913152 0.6635844
#>         X8 X9 X10
#> 1 7.895351  4   4
#> 2 2.387090  4   4
#> 3 1.356079  4   4
#> 4 4.287998  3   4
#> 5 1.716021  2   3
#> 
#> [[2]]
#>       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#>  [1,]  1.0  0.0  0.0  0.3  0.0  0.3  0.0  0.3  0.0   0.3
#>  [2,]  0.0  1.0  0.0  0.0  0.0  0.0  0.7  0.0  0.7   0.0
#>  [3,]  0.0  0.0  1.0  0.0  0.5  0.0  0.0  0.0  0.0   0.0
#>  [4,]  0.3  0.0  0.0  1.0  0.0  0.3  0.0  0.3  0.0   0.3
#>  [5,]  0.0  0.0  0.5  0.0  1.0  0.0  0.0  0.0  0.0   0.0
#>  [6,]  0.3  0.0  0.0  0.3  0.0  1.0  0.0  0.3  0.0   0.3
#>  [7,]  0.0  0.7  0.0  0.0  0.0  0.0  1.0  0.0  0.7   0.0
#>  [8,]  0.3  0.0  0.0  0.3  0.0  0.3  0.0  1.0  0.0   0.3
#>  [9,]  0.0  0.7  0.0  0.0  0.0  0.0  0.7  0.0  1.0   0.0
#> [10,]  0.3  0.0  0.0  0.3  0.0  0.3  0.0  0.3  0.0   1.0
#> 
#> [[3]]
#>  [1]  7  5  9  4  3  6  8  1 10  2
```

Additionally, the gauss_gen function, which is used in CopulaSim,
generates the latent Gaussian variables linked by the correlation matrix
R. Its only arguments are the correlation matrix R and the number of
observations.

``` r
latent_data <- gauss_gen(R,10)
```

<table class="table" style="font-size: 8px; width: auto !important; margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:right;">
X1
</th>
<th style="text-align:right;">
X2
</th>
<th style="text-align:right;">
X3
</th>
<th style="text-align:right;">
X4
</th>
<th style="text-align:right;">
X5
</th>
<th style="text-align:right;">
X6
</th>
<th style="text-align:right;">
X7
</th>
<th style="text-align:right;">
X8
</th>
<th style="text-align:right;">
X9
</th>
<th style="text-align:right;">
X10
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.90
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.26
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.63
</td>
</tr>
<tr>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.66
</td>
<td style="text-align:right;">
0.65
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.44
</td>
<td style="text-align:right;">
0.55
</td>
<td style="text-align:right;">
0.27
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.37
</td>
</tr>
<tr>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.33
</td>
<td style="text-align:right;">
0.40
</td>
<td style="text-align:right;">
0.68
</td>
<td style="text-align:right;">
0.19
</td>
<td style="text-align:right;">
0.25
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.59
</td>
<td style="text-align:right;">
0.01
</td>
<td style="text-align:right;">
0.00
</td>
</tr>
<tr>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.45
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.30
</td>
<td style="text-align:right;">
0.82
</td>
<td style="text-align:right;">
0.43
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.49
</td>
</tr>
<tr>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.51
</td>
<td style="text-align:right;">
0.85
</td>
<td style="text-align:right;">
0.69
</td>
<td style="text-align:right;">
0.46
</td>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.70
</td>
</tr>
<tr>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.97
</td>
<td style="text-align:right;">
0.88
</td>
<td style="text-align:right;">
0.39
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.98
</td>
<td style="text-align:right;">
0.31
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.84
</td>
<td style="text-align:right;">
0.23
</td>
</tr>
<tr>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.56
</td>
<td style="text-align:right;">
0.95
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
0.41
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.74
</td>
<td style="text-align:right;">
0.36
</td>
</tr>
<tr>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.04
</td>
<td style="text-align:right;">
0.75
</td>
<td style="text-align:right;">
0.71
</td>
<td style="text-align:right;">
0.92
</td>
<td style="text-align:right;">
0.91
</td>
<td style="text-align:right;">
0.81
</td>
<td style="text-align:right;">
0.11
</td>
<td style="text-align:right;">
0.12
</td>
</tr>
<tr>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.63
</td>
<td style="text-align:right;">
0.62
</td>
<td style="text-align:right;">
0.73
</td>
<td style="text-align:right;">
0.99
</td>
<td style="text-align:right;">
0.89
</td>
<td style="text-align:right;">
0.94
</td>
<td style="text-align:right;">
1.00
</td>
<td style="text-align:right;">
0.58
</td>
<td style="text-align:right;">
0.50
</td>
</tr>
<tr>
<td style="text-align:right;">
0.16
</td>
<td style="text-align:right;">
0.08
</td>
<td style="text-align:right;">
0.05
</td>
<td style="text-align:right;">
0.79
</td>
<td style="text-align:right;">
0.14
</td>
<td style="text-align:right;">
0.36
</td>
<td style="text-align:right;">
0.35
</td>
<td style="text-align:right;">
0.67
</td>
<td style="text-align:right;">
0.49
</td>
<td style="text-align:right;">
0.44
</td>
</tr>
</tbody>
</table>

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-liu" class="csl-entry">

Liu. 2009. “The Nonparanormal: Semiparametric Estimation of High
Dimensional Undirected Graphs.” *Journal of Machine Learning Research*.

</div>

<div id="ref-tmj" class="csl-entry">

Tomilina, Mazo, and Jaffrézic. 2024. “Mixed Copula-Based Models for
Semi-Parametric Network Inference: An Application to Multi-Omics Data.”

</div>

<div id="ref-ICGC" class="csl-entry">

Zhang, J., R. Bajari, and Andric D. 2019. “The International Cancer
Genome Consortium Data Portal.” *Nat Biotechnol.*
<https://doi.org/doi:10.1038/s41587-019-0055-9>.

</div>

</div>
