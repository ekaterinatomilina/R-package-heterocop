# Heterocop

Simulation and estimation in copula models for heterogeneous data.

## Description

Let $(X_1, ..., X_d)$ a vector of $d$ variables. One can express their joint cumulative distribution function (CDF) as a Gaussian copula, that is: $\begin{align*} F(x_1, ..., x_d) &=& C_{\Sigma}(F_1(x_1), ..., F_d(x_d))\\ &=&\Phi_{\Sigma}(\Phi^{-1}(F_1(x_1)), ..., \Phi^{-1}(F_d(x_d))) \end{align*}$

where $F_d$ denotes the marginal CDF of the variable $X_d$, $\Sigma$ the correlation matrix of the Gaussian copula, $\Phi$ the standard Gaussian CDF.

The **heterocop** package enables the user to simulate data distributed via the Gaussian copula as well as to estimate its correlation matrix.

### Estimation

The estimation is made via the pairwise maximum likelihood estimation method. Indeed, the pairwise correlation coefficient between variables $X_j$ and $X_k$ developed by Mazo et al. (2022) is used, which maximizes the pairwise-likelihood expressed as below:

$$
L(\theta)=\dfrac{1}{n}\sum_{a\in\mathcal{A}}\sum_{i=1}^n l_a(X_i;\theta)
$$

for $\theta \in \Theta$

The density $l_a$ takes the following form depending on the nature of the variables.

When both variables are continuous, blablabla

### Simulation

For the simulation, the generalized inverse method is used. It follows the algorithm below. blablabla
