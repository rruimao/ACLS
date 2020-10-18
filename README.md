# ACLS
Adaptive Capped Least Squares
## Description
This package includes two methods applied to minimize the adaptive resistant loss: randomized gradient descent method and gradient descent method with initials obtained from CPLEX.

Suppose we observe data vectors  $\{(x_i, y_i) \}_{i=1}^n$ that follow a linear model $y_i=x_i^{\text{T}} \beta^* +\epsilon_i, \ \ i=1,\ldots, n, $, where $y_i$ is a univariate response,  $x_i$ is a $d$-dimensional predictor, $\beta^*$ denotes the vector of regression coefficients, and $\epsilon_i$ is a random error. We propose the adpative resistant loss, $\ell(x)=x^2/2$ if $|x| \leq \tau$; $\tau^2/2,$ if $|x|>\tau$, where $\tau=\tau(n)>0$ is referred to as the resistance parameter. The proposed methods are applied to find $\beta$ that minimizes $\mathcal{L}(\beta)= n^{-1} \sum \ell(y_i-x _i^\text{T} \beta )$.

## Installtation
Install **ACLS** from GitHub:
``` R
install.packages("devtools")
library(devtools)
devtools::install_github("rruimao/ACLS")
library(ACLS)
``` 

Install **Rcplex** if want to use CPLEX to obtain the intial or the solution to the problem:

See the **INSTALL** file from https://cran.r-project.org/web/packages/Rcplex/index.html.


## Functions
There are three functions in this package:

-**GD**: Gradient descent method with user-defined initial.
-**RGD**: Randomized gradient descent method.
-**cplexcoef**: Get the coefficients for the CPLEX.


