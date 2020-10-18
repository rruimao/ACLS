# ACLS
Adaptive Capped Least Squares
## Description
This package includes two methods applied to minimize the adaptive resistant loss: randomized gradient descent method and gradient descent method with initials obtained from CPLEX.

Suppose we observe data vectors  $\{(x_i, y_i) \}_{i=1}^n$ that follow a linear model $y_i=x_i^{\text{T}} \beta^* +\epsilon_i, \ \ i=1,\ldots, n, $, where $y_i$ is a univariate response,  $x_i$ is a $d$-dimensional predictor, $\beta^*$ denotes the vector of regression coefficients, and $\epsilon_i$ is a random error. We propose the adpative resistant loss,$\ell(x)=x^2/2$ if $|x| \leq \tau$; $\tau^2/2,$ if $|x|>\tau,$

