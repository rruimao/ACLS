# ACLS
Adaptive Capped Least Squares
## Description
This package includes two methods applied to minimize the adaptive resistant loss: randomized gradient descent method and gradient descent method with initials obtained from CPLEX.

Suppose we observe data vectors  $\{(x_i, y_i) \}_{i=1}^n$ that follow a linear model $y_i=x_i^{\text{T}} \beta^* +\epsilon_i, \ \ i=1,\ldots, n, $, where $y_i$ is a univariate response,  $x_i$ is a $d$-dimensional predictor, $\beta^*$ denotes the vector of regression coefficients, and $\epsilon_i$ is a random error. We propose the adpative resistant loss, $\ell(x)=x^2/2$ if $|x| \leq \tau$; $\tau^2/2,$ if $|x|>\tau$, where $\tau=\tau(n)>0$ is referred to as the resistance parameter. The proposed methods are applied to find $\beta$ that minimizes $\mathcal{L}(\beta)= n^{-1} \sum \ell(y_i-x _i^\text{T} \beta )$.

## Installation
Install **ACLS** from GitHub:
``` R
install.packages("devtools")
library(devtools)
devtools::install_github("rruimao/ACLS/ACLS")
library(ACLS)
``` 

Install **Rcplex** if want to use CPLEX to obtain the intial or the solution to the problem:

See the **INSTALL** file from https://cran.r-project.org/web/packages/Rcplex/index.html.


## Functions
There are four functions in this package:

- **GD**: Gradient descent method with user-defined initial.
- **RGD**: Randomized gradient descent method.
- **picksamples**: Pick subsamples to obtain the initial.
- **cplexcoef**: Get the coefficients for the problem solved by CPLEX.

## Examples
We present examples of random generated $\{(x_i, y_i) \}_{i=1}^n$ with $n=50$, $d=5$ and $\beta^* =(0,3,4,1,2,0)^{\text{T}}$.

``` R
# n: sample size; d: dimensionality
n<-50
d<-5
a<-matrix(data = rnorm(1, 0, 1), nrow = n, ncol = d)
x_0<-matrix(1L,nrow=n,ncol=1)
X=cbind(x_0,a)
beta_true<-c(0,3,4,1,2,0)
eps<-matrix(rnorm(n));
#Genarate response Y using true coefficient beta_true
Y<-X %*% beta_true+eps
```

### Randomized gradient descent
We randomly generate 10 initials $\beta^0 \sim \text{Unif}(\mathbb{B}_2(\tau))$, where $\text{Unif}(\mathbb{B}_2(\tau))$ is a uniform distribution on the $\ell_2$-ball $\mathbb{B}_2(\tau)=\{x: \|x\|_2 \leq \tau \}$. This method finds the initial that provides the smallest adaptive resistant loss.

``` R
tau<-sqrt(n)/log(log(n))
iter<-10
beta.rgd<-RGD(X,Y,tau,iter)
```

### Gradient descent method with initials obtained from CPLEX
We first sample $30\%$ data and use CPLEX to find the optimal solution as the initial, then we apply gradient descent method with this initial.
``` R
n_ratio<-0.3
Sample<-picksamples(X,Y,n_ratio)
Coef_sp<-cplexcoef(Sample$X_sp,Sample$Y_sp,tau)
ans_sp<-Rcplex(Coef_sp$f,Coef_sp$A,Coef_sp$b,Coef_sp$Q,Coef_sp$lb,Coef_sp$ub,vtype=Coef_sp$vtype)
beta_0<-ans$xopt[(3*n+1):(3*n+d+1)]
beta.gdc<-GD(beta_0,tau,X,Y)
```

