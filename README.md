# ACLS
Adaptive Capped Least Squares
## Description
This package includes two methods applied to minimize the adaptive capped least squares loss: randomized gradient descent method and gradient descent method with initials obtained from CPLEX.

Suppose we observe data vectors  (x<sub>i</sub>i,y<sub>i</sub>) that follow a linear model $y_i=x_i^{\text{T}} \beta^* +\epsilon_i, \ \ i=1,\ldots, n, $, where $y_i$ is a univariate response,  $x_i$ is a $d$-dimensional predictor, $\beta^*$ denotes the vector of regression coefficients, and $\epsilon_i$ is a random error. We propose the adpative resistant loss, $\ell(x)=x^2/2$ if $|x| \leq \tau$; $\tau^2/2,$ if $|x|>\tau$, where $\tau=\tau(n)>0$ is referred to as the resistance parameter. The proposed methods are applied to find $\beta$ that minimizes $\mathcal{L}(\beta)= n^{-1} \sum \ell(y_i-x _i^\text{T} \beta )$.

## Installation
Install **ACLS** from GitHub:
``` R
install.packages("devtools")
library(devtools)
devtools::install_github("rruimao/ACLS/ACLS")
library(ACLS)
``` 

Install **Rcplex** if the user wants to use CPLEX to obtain the intial or the solution to the problem:

See the **INSTALL** file from https://cran.r-project.org/web/packages/Rcplex/index.html.


## Functions
There are four functions in this package:

- **GD**: Gradient descent method with user-defined initial.
- **RGD**: Randomized gradient descent method.
- **picksamples**: Pick subsamples to obtain the initial.
- **cplexcoef**: Get the coefficients for the problem solved by CPLEX.

## Examples
We present two examples: random generated data with $y$-outliers and random generated data with $x$-outliers and $y$-outliers. 



### First example: random generated data with $y$-outliers
we generate contaminated random errors $\epsilon_i$ from a mixture of normal distribution $0.9 \mathcal{N}(0,1)+0.1 \mathcal{N}(10,1)$ and $x_i$'s are independently and identically distributed (i.i.d.) as $\mathcal{N}(0,I_d)$ where $I_d$ is an identity matrix (using mvrnorm function from MASS package). We set $\beta^* =(0,3,4,1,2,0)^{\text{T}}$ to generate $y_i$.
``` R
# n: sample size; d: dimensionality
n<-50
d<-5
mu<-matrix(0L,nrow=d,ncol=1)
Sig<-diag(d)
install.packages("MASS")
library("MASS")
a<-mvrnorm(n, mu, Sig, tol = 1e-06, empirical = FALSE)
x_0<-matrix(1L,nrow=n,ncol=1)
X=cbind(x_0,a)
beta_true<-c(0,3,4,1,2,0)
components<-sample(1:2,prob=c(0.9,0.1),size=n,replace=TRUE)
mus<-c(0,10)
sds<-c(1,1)
eps<-matrix(rnorm(n,mean=mus[components],sd=sds[components]))
#Genarate response Y using true coefficient beta_true
Y<-X %*% beta_true+eps
```

We first apply our proposed methods, randomized gradient descent (denoted by ACLS) and gradient descent with intials obtained from CPLEX (denoted by ACLS-h), to the data. We then compare mean square errors (MSEs) of these two methods with MSEs of ordinary least squares method (OLS), Huber method with adaptive resistant parameter (denoted by AHR) and least trimmed squares method (LTS). Results obtained from the CPLEX on the whole dataset (denoted by ACLS-C) are treated as benchmark. 


### Randomized gradient descent
We randomly generate 10 initials $\beta^0 \sim \text{Unif}(\mathbb{B}_2(\tau))$, where $\text{Unif}(\mathbb{B}_2(\tau))$ is a uniform distribution on the $\ell_2$-ball $\mathbb{B}_2(\tau)=\{x: \|x\|_2 \leq \tau \}$. This method finds the initial that provides the smallest adaptive capped least squares loss.

``` R
tau<-sqrt(n)/log(log(n))
iter<-10
beta.rgd<-RGD(X,Y,tau,iter)
```

### Gradient descent method with initials obtained from CPLEX 
We first sample $30\%$ data and use CPLEX to find the optimal solution as the initial, then we apply gradient descent method with this initial.
``` R
library("Rcplex")
n_ratio<-0.3
n_sp<-round(n*n_ratio)
Sample<-picksamples(X,Y,n_ratio)
Coef_sp<-cplexcoef(Sample$X_sp,Sample$Y_sp,tau)
ans_sp<-Rcplex(Coef_sp$f,Coef_sp$A,Coef_sp$b,Coef_sp$Q,Coef_sp$lb,Coef_sp$ub,vtype=Coef_sp$vtype)
beta_0<-ans_sp$xopt[(3*n_sp+1):(3*n_sp+d+1)]
beta.gdc<-GD(beta_0,tau,X,Y)
```

### Other methods
``` R
install.packages("robustbase")
library("robustbase")
install.packages("robustreg")
library("robustreg")
beta.OLS<-solve(t(X)%*%X)%*%t(X)%*%Y
Z<-data.frame(X,Y)
colnames(Z)<-c("X1","X2","X3","X4","X5","X6","Y")
Z.Huber<-robustRegH(Y~X2+X3+X4+X5+X6,data=Z,tune=tau)
beta.Huber<-Z.Huber$coefficients
Z.LTS<-ltsReg(X[,-1],Y,intercept=TRUE,adjust=TRUE)
beta.LTS<-Z.LTS$coefficients
Coef<-cplexcoef(X,Y,tau)
ans<-Rcplex(Coef$f,Coef$A,Coef$b,Coef$Q,Coef$lb,Coef$ub,vtype=Coef$vtype)
beta.cplex<-ans$xopt[(3*n+1):(3*n+d+1)]
```

We summerize the MSEs of all methods in the following table.

|    |OLS | AHR |  LTS | ACLS | ACLS-h | ACLS-C |
| :---         |     :---:      |          ---: |          ---: |          ---: |          ---: |          ---: |
| MSE   | 2.8349     | 2.8349  |0.2096|0.2096|0.2096|0.2096|

### Second example: random generated data with $x$-outliers and $y$-outliers
we generate contaminated random errors $\epsilon_i$ from a mixture of normal distribution $0.9 \mathcal{N}(0,1)+0.1 \mathcal{N}(10,1)$ and $x_i$'s are independently and identically distributed (i.i.d.) as $\mathcal{N}(0,I_d)$ where $I_d$ is an identity matrix (using mvrnorm function from MASS package). We then add a random perturbation vector  $z_i \sim \mathcal{N}(10\times 1_{d-1},I_{d-1} )$
	to each covariate $x_i$ in the contaminated samples. We set $\beta^* =(0,3,4,1,2,0)^{\text{T}}$ and use uncontaminated $x_i$ to generate $y_i$.
	
``` R
# n: sample size; d: dimensionality
n<-50
d<-5
mu<-matrix(0L,nrow=d,ncol=1)
Sig<-diag(d)
install.packages("MASS")
library("MASS")
a<-mvrnorm(n, mu, Sig, tol = 1e-06, empirical = FALSE)
x_0<-matrix(1L,nrow=n,ncol=1)
X=cbind(x_0,a)
beta_true<-c(0,3,4,1,2,0)
X_index<-matrix(sample(1:n,round(0.1*n)))
R_new<-matrix(0L,nrow=n,ncol=d)
eps_2<-matrix(rnorm(n))
  for (i in 1:dim(X_index)[1])
  {
    j<-X_index[i]
    R_new[j,]<-rnorm(d, mean=10, sd=1)
    eps_2[j]<-rnorm(1,mean=10,sd=1)
  }
  R_new<-a+R_new
  X_new<-cbind(x_0,R_new)
#Genarate response Y using true coefficient beta_true
Y_2<-X %*% beta_true+eps_2
```	

We use the same code in the first example replacing X, Y with X_new, Y_2 respectively to get estimators for all methods. We also collect the MSEs in the following table.

|    |OLS | AHR |  LTS | ACLS | ACLS-h | ACLS-C |
| :---         |     :---:      |          ---: |          ---: |          ---: |          ---: |          ---: |
| MSE   | 16    | 16    |1|1|1|1|
